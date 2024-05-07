/*
 *  Copyright (c) 2012-2017, Jyri J. Virkki
 *  All rights reserved.
 *
 *  This file is under BSD license. See LICENSE file.
 */

/*
 * Refer to bloom.h for documentation on the public interfaces.
 */

#include <assert.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "bloom.h"
#include "murmurhash2.h"

#define MAKESTRING(n) STRING(n)
#define STRING(n) #n

// extern void (*RedisModule_Free)(void *ptr);
// extern void * (*RedisModule_Calloc)(size_t nmemb, size_t size);

#define BLOOM_CALLOC calloc
#define BLOOM_FREE free

/*
#ifndef BLOOM_CALLOC
#define BLOOM_CALLOC calloc
#define BLOOM_FREE free
#endif
*/

#define MODE_READ 0
#define MODE_WRITE 1

#define LN2 (0.693147180559945)

inline static int test_bit_set_bit(unsigned char *buf, uint64_t x, int mode) {
    uint64_t byte = x >> 3;
    uint8_t mask = 1 << (x % 8);
    uint8_t c = buf[byte]; // expensive memory access

    if (c & mask) {
        return 1;
    } else {
        if (mode == MODE_WRITE) {
            buf[byte] = c | mask;
        }
        return 0;
    }
}

bloom_hashval bloom_calc_hash(const void *buffer, int len) {
    bloom_hashval rv;
    rv.a = murmurhash2(buffer, len, 0x9747b28c);
    rv.b = murmurhash2(buffer, len, rv.a);
    return rv;
}

bloom_hashval bloom_calc_hash64(const void *buffer, int len) {
    bloom_hashval rv;
    rv.a = MurmurHash64A_Bloom(buffer, len, 0xc6a4a7935bd1e995ULL);
    rv.b = MurmurHash64A_Bloom(buffer, len, rv.a);
    return rv;
}

// This function is defined as a macro because newer filters use a power of two
// for bit count, which is must faster to calculate. Older bloom filters don't
// use powers of two, so they are slower. Rather than calculating this inside
// the function itself, we provide two variants for this. The calling layer
// already knows which variant to call.
//
// modExp is the expression which will evaluate to the number of bits in the
// filter.
#define CHECK_ADD_FUNC(T, modExp)                                                                  \
    T i;                                                                                           \
    int found_unset = 0;                                                                           \
    const register T mod = modExp;                                                                 \
    for (i = 0; i < bloom->hashes; i++) {                                                          \
        T x = ((hashval.a + i * hashval.b)) % mod;                                                 \
        if (!test_bit_set_bit(bloom->bf, x, mode)) {                                               \
            if (mode == MODE_READ) {                                                               \
                return 0;                                                                          \
            }                                                                                      \
            found_unset = 1;                                                                       \
        }                                                                                          \
    }                                                                                              \
    if (mode == MODE_READ) {                                                                       \
        return 1;                                                                                  \
    }                                                                                              \
    return found_unset;

static int bloom_check_add32(struct bloom *bloom, bloom_hashval hashval, int mode) {
    CHECK_ADD_FUNC(uint32_t, (1 << bloom->n2));
}

static int bloom_check_add64(struct bloom *bloom, bloom_hashval hashval, int mode) {
    CHECK_ADD_FUNC(uint64_t, (1LLU << bloom->n2));
}

// This function is used for older bloom filters whose bit count was not
// 1 << X. This function is a bit slower, and isn't exposed in the API
// directly because it's deprecated
static int bloom_check_add_compat(struct bloom *bloom, bloom_hashval hashval, int mode) {
    CHECK_ADD_FUNC(uint64_t, bloom->bits)
}
// 这个`calc_bpe`函数是在布隆过滤器的上下文中用来计算“每个元素的比特数”（Bits Per Element），这是为了达到预定的错误率（false positive rate）。函数的输入是期望的错误率（`error`），输出是为了达到这个错误率所需要的每个元素的比特数。
//
// 布隆过滤器是一种空间效率很高的数据结构，用于测试一个元素是否在一个集合中。它允许有一定的错误率（false positives），但不会产生错误的否定（false negatives）。换句话说，布隆过滤器可能会告诉你一个元素在集合中，即使它不在；但如果它告诉你一个元素不在集合中，那么这个元素确实不在集合中。
//
// 函数中的`denom`（分母）是$\ln(2)^2$的值。这个值来自于布隆过滤器的最优比特数（每个元素）的计算公式，该公式旨在最小化在给定错误率下的空间使用。具体来说，最优的比特数（每个元素）可以通过以下公式计算：
//
// $$
// \text{BPE} = \frac{-\ln(\text{error rate})}{(\ln(2))^2}
// $$
//
// 这里：
// - $\ln$是自然对数。
// - `error rate`是期望的错误率（false positive rate）。
//
// 函数`calc_bpe`实现了上述公式。它首先计算错误率的自然对数（`log(error)`），然后用这个值除以$\ln(2)^2$，最后取结果的绝对值（以确保结果为正，因为错误率（`error`）应该是小于1的正数，所以其对数是负的）。
// 
// 这个计算结果告诉我们，为了达到指定的错误率，我们需要为布隆过滤器中的每个元素分配多少比特。这个值是设计布隆过滤器时的关键参数之一，它影响着布隆过滤器的空间效率和准确性。
// 假设计算结果为 bpe = 9.6。这意味着为了达到1%的误报率，
// 我们需要为布隆过滤器中的每个元素分配大约10位。由于位数不能为小数，我们通常会取整数位，即在实际应用中我们会选择10位或向上取整到接近的整数值
// 根据计算出的 bpe 值，如果我们预期要存储1000个元素，那么布隆过滤器的总位数应该是 1000 * 10 = 10000位
static double calc_bpe(double error) {
    static const double denom = 0.480453013918201; // ln(2)^2
    double num = log(error);

    double bpe = -(num / denom);
    if (bpe < 0) {
        bpe = -bpe;
    }
    return bpe;
}

int bloom_init(struct bloom *bloom, uint64_t entries, double error, unsigned options) {
    if (entries < 1 || error <= 0 || error >= 1.0) {
        return 1;
    }
    // 容错率
    bloom->error = error;
    // 位数
    bloom->bits = 0;
    // 预估的元素个数
    bloom->entries = entries;
    // 每个元素占用的位数
    bloom->bpe = calc_bpe(error);

    uint64_t bits;

    if (options & BLOOM_OPT_ENTS_IS_BITS) {
        // Size is determined by the number of bits
        if (/* entries == 0 || */ entries > 64) {
            return 1;
        }

        bloom->n2 = entries;
        bits = 1LLU << bloom->n2;
        bloom->entries = bits / bloom->bpe;

    } else if (options & BLOOM_OPT_NOROUND) {
        // Don't perform any rounding. Conserve memory instead
        bits = (uint64_t)(entries * bloom->bpe);

        // Guard against very small 'bpe'. Have at least one bit in the filter.
        if (bits == 0) {
            bits = 1;
        }
        bloom->bits = bits;
        bloom->n2 = 0;

    } else {
        // 数组的位数
        double bn2 = logb(entries * bloom->bpe);
        if (bn2 > 63 || bn2 == INFINITY) {
            return 1;
        }
        bloom->n2 = bn2 + 1;
        bits = 1LLU << bloom->n2;

        // Determine the number of extra bits available for more items. We rounded
        // up the number of bits to the next-highest power of two. This means we
        // might have up to 2x the bits available to us.
        size_t bitDiff = bits - (entries * bloom->bpe);
        // The number of additional items we can store is the extra number of bits
        // divided by bits-per-element
        size_t itemDiff = bitDiff / bloom->bpe;
        bloom->entries += itemDiff;
    }
    // 将位数转换为字节数
    if (bits % 64) {
        // 如果bits不是64的整倍数，意味着最后一个64位块不会完全用满，这时通过((bits / 64) + 1) * 8计算出需要额外的一个64位块来存储剩余的位，并且将总字节大小向上调整到完整的64位块
        bloom->bytes = ((bits / 64) + 1) * 8;
    } else {
        bloom->bytes = bits / 8;
    }
    // 更新实际的位数
    bloom->bits = bloom->bytes * 8;

    bloom->force64 = (options & BLOOM_OPT_FORCE64);
    bloom->hashes = (int)ceil(LN2 * bloom->bpe); // ln(2)
    bloom->bf = (unsigned char *)BLOOM_CALLOC(bloom->bytes, sizeof(unsigned char));
    if (bloom->bf == NULL) {
        return 1;
    }

    return 0;
}

int bloom_check_h(const struct bloom *bloom, bloom_hashval hash) {
    if (bloom->n2 > 0) {
        if (bloom->force64 || bloom->n2 > 31) {
            return bloom_check_add64((void *)bloom, hash, MODE_READ);
        } else {
            return bloom_check_add32((void *)bloom, hash, MODE_READ);
        }
    } else {
        return bloom_check_add_compat((void *)bloom, hash, MODE_READ);
    }
}

int bloom_check(const struct bloom *bloom, const void *buffer, int len) {
    return bloom_check_h(bloom, bloom_calc_hash(buffer, len));
}

int bloom_add_h(struct bloom *bloom, bloom_hashval hash) {
    if (bloom->n2 > 0) {
        if (bloom->force64 || bloom->n2 > 31) {
            return !bloom_check_add64(bloom, hash, MODE_WRITE);
        } else {
            return !bloom_check_add32(bloom, hash, MODE_WRITE);
        }
    } else {
        return !bloom_check_add_compat(bloom, hash, MODE_WRITE);
    }
}

int bloom_add(struct bloom *bloom, const void *buffer, int len) {
    return bloom_add_h(bloom, bloom_calc_hash(buffer, len));
}

void bloom_free(struct bloom *bloom) { BLOOM_FREE(bloom->bf); }

const char *bloom_version() { return MAKESTRING(BLOOM_VERSION); }

// Returns 0 on success
int bloom_validate_integrity(struct bloom *bloom) {
    if (bloom->error <= 0 || bloom->error >= 1.0 ||
        (bloom->n2 != 0 && bloom->bits < (1ULL << bloom->n2)) || bloom->bits == 0 ||
        bloom->bits != bloom->bytes * 8 || bloom->hashes != (int)ceil(LN2 * bloom->bpe)) {
        return 1;
    }

    return 0;
}
