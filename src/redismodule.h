#ifndef REDIS_MODULE_H
#define REDIS_MODULE_H

#include <stdlib.h>
#include <sds.h>
#define RedisModule_Alloc malloc
#define RedisModule_Free free
#define RedisModule_Realloc realloc
#define RedisModule_Calloc calloc

#define REDISMODULE_OK 0
#define REDISMODULE_ERR 1

#define RedisModuleString sds
#endif /* REDIS_MODULE_H */
