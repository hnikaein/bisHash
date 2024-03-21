#ifndef HASH_UTILS_H
#define HASH_UTILS_H


#define SKETCH_MODE_NO_CONVERSION           0
#define SKETCH_MODE_WITH_C_T_CONVERSION     1
#define SKETCH_MODE_WITH_G_A_CONVERSION     2

#define NUC_INT_A                           0
#define NUC_INT_C                           1
#define NUC_INT_G                           2
#define NUC_INT_T                           3
#define NUC_INT_OTHER                       4


int zigma_hash(const char *seq_str, int start, int end, int prev_hash, int mode, int seed, int seed_k, int p);

#endif //HASH_UTILS_H
