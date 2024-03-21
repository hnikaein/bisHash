#ifndef CONFIGS_H
#define CONFIGS_H

extern int match_score, mismath_penalty, gap_open_penalty, gap_extend_penalty;

#define THREADS_COUNT                       1
#define MAX_READ_LEN                        2000
#define ALT_MATCHS                          30
#define ALT_RATIO_L1                        0.8
#define ALT_RATIO_L2                        0.98

#define CHUNK_SIZE                          7000
#define CHUNK_OVERLAP                       MAX_READ_LEN

#define MAX_BASENUMBER                      8
#define FAMILY_DECOMPOSE_LETTERS            5
#define KMER_LENGTH                         20
#define BIG_PRIME_NUMBER                    49979693 // 16777215 14868587 16777259
#define MIN_HASH_COUNT                      2

#define DEFAULT_PENALTY_MATCH               5
#define DEFAULT_PENALTY_MISMATCH            5
#define DEFAULT_PENALTY_GAP_OPEN            5
#define DEFAULT_PENALTY_GAP_EXTEND          3
#define DEFAULT_READ_INDEX                  true
#define DEFAULT_WRITE_INDEX                 true

//#define array_len(x)                ((sizeof(x)) / (sizeof((x)[0])))

#endif //CONFIGS_H
