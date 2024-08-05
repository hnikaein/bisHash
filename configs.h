#ifndef CONFIGS_H
#define CONFIGS_H

#define THREADS_COUNT                       1
#define MAX_READ_LEN                        2000
#define ALT_MATCHS                          10
#define ALT_RATIO_L1                        0.8
#define ALT_RATIO_L2                        0.98

#define CHUNK_SIZE                          7000
#define CHUNK_OVERLAP                       MAX_READ_LEN

#define MAX_BASENUMBER                      5
#define FAMILY_DECOMPOSE_LETTERS            5
#define KMER_LENGTH                         20
#define BIG_PRIME_NUMBER                    16777215 // 49979693 16777215 14868587 16777259
#define MIN_HASH_COUNT                      2

#define DEFAULT_PENALTY_MATCH               5
#define DEFAULT_PENALTY_MISMATCH            5
#define DEFAULT_PENALTY_GAP_OPEN            5
#define DEFAULT_PENALTY_GAP_EXTEND          3
#define DEFAULT_READ_INDEX                  true
#define DEFAULT_WRITE_INDEX                 true


struct PenaltyConfig {
    int match_score = DEFAULT_PENALTY_MATCH;
    int mismath_penalty = DEFAULT_PENALTY_MISMATCH;
    int gap_open_penalty = DEFAULT_PENALTY_GAP_OPEN;
    int gap_extend_penalty = DEFAULT_PENALTY_GAP_EXTEND;
};
#endif //CONFIGS_H
