#ifndef BISHASH_PARAMS_H
#define BISHASH_PARAMS_H

#include "configs.h"

struct BisHashArgs {
    bool read_index = DEFAULT_READ_INDEX, write_index = DEFAULT_WRITE_INDEX;
    double alt_matchs_ratio = ALT_RATIO_L1;
    int log_level, threads_count = THREADS_COUNT;
    int family_decompose_letters = FAMILY_DECOMPOSE_LETTERS, kmer_length = KMER_LENGTH, chunk_size = CHUNK_SIZE, chunk_overlap = CHUNK_OVERLAP;
    int from_read = 0, to_read = 0;
    char *ref_file_name, *output_file_name, *reads_file_name;
    PenaltyConfig penalty_config{};
};

BisHashArgs read_args(int argc, char *argv[]);

#endif
