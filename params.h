#ifndef BISHASH_PARAMS_H
#define BISHASH_PARAMS_H

#include <climits>
#include "configs.h"

struct BisHashArgs {
    int log_level = 0, threads_count = THREADS_COUNT;
    int family_decompose_letters = FAMILY_DECOMPOSE_LETTERS, kmer_length = KMER_LENGTH, chunk_size = CHUNK_SIZE, chunk_overlap = CHUNK_OVERLAP;
    int from_read = 0, to_read = INT_MAX;
    double alt_matchs_ratio = ALT_RATIO_L1;
    bool read_index = DEFAULT_READ_INDEX, write_index = DEFAULT_WRITE_INDEX, only_create_index = false;
    char *ref_file_name = nullptr, *output_file_name = nullptr, *reads_file_name = nullptr;
    PenaltyConfig penalty_config{};
};

BisHashArgs read_args(int argc, char *argv[]);

#endif
