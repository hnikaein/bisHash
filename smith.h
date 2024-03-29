#ifndef BISHASH_SMITH_H
#define BISHASH_SMITH_H

#include "bio_utils/sequence.h"
#include "configs.h"

void
align_chunk_reads_phase2(const Sequence &chunk, const std::vector<std::pair<int, bool>> &chunk_reads, const std::vector<Sequence> &reads,
                         unsigned long chr_size, PenaltyConfig penalty_config,
                         std::vector<std::vector<std::pair<int, char *>>> &output_map, std::vector<int> &output_map_least_penalty);

#endif
