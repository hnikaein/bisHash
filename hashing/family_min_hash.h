#ifndef FAMILY_MIN_HASH_H
#define FAMILY_MIN_HASH_H

#include <map>
#include "hash_utils.h"

class FamilyMinHash {
public:
    FamilyMinHash(int decompose_letters, int kmer_length, int max_base_number, int p, int min_hash_count);

    std::vector<std::tuple<int, int, int>>
    get_sketch(const char *seq_str, int seq_size, int mode = SKETCH_MODE_NO_CONVERSION) const;

private:
    int decompose_letters;
    int kmer_length;
    int p;
    int min_hash_count;
    int seed;
    int seed_k;
    int seed_dl;
};

#endif //FAMILY_MIN_HASH_H