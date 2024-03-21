#include <vector>
#include "family_min_hash.h"

using namespace std;

FamilyMinHash::FamilyMinHash(int decompose_letters, int kmer_length, int max_base_number, int p, int min_hash_count) :
        decompose_letters(decompose_letters), kmer_length(kmer_length), seed(max_base_number), p(p),
        min_hash_count(min_hash_count) {
    seed_k = 1;
    for (int i = 0; i < kmer_length - decompose_letters - 1; ++i) {
        seed_k *= seed;
        seed_k %= p;
    }
    seed_dl = 1;
    for (int i = 0; i < decompose_letters - 1; ++i) {
        seed_dl *= seed;
        seed_dl %= p;
    }
}

vector<tuple<int, int, int>> FamilyMinHash::get_sketch(const char *seq_str, const int seq_size, int mode) const {
    map<int, tuple<int, int, int>> sketch;
    int hash_value = 0;
    int hash_value_no = 0;
    int decompose_letters_hash = 0;

    for (int i = 0; i < seq_size - kmer_length + 1; ++i) {
        decompose_letters_hash = zigma_hash(seq_str, i, i + decompose_letters,
                                            decompose_letters_hash,
                                            mode,
                                            seed,
                                            seed_dl,
                                            p
        );

        hash_value = zigma_hash(seq_str, i + decompose_letters, i + kmer_length,
                                hash_value,
                                mode,
                                seed,
                                seed_k,
                                p
        );
        hash_value_no = zigma_hash(seq_str, i + decompose_letters, i + kmer_length,
                                   hash_value_no,
                                   SKETCH_MODE_NO_CONVERSION,
                                   seed,
                                   seed_k,
                                   p
        );
        if (!sketch.count(decompose_letters_hash))
            sketch[decompose_letters_hash] = make_tuple(1, hash_value, hash_value_no);
        else {
            auto prev_sketch = sketch[decompose_letters_hash];
            sketch[decompose_letters_hash] = make_tuple(get<0>(prev_sketch) + 1, min(hash_value, get<1>(prev_sketch)),
                                                        min(hash_value_no, get<2>(prev_sketch)));
        }
    }
    vector<tuple<int, int, int>> nsketch;
    for (const auto &sketch_item: sketch) {
        decompose_letters_hash = sketch_item.first;
        int hash_count;
        tie(hash_count, hash_value, hash_value_no) = sketch_item.second;
        if (hash_count >= min_hash_count)
            nsketch.emplace_back(decompose_letters_hash, hash_value, hash_value_no);
    }
    return nsketch;
}