#include "hash_utils.h"

char inline get_nuc_value(char ch) {
    switch (ch) {
        case 'A':
//        case 'a':
            return NUC_INT_A;
        case 'C':
//        case 'c':
            return NUC_INT_C;
        case 'G':
//        case 'g':
            return NUC_INT_G;
        case 'T':
//        case 't':
            return NUC_INT_T;
        default:
            return NUC_INT_OTHER;
    }
}

char inline get_nuc_value_with_mode(char ch, int mode) {
    char ch_value = get_nuc_value(ch);
    if (mode == SKETCH_MODE_NO_CONVERSION)
        return ch_value;
    if (mode == SKETCH_MODE_WITH_C_T_CONVERSION) {
        if (ch_value == NUC_INT_C)
            return NUC_INT_T;
        return ch_value;
    }
    if (mode == SKETCH_MODE_WITH_G_A_CONVERSION) {
        if (ch_value == NUC_INT_G)
            return NUC_INT_A;
        return ch_value;
    }
    return ch_value;
}

int zigma_hash(const char *seq_str, const int start, const int end, const int prev_hash, const int mode, const int seed,
               const int seed_k, const int p) {
    int res = 0;
    if (prev_hash == 0)
        for (int i = start; i < end; ++i) {
            res *= seed;
            res += get_nuc_value_with_mode(seq_str[i], mode);
            res %= p;
        }
    else {
        res = prev_hash;
        res -= get_nuc_value_with_mode(seq_str[start - 1], mode) * seed_k;
        res %= p;
        res *= seed;
        res += get_nuc_value_with_mode(seq_str[end - 1], mode);
        res %= p;
    }
    return res;
}