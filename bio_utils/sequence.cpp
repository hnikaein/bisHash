#include "sequence.h"
#include "../utils/logger.h"
#include <cstring>
#include <sys/stat.h>

using namespace std;

Sequence::Sequence(const char *const seq_str, const unsigned long seq_str_len, const char *name,
                   const char *quality_str, const int chr_num, const unsigned long chr_pos, const bool create_reverse) :
        seq_str(seq_str), size(seq_str_len), name(name), quality_str(quality_str), chr_num(chr_num),
        chr_pos(chr_pos) {
    if (create_reverse) {
        reverse_seq_str = new char[seq_str_len + 1];
        get_reversed(reverse_seq_str);
        delete_flag |= 0b1000;
    } else
        reverse_seq_str = nullptr;
}

Sequence::Sequence(Sequence &&sequence) noexcept: name(sequence.name), quality_str(sequence.quality_str),
                                                  seq_str(sequence.seq_str), size(sequence.size),
                                                  chr_num(sequence.chr_num), chr_pos(sequence.chr_pos),
                                                  reverse_seq_str(sequence.reverse_seq_str),
                                                  delete_flag(sequence.delete_flag) {
    sequence.name = sequence.quality_str = sequence.seq_str = nullptr; // TODO why?
}

Sequence::~Sequence() {
    if (delete_flag & 0b0100)
        free((void *) name);
    if (delete_flag & 0b0010)
        delete[] seq_str;
    if (delete_flag & 0b0001)
        delete[] quality_str;
    if (delete_flag & 0b1000)
        delete[] reverse_seq_str;
}

void Sequence::get_reversed(char *destination) const {
    auto *new_seq_str = destination - 1 + size;
    for (int i = 0; i < size; ++i)
        switch (seq_str[i]) {
            case 'A':
//            case 'a':
                new_seq_str[-i] = 'T';
                break;
            case 'C':
//            case 'c':
                new_seq_str[-i] = 'G';
                break;
            case 'G':
//            case 'g':
                new_seq_str[-i] = 'C';
                break;
            case 'T':
//            case 't':
                new_seq_str[-i] = 'A';
                break;
            default:
                new_seq_str[-i] = seq_str[i];
                break;
        }
    destination[size] = '\0';
}

[[maybe_unused]] Sequence Sequence::get_reversed() const {
    auto *new_seq_str = new char[size + 1];
    get_reversed(new_seq_str);
    Sequence sequence(new_seq_str, size, name, quality_str, -1 * chr_num, chr_pos);
    sequence.delete_flag = 0b010;
    return sequence;
}

vector<Sequence> Sequence::chunkenize_big_sequence(const vector<Sequence> &seqs, unsigned int chunk_size,
                                                   bool with_reverse, int chunk_diff) {
    if (chunk_diff == 0)
        chunk_diff = static_cast<int>(chunk_size / 2);
    vector<Sequence> chunks;
    for (const auto &seq: seqs) {
        for (int reverse = 0; reverse < (with_reverse ? 2 : 1); reverse++) {
            auto chr_num_with_reverse_applied = reverse ? -1 * seq.chr_num : seq.chr_num;
            auto seq_str_with_reverse_applied = reverse ? seq.reverse_seq_str : seq.seq_str;
            for (int e1 = 0, chunk_i = 0;; e1 += chunk_diff, chunk_i++) {
#ifdef _DEBUG
                auto chunk_name = strdup(
                        Logger::formatString("%s_%08d_%08d_%lu", seq.name, chr_num_with_reverse_applied,
                                             chunk_i, chunk_size).c_str());
#else
                char *chunk_name = nullptr;
#endif
                Sequence s1(seq_str_with_reverse_applied + e1,
                            min(static_cast<unsigned>(chunk_size), static_cast<unsigned>(seq.size - e1)),
                            chunk_name, nullptr, chr_num_with_reverse_applied, reverse ? seq.size - e1 + 1 : e1 + 1);
#ifdef _DEBUG
                s1.delete_flag = 0b100;
#endif
                chunks.push_back(std::move(s1));
                if (e1 + chunk_size >= seq.size)
                    break;
            }
        }
    }
    return chunks;
}

[[maybe_unused]] int Sequence::write_to_file(const char *file_name, const bool append, const bool force_write) const {
    if (!force_write) {
        struct stat st{};
        if (stat(file_name, &st) == 0)
            return 1;
    }
    auto file = fopen(file_name, append ? "a" : "w");
    char tempch = '>', endlch = '\n';
    fwrite(&tempch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(name, static_cast<size_t>(strlen(name)), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fwrite(seq_str, static_cast<size_t>(size), sizeof(char), file);
    fwrite(&endlch, static_cast<size_t>(1), sizeof(char), file);
    fclose(file);
    return 0;
}

