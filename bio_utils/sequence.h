#include <vector>
#include <string>


#ifndef SEQUENCE_H
#define SEQUENCE_H

class Sequence {
public:
    static std::vector<Sequence> chunkenize_big_sequence(const std::vector<Sequence> &seqs, unsigned int chunk_size,
                                                         bool with_reverse = true, int chunk_diff = 0);

    Sequence(const char *seq_str, unsigned long seq_str_len, const char *name = nullptr,
             const char *quality_str = nullptr, int chr_num = 0, unsigned long chr_pos = 0,
             bool create_reverse = false);

    Sequence(Sequence &&sequence) noexcept;

    void create_reverse();

    void get_reversed(char *destination = nullptr) const;

    [[maybe_unused]] [[nodiscard]] Sequence get_reversed() const;

    [[maybe_unused]] int write_to_file(const char *file_name, bool append = false, bool force_write = true) const;

    ~Sequence();

    const char *seq_str, *quality_str;
    char *reverse_seq_str;
    int chr_num;
    unsigned long size, chr_pos;
    int delete_flag = 0b000;
    const char *name;
};

#endif //SEQUENCE_H
