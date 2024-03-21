#include "sequence.h"
#include <vector>
#include <tuple>

#ifndef BIO_READER_H
#define BIO_READER_H

enum [[maybe_unused]] FileType {
    FASTA, FASTQ, SAM
};

[[maybe_unused]] std::tuple<std::vector<Sequence>, char *> read_sequences_from_file(const char *file_name);

std::tuple<std::vector<Sequence>, char *>
read_sequences_from_file(const char *file_name, const FileType &file_type, bool create_reverse = false);

std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>, char *> read_fasta(const char *file_name);

std::tuple<std::vector<char *>, std::vector<char *>, std::vector<int>, std::vector<char *>, char *>
read_fastq(const char *file_name);

#endif //BIO_READER_H
