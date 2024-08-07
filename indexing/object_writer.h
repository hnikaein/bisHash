#include <vector>
#include <map>
#include <cstdio>

#ifndef VECTOR_WRITER_H
#define VECTOR_WRITER_H

[[maybe_unused]] void write_to_file(char const *file_name, const std::vector<int> *data, int size);

void write_data(FILE *file, std::vector<std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>>> **datas,
                int data_size, const int &chunks_count);

void write_data(const char *file_name, std::vector<std::pair<std::map<int, std::vector<int>>, std::map<int, std::vector<int>>>> **datas,
                int data_size, const int &chunks_count);

[[maybe_unused]] std::vector<int> *read_vectors_from_file(char const *file_name);

std::vector<std::pair<std::vector<int *>, std::vector<int *>>> *read_data(char const *file_name, int *&all_file, int &chunks_count);

std::vector<std::pair<std::vector<int *>, std::vector<int *>>> *read_data(int *all_file, int &chunks_count);

#endif //VECTOR_WRITER_H
