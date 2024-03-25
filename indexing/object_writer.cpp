#include "object_writer.h"
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <stdexcept>

using namespace std;

[[maybe_unused]] void write_to_file(const char *file_name, vector<vector<tuple<int, int, int>>> **datas, int size) {
    auto file = fopen(file_name, "wb");
    const int MYBUFSIZE = 10 * BUFSIZ;
    int buffer[MYBUFSIZE], write_size = 0;
    buffer[write_size++] = size;
    for (int i = 0; i < size; ++i) {
        auto data = datas[i];
        auto data_size = static_cast<int>(data->size());
        memcpy(buffer + (write_size++), &data_size, sizeof(int));
        for (auto &data_item: *data) {
            int data_item_size = static_cast<int>(data_item.size());
            if (write_size + data_item_size * 3 >= 8 * BUFSIZ) {
                fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
                write_size = 0;
            }
            memcpy(buffer + (write_size++), &data_item_size, sizeof(int));
            for (auto &data_item_itam: data_item) {
                buffer[write_size++] = get<0>(data_item_itam);
                buffer[write_size++] = get<1>(data_item_itam);
                buffer[write_size++] = get<2>(data_item_itam);
            }
        }
    }
    fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
    fclose(file);
}

[[maybe_unused]] vector<vector<tuple<int, int, int>>> *read_vector_of_maps_from_file(char const *file_name) {
    struct stat st{};
    auto file = fopen(file_name, "rb");
    if (!file || stat(file_name, &st) != 0)
        throw runtime_error(string("file ") + file_name + " is not present or not readable");
    auto all_file = new int[st.st_size / 4];
    fread(all_file, sizeof(int), static_cast<size_t>(st.st_size / 4), file);
    fclose(file);
    auto buffer_p = all_file;
    int data_size = *(buffer_p++);
    auto result = new vector<vector<tuple<int, int, int>>>[data_size];
    for (int i = 0; i < data_size; ++i) {
        int vec_size = *(buffer_p++);
        result[i].resize(vec_size);
        for (int j = 0; j < vec_size; ++j) {
            int map_size = *(buffer_p++);
            for (int k = 0; k < map_size; ++k) {
                int first = *(buffer_p++);
                int second = *(buffer_p++);
                int third = *(buffer_p++);
                result[i][j].emplace_back(first, second, third);
            }
        }
    }
    delete[] all_file;
    return result;
}


[[maybe_unused]] void write_to_file(const char *const file_name, const vector<int> *data, const int size) {
    auto file = fopen(file_name, "wb");
    const int MYBUFSIZE = 10 * BUFSIZ;
    int buffer[MYBUFSIZE], write_size = 0, last_data_size = -1, additional_zeros = 0;
    buffer[write_size++] = size;
    for (int i = 0; i < size; ++i) {
        auto data_size = static_cast<int>(data[i].size());
        if (last_data_size == 0) {
            if (last_data_size == data_size) {
                additional_zeros++;
                continue;
            } else {
                memcpy(buffer + (write_size++), &additional_zeros, sizeof(int));
                additional_zeros = 0;
            }
        }
        last_data_size = data_size;
        if (write_size + data_size >= MYBUFSIZE) {
            fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
            write_size = 0;
            if (data_size >= MYBUFSIZE) {
                fwrite(&data_size, static_cast<size_t>(1), sizeof(int), file);
                fwrite(&data[i][0], static_cast<size_t>(data_size), sizeof(int), file);
            }
        }
        if (data_size < MYBUFSIZE) {
            memcpy(buffer + (write_size++), &data_size, sizeof(int));
            if (data_size) {
                memcpy(buffer + write_size, &data[i][0], data_size * sizeof(int));
                write_size += data_size;
            }
        }
    }
    fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
    fclose(file);
}

[[maybe_unused]] vector<int> *read_vectors_from_file(char const *file_name) {
    struct stat st{};
    auto file = fopen(file_name, "rb");
    if (!file || stat(file_name, &st) != 0)
        throw runtime_error(string("file ") + file_name + " is not present or not readable");
    auto all_file = new int[st.st_size / 4];
    fread(all_file, sizeof(int), static_cast<size_t>(st.st_size / 4), file);
    fclose(file);
    auto buffer_p = all_file;
    int size = *(buffer_p++);
    auto result = new vector<int>[size];
    for (int i = 0; i < size; ++i) {
        int vec_size = *(buffer_p++);
        if (vec_size == 0) {
            int additional_zeros = *(buffer_p++);
            i += additional_zeros;
        } else {
            result[i].assign(buffer_p, buffer_p + vec_size);
            buffer_p += vec_size;
        }
    }
    delete[] all_file;
    return result;
}
