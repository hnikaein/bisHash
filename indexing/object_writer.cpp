#include "object_writer.h"
#include <sys/stat.h>
#include <stdexcept>
#include <cstring>

using namespace std;

const int MYBUFSIZE = 10 * BUFSIZ;

inline void fill_dest_from_src(const int *src, const int &src_size, int &src_read, int *dest, const int &dest_size, int &dest_filled) {
    int remained = min(src_size - src_read, dest_size - dest_filled);
    memcpy(dest + dest_filled, src + src_read, remained * sizeof(int));
    dest_filled += remained;
    src_read += remained;
}

inline void check_buffer(FILE *file, const int *buffer, int &write_size) {
    if (write_size >= 9 * BUFSIZ) {
        fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
        write_size = 0;
    }
}

void write_data(const char *file_name, vector<pair<map<int, vector<int>>, map<int, vector<int>>>> **datas, int datas_size,
                const int &chunks_count) {
    auto file = fopen(file_name, "wb");
    write_data(file, datas, datas_size, chunks_count);
    fclose(file);
}

void write_data(FILE *file, vector<pair<map<int, vector<int>>, map<int, vector<int>>>> **datas, int datas_size,
                const int &chunks_count) {
    int buffer[MYBUFSIZE], write_size = 0;
    buffer[write_size++] = chunks_count;
    buffer[write_size++] = datas_size;
    for (int i = 0; i < datas_size; ++i) {
        auto data = datas[i];
        auto data_size = static_cast<int>(data->size());
        buffer[write_size++] = data_size;
        for (const auto &data_maps: *data) {
            for (const auto *data_map: {&data_maps.first, &data_maps.second}) {
                int data_map_size = static_cast<int>(data_map->size());
                buffer[write_size++] = data_map_size;
                check_buffer(file, buffer, write_size);
                for (const auto &data_item_pair: *data_map) {
                    buffer[write_size++] = data_item_pair.first;
                    buffer[write_size++] = static_cast<int>(data_item_pair.second.size());
                    check_buffer(file, buffer, write_size);
                    for (int j = 0; j < data_item_pair.second.size();) {
                        fill_dest_from_src(&data_item_pair.second[0], static_cast<int>(data_item_pair.second.size()), j, buffer, MYBUFSIZE,
                                           write_size);
                        check_buffer(file, buffer, write_size);
                    }
                }
            }
        }
    }
    fwrite(buffer, static_cast<size_t>(write_size), sizeof(int), file);
}

vector<pair<vector<int *>, vector<int *>>> *read_data(char const *file_name, int *&all_file, int &chunks_count) {
    auto file = fopen(file_name, "rb");
    struct stat st{};
    if (!file || stat(file_name, &st) != 0)
        throw runtime_error(string("file ") + file_name + " is not present or not readable");
    all_file = new int[st.st_size / 4];
    fread(all_file, sizeof(int), static_cast<size_t>(st.st_size / 4), file);
    fclose(file);
    auto result = read_data(all_file, chunks_count);
    return result;
}

vector<pair<vector<int *>, vector<int *>>> *read_data(int *all_file, int &chunks_count) {
    auto buffer_p = all_file;
    chunks_count = *(buffer_p++);
    int data_size = *(buffer_p++);
    auto result = new vector<pair<vector<int *>, vector<int *>>>[data_size];
    for (int i = 0; i < data_size; ++i) {
        int vec_size = *(buffer_p++);
        result[i].resize(vec_size);
        for (int j = 0; j < vec_size; ++j) {
            auto data_map = &result[i][j].first;
            for (int l = 0; l < 2; ++l) {
                if (l == 1)
                    data_map = &result[i][j].second;
                int map_size = *(buffer_p++);
                (*data_map).resize(map_size);
                for (int k = 0; k < map_size; ++k) {
                    (*data_map)[k] = buffer_p;
                    buffer_p += buffer_p[1] + 2;
                }
            }
        }
    }
    return result;
}


[[maybe_unused]] void write_to_file(const char *const file_name, const vector<int> *data, const int size) {
    auto file = fopen(file_name, "wb");
    int buffer[MYBUFSIZE], write_size = 0, last_data_size = -1, additional_zeros = 0;
    buffer[write_size++] = size;
    for (int i = 0; i < size; ++i) {
        auto data_size = static_cast<int>(data[i].size());
        if (last_data_size == 0) {
            if (last_data_size == data_size) {
                additional_zeros++;
                continue;
            } else {
                buffer[write_size++] = additional_zeros;
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
            buffer[write_size++] = data_size;
            if (data_size) {
                copy(&data[i][0], &data[i][data_size], buffer + write_size);
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
