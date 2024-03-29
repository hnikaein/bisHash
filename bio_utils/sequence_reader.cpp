#include "sequence_reader.h"
#include <sys/stat.h>
#include <cstring>
#include <stdexcept>

using namespace std;

[[maybe_unused]] tuple<vector<Sequence>, char *> read_sequences_from_file(const char *file_name) {
    if (tolower(file_name[strlen(file_name) - 1]) == 'a')
        return read_sequences_from_file(file_name, FASTA);
    else
        return read_sequences_from_file(file_name, FASTQ);
}

inline char upper_nuc(char nuc) {
    switch (nuc) {
        case 'a':
            return 'A';
        case 'c':
            return 'C';
        case 'g':
            return 'G';
        case 't':
            return 'T';
        default:
            return nuc;
    }
}

tuple<vector<Sequence>, char *>
read_sequences_from_file(const char *file_name, const FileType &file_type, const bool create_reverse) {
    vector<char *> names, seqs, quality;
    vector<int> lens;
    char *should_be_deleted = nullptr;
    if (file_type == FASTA)
        tie(names, seqs, lens, should_be_deleted) = read_fasta(file_name);
    else if (file_type == FASTQ)
        tie(names, seqs, lens, quality, should_be_deleted) = read_fastq(file_name);
    auto l = vector<Sequence>();
    l.reserve(names.size());
    if (quality.empty())
        for (int i = 0; i < names.size(); i++) {
            for (int j = 0; j < lens[i]; ++j)
                seqs[i][j] = upper_nuc(seqs[i][j]);
            l.emplace_back(seqs[i], lens[i], names[i], nullptr, i + 1, 1, create_reverse);
        }
    else
        for (int i = 0; i < names.size(); i++) {
            for (int j = 0; j < lens[i]; ++j)
                seqs[i][j] = upper_nuc(seqs[i][j]);
            l.emplace_back(seqs[i], lens[i], names[i], quality[i], i + 1, 1, create_reverse);
        }
    return make_tuple(std::move(l), should_be_deleted);
}


tuple<vector<char *>, vector<char *>, vector<int>, char *> read_fasta(const char *const file_name) {
    vector<char *> seqs, names;
    vector<int> lens;
    struct stat st{};
    auto file = fopen(file_name, "r");
    if (!file || stat(file_name, &st) != 0)
        throw runtime_error(string("file ") + file_name + " is not present or not readable");

    auto result = new char[st.st_size];
    bool name_state = true;
    long result_i = 0, last_result_i = 0;
    char buffer[BUFSIZ], last_ch, ch = '\0';
    size_t read_len;

    fread(buffer, 1, 1, file);
    names.push_back(result);

    while ((read_len = fread(buffer, 1, BUFSIZ, file))) {
        for (int i = 0; i < read_len; i++) {
            last_ch = ch;
            ch = buffer[i];
            if (ch == '>') {
                result[result_i++] = '\0';
                name_state = true;
                names.push_back(result + result_i);
                lens.push_back(static_cast<int>(result_i - last_result_i - 1));
                continue;
            } else if (ch == '\n' || ch == '\r') {
                if (last_ch != '\n' && last_ch != '\r' && name_state) {
                    result[result_i++] = '\0';
                    name_state = false;
                    seqs.push_back(result + result_i);
                    last_result_i = result_i;
                }
                continue;
            }
            result[result_i++] = ch;
        }
    }
    result[result_i] = '\0';
    lens.push_back(static_cast<int>(result_i - last_result_i));
    fclose(file);
    return make_tuple(std::move(names), std::move(seqs), std::move(lens), result);
}

tuple<vector<char *>, vector<char *>, vector<int>, vector<char *>, char *> read_fastq(const char *const file_name) {
    vector<char *> seqs, names, quality;
    vector<int> lens;
    auto file = fopen(file_name, "r");
    if (!file)
        return make_tuple(std::move(names), std::move(seqs), std::move(lens), std::move(quality), nullptr);
    struct stat st{};
    if (stat(file_name, &st) != 0)
        throw runtime_error(string("file ") + file_name + " is not present or not readable");

    auto result = new char[st.st_size];

    int state = 0;
    int result_i = 0, last_result_i = 0;
    char buffer[BUFSIZ], last_ch, ch = '\0';
    size_t read_len;
    names.push_back(result + 1);


    while ((read_len = fread(buffer, 1, BUFSIZ, file))) {
        for (int i = 0; i < read_len; i++) {
            last_ch = ch;
            ch = buffer[i];
            if (ch == '\n' || ch == '\r') {
                if (last_ch == '\n' || last_ch == '\r')
                    continue;
                if (state == 0 || state == 1 || state == 3)
                    result[result_i++] = '\0';
                state = (state + 1) % 4;
                if (state == 0)
                    names.push_back(result + result_i + 1);
                if (state == 1) {
                    seqs.push_back(result + result_i);
                    last_result_i = result_i;
                }
                if (state == 2)
                    lens.push_back(result_i - last_result_i - 1);
                if (state == 3)
                    quality.push_back(result + result_i);
            } else if (state != 2)
                result[result_i++] = ch;
        }
    }
    fclose(file);
    result[result_i] = '\0';
    if (names.size() > lens.size()) // in case of a trailing newline
        names.pop_back();
    return make_tuple(std::move(names), std::move(seqs), std::move(lens), std::move(quality), result);
}
