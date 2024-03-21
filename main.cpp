#include <cstring>
#include <algorithm>
#include <getopt.h>
#include <map>
#include <cmath>
#include "configs.h"
#include "utils/multiproc.h"
#include "utils/logger.h"
#include "utils/time_profile.h"
#include "bio_utils/bio_reader.h"
#include "hashing/family_min_hash.h"
#include "indexing/object_writer.h"

using namespace std;

extern Logger *logger;

FamilyMinHash *family_min_hash;

vector<Sequence> ref_genome, chunks, reads;
vector<vector<tuple<int, int, int>>> chunks_sketchs_CT, chunks_sketchs_GA;
vector<vector<int>> reads_chunks;
vector<vector<pair<int, bool>>> chunks_reads;
vector<vector<pair<int, char *>>> output_map;
vector<int> output_map_least_penalty;

pthread_mutex_t output_map_inserting;
FILE *output_file;
string config_str;

bool read_index = DEFAULT_READ_INDEX, write_index = DEFAULT_WRITE_INDEX;
double alt_matchs_ratio = ALT_RATIO_L1;
int log_level = Logger::INFO, threads_count = THREADS_COUNT;
int match_score = DEFAULT_PENALTY_MATCH, mismath_penalty = DEFAULT_PENALTY_MISMATCH, gap_open_penalty = DEFAULT_PENALTY_GAP_OPEN, gap_extend_penalty = DEFAULT_PENALTY_GAP_EXTEND;
int family_decompose_letters = FAMILY_DECOMPOSE_LETTERS, kmer_length = KMER_LENGTH, chunk_size = CHUNK_SIZE, chunk_overlap = CHUNK_OVERLAP;
int from_read = 0, to_read = 0;
int log_level_power;
int chunks_size;
char *ref_file_name, *output_file_name, *reads_file_name;
char *ref_genome_should_be_deleted, *reads_should_be_deleted;


void read_args(int argc, char *argv[]) {
    static struct option long_options[] =
            {
                    {"threads",                  required_argument, nullptr, 't'},
                    {"ref",                      required_argument, nullptr, 'r'},
                    {"query",                    required_argument, nullptr, 'q'},
                    {"output",                   required_argument, nullptr, 'o'},
                    {"log",                      required_argument, nullptr, 'l'},
                    {"alt-match-ratio",          required_argument, nullptr, 'm'},
                    {"no-read-index",            no_argument,       nullptr, 'n'},
                    {"no-write-index",           no_argument,       nullptr, 'w'},
                    {"match-score",              required_argument, nullptr, 'A'},
                    {"mismatch-penalty",         required_argument, nullptr, 'B'},
                    {"gap-open-penalty",         required_argument, nullptr, 'O'},
                    {"gap-extend-penalty",       required_argument, nullptr, 'E'},
                    {"from-read",                required_argument, nullptr, 'F'},
                    {"to-read",                  required_argument, nullptr, 'T'},
                    {"family-decompose-letters", required_argument, nullptr, 'D'},
                    {"kmer-length",              required_argument, nullptr, 'K'},
                    {"chunk-size",               required_argument, nullptr, 'S'},
                    {"chunk-overlap",            required_argument, nullptr, 'V'},
            };

    int option_index = 0, c;
    while ((c = getopt_long(argc, argv, "t:r:q:o:m:l:nwA:B:O:E:F:T:D:K:S:V:", long_options, &option_index)) >= 0)
        switch (c) {
            case 't':
                threads_count = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'r':
                ref_file_name = strdup(optarg);
                break;
            case 'q':
                reads_file_name = strdup(optarg);
                break;
            case 'o':
                output_file_name = strdup(optarg);
                break;
            case 'm':
                alt_matchs_ratio = strtod(optarg, nullptr);
                break;
            case 'l':
                log_level = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'n':
                read_index = false;
                break;
            case 'w':
                write_index = false;
                break;
            case 'A':
                match_score = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'B':
                mismath_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'O':
                gap_open_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'E':
                gap_extend_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'F':
                from_read = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'T':
                to_read = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'D':
                family_decompose_letters = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'K':
                kmer_length = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'S':
                chunk_size = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'V':
                chunk_overlap = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            default:
                break;
        }
    family_min_hash = new FamilyMinHash(family_decompose_letters, kmer_length, MAX_BASENUMBER, BIG_PRIME_NUMBER,
                                        MIN_HASH_COUNT);
    log_level_power = max(static_cast<int>(pow(10, 8 - log_level)), 1);
    logger = new Logger(log_level);
    optind = 1;
}

template<typename T>
void remove_vector(vector<T> &vec) {
    vector<T>().swap(vec);
}

// Calculate the chunk sketch and add the chunk id to hashes that exist in chunk sketch
int make_chunk_sketch(const int chunk_i) {
    chunks_sketchs_CT[chunk_i] = family_min_hash->get_sketch(chunks[chunk_i].seq_str,
                                                             static_cast<int>(chunks[chunk_i].size),
                                                             SKETCH_MODE_WITH_C_T_CONVERSION);
    chunks_sketchs_GA[chunk_i] = family_min_hash->get_sketch(chunks[chunk_i].seq_str,
                                                             static_cast<int>(chunks[chunk_i].size),
                                                             SKETCH_MODE_WITH_G_A_CONVERSION);
    return 1;
}

void read_chunks() {
    tie(ref_genome, ref_genome_should_be_deleted) = read_sequences_from_file(ref_file_name, FASTA, true);
    chunks = Sequence::chunkenize_big_sequence(ref_genome, chunk_size, true, chunk_size - chunk_overlap);
}

void prepare_ref_sketch() {
    auto index_file_name = config_str + ".cbhs";
    add_time();
    if (read_index)
        try {
            auto indexes = read_vector_of_maps_from_file(index_file_name.c_str());
            chunks_sketchs_CT = std::move(indexes[0]);
            chunks_sketchs_GA = std::move(indexes[1]);
            chunks_size = static_cast<int>(chunks_sketchs_CT.size());
            logger->info("index read completed");
        } catch (...) {
            read_index = false;
        }
    if (!read_index) {
        read_chunks();
        chunks_size = static_cast<int>(chunks.size());
        chunks_sketchs_CT.resize(chunks_size);
        chunks_sketchs_GA.resize(chunks_size);
        multiproc(threads_count, make_chunk_sketch, chunks_size);
    }
    add_time();
    if (!read_index && write_index) {
        vector<vector<tuple<int, int, int>>> *indexes[] = {&chunks_sketchs_CT, &chunks_sketchs_GA};
        write_to_file(index_file_name.c_str(), indexes, 2);
    }
    add_time();
    logger->info("ref sketchs prepared: %sms: %d records", get_times_str(true), chunks_size);
}

int find_read_chunks(const int read_i) {
    auto read = &reads[read_i];

    auto read_sketch_ct = family_min_hash->get_sketch(read->seq_str, static_cast<int>(read->size),
                                                      SKETCH_MODE_WITH_C_T_CONVERSION);
    auto read_sketch_ga = family_min_hash->get_sketch(read->seq_str, static_cast<int>(read->size),
                                                      SKETCH_MODE_WITH_G_A_CONVERSION);
    vector<int> max_sim_i(ALT_MATCHS);
    vector<int> max_sim(ALT_MATCHS);
    int last_inserted_j = 0;
    auto read_sketch = &read_sketch_ct;
    auto chunks_sketchs = &chunks_sketchs_CT;
    auto chunks_sketchs_i = -1;
    for (int sketch_type = 0; sketch_type < 2; ++sketch_type) {
        if (sketch_type == 1) {
            read_sketch = &read_sketch_ga;
            chunks_sketchs = &chunks_sketchs_GA;
        }
        for (const auto &chunk_sketch: *chunks_sketchs) {
            chunks_sketchs_i++;
            if (max_sim.size() > 3 * ALT_MATCHS) {
                max_sim.resize(2 * ALT_MATCHS);
                max_sim.reserve(4 * ALT_MATCHS);
                max_sim_i.resize(2 * ALT_MATCHS);
                max_sim_i.reserve(4 * ALT_MATCHS);
            }
            int sketch_sim = 0;
            int chunk_sketch_i = 0;
            for (const auto &read_sketch_i: *read_sketch)
                while (chunk_sketch_i < chunk_sketch.size()) {
                    if (get<0>(chunk_sketch[chunk_sketch_i]) > get<0>(read_sketch_i))
                        break;
                    if (get<0>(chunk_sketch[chunk_sketch_i]) == get<0>(read_sketch_i))
                        if (get<1>(chunk_sketch[chunk_sketch_i]) == get<2>(read_sketch_i) ||
                            get<2>(chunk_sketch[chunk_sketch_i]) == get<2>(read_sketch_i))
                            sketch_sim += 1;
                    chunk_sketch_i++;
                }
            if (sketch_sim <= max_sim[ALT_MATCHS - 1])
                continue;
            if (max_sim_i[last_inserted_j] == chunks_sketchs_i - 1) {
                if (max_sim[last_inserted_j] < sketch_sim) {
                    max_sim_i.erase(max_sim_i.begin() + last_inserted_j);
                    max_sim.erase(max_sim.begin() + last_inserted_j);
                } else {
                    if (max_sim[last_inserted_j] == sketch_sim)
                        max_sim_i[last_inserted_j] = chunks_sketchs_i;
                    continue;
                }
            }
            for (int j = 0; j < ALT_MATCHS; ++j)
                if (sketch_sim > max_sim[j]) {
                    max_sim_i.insert(max_sim_i.begin() + j, chunks_sketchs_i);
                    max_sim.insert(max_sim.begin() + j, sketch_sim);
                    last_inserted_j = j;
                    break;
                }
        }
    }
    int idx = 0;
    if (max_sim[0] != 0)
        for (; idx < max_sim.size(); ++idx)
            if (alt_matchs_ratio * max_sim[0] > max_sim[idx])
                break;
    max_sim_i.resize(min(idx, ALT_MATCHS));
    reads_chunks[read_i] = std::move(max_sim_i);
    return 1;
}

inline int write_cigar_block(int last_cigar_char_len, char last_cigar_char, char *cigar) {
    if (last_cigar_char_len == 1)
        return sprintf(cigar, "%c", last_cigar_char);
    else
        return sprintf(cigar, "%d%c", last_cigar_char_len, last_cigar_char);
}

void align_chunk_reads_phase2(int chunk_i) {
    char reverse_cigar[2 * MAX_READ_LEN], cigar[2 * MAX_READ_LEN], output_buf[5 * MAX_READ_LEN];
    vector<vector<int>> penalty(MAX_READ_LEN + 10), operation(MAX_READ_LEN + 10), start(MAX_READ_LEN + 10);

    add_time();
    const auto &chunk = chunks[chunk_i];
    int operations[3] = {1, -1, 0};

    penalty[0].resize(chunk_size + 10);
    operation[0].resize(chunk_size + 10);
    start[0].resize(chunk_size + 10);
    for (int i = 0; i < chunk.size + 1; ++i) {
        penalty[0][i] = 0;
        operation[0][i] = 0;
        start[0][i] = i;
    }

    for (const auto &chunk_read: chunks_reads[chunk_i]) {
        bool is_ct = chunk_read.second;
        int read_i = chunk_read.first;
        const auto &read = reads[read_i];
        for (int i = 1; i < read.size + 1; ++i) {
            if (penalty[i].empty()) {
                penalty[i].resize(chunk_size + 10);
                operation[i].resize(chunk_size + 10);
                start[i].resize(chunk_size + 10);
            }
            penalty[i][0] = gap_open_penalty + gap_extend_penalty * (i - 1);
            operation[i][0] = 1;
            start[i][0] = 0;
        }
        for (int i = 1; i < read.size + 1; ++i)
            for (int j = 1; j < chunk.size + 1; ++j) {
                int m_penalty;
                if (
                        read.seq_str[i - 1] == chunk.seq_str[j - 1] // XXX check upper lower case
                        ||
                        (is_ct and (read.seq_str[i - 1] == 'T'
//                        || read.seq_str[i - 1] == 't'
                        ) and
                         (chunk.seq_str[j - 1] == 'C'
//                         || chunk.seq_str[j - 1] == 'c'
                        ))
                        ||
                        (not is_ct and (read.seq_str[i - 1] == 'A'
//                        || read.seq_str[i - 1] == 'a'
                        ) and
                         (chunk.seq_str[j - 1] == 'G'
//                         || chunk.seq_str[j - 1] == 'g'
                        ))
                        )
                    m_penalty = -match_score;
                else
                    m_penalty = mismath_penalty;
                int penalties[3] = {
                        penalty[i - 1][j] + (operation[i - 1][j] == 0 ? gap_open_penalty : gap_extend_penalty),
                        penalty[i][j - 1] + (operation[i][j - 1] == 0 ? gap_open_penalty : gap_extend_penalty),
                        penalty[i - 1][j - 1] + m_penalty};
                int starts[3] = {start[i - 1][j], start[i][j - 1], start[i - 1][j - 1]};
                auto min_penalties_i = penalties[2] <= penalties[1] ?
                                       (penalties[2] <= penalties[0] ? 2 : 0) :
                                       (penalties[1] <= penalties[0] ? 1 : 0);
                penalty[i][j] = penalties[min_penalties_i];
                operation[i][j] = operations[min_penalties_i];
                start[i][j] = starts[min_penalties_i];
            }
        auto min_penalty = 10000000, min_j = 10000000;
        for (int j = 0; j < chunk.size + 1; ++j)
            if (penalty[read.size][j] < min_penalty) {
                min_penalty = penalty[read.size][j];
                min_j = j;
            }
        if (output_map_least_penalty[read_i] * ALT_RATIO_L2 < min_penalty)
            continue;
        auto read_pos = chunk.chr_pos + (chunk.chr_num > 0 ? start[read.size][min_j] : -min_j);
        auto read_chr = abs(chunk.chr_num);
        add_time();

        int read_len = static_cast<int>(read.size), reverse_cigar_pos = 0, cigar_pos = 0;
        while (read_len > 0) {
            switch (operation[read_len][min_j]) {
                case -1:
                    reverse_cigar[reverse_cigar_pos++] = 'D';
                    min_j -= 1;
                    break;
                case 1:
                    reverse_cigar[reverse_cigar_pos++] = 'I';
                    read_len -= 1;
                    break;
                case 0:
                    reverse_cigar[reverse_cigar_pos++] = 'M';
                    min_j -= 1;
                    read_len -= 1;
                    break;
                default:
                    break;
            }
        }
        int last_cigar_char_len = 1;
        char last_cigar_char = reverse_cigar[reverse_cigar_pos - 1];
        for (int i = reverse_cigar_pos - 2; i > -1; --i) {
            if (reverse_cigar[i] == last_cigar_char) {
                last_cigar_char_len++;
                continue;
            }
            cigar_pos += write_cigar_block(last_cigar_char_len, last_cigar_char, cigar + cigar_pos);
            last_cigar_char_len = 1;
            last_cigar_char = reverse_cigar[i];
        }
        cigar_pos += write_cigar_block(last_cigar_char_len, last_cigar_char, cigar + cigar_pos);
        cigar[cigar_pos] = '\0';
        add_time();
        sprintf(output_buf, "%s\t%d\t%d\t%lu\t%s\n", read.name, read_i, read_chr, read_pos, cigar);
        char *output_buf_copy = strdup(output_buf);
        pthread_mutex_lock(&output_map_inserting);
        if (output_map_least_penalty[read_i] * ALT_RATIO_L2 > min_penalty) {
            output_map[read_i].emplace_back(min_penalty, output_buf_copy);
            if (output_map_least_penalty[read_i] > min_penalty)
                output_map_least_penalty[read_i] = min_penalty;
        }
        pthread_mutex_unlock(&output_map_inserting);
    }
    remove_vector(chunks_reads[chunk_i]);
    add_time();
    logger->debug("align_chunk_reads_phase2 %s ms", get_times_str(true));
}

inline int align_chunk_reads_phase1(int chunk_i) {
    if (chunk_i % log_level_power == 0)
        logger->info("aligning chunk %d", chunk_i);
    if (chunks_reads[chunk_i].empty())
        return 0;
    add_time();
    align_chunk_reads_phase2(chunk_i);
    add_time();
    logger->debug("total align_chunk_reads_phase2 time: %d ms", last_time());
    return 1;
}

int main(int argc, char *argv[]) {
    read_args(argc, argv);
    config_str = Logger::formatString("%s_%d_%d_%d_%d_%d", ref_file_name, family_decompose_letters,
                                      kmer_length, MIN_HASH_COUNT, chunk_size, chunk_overlap);
    const char *config = config_str.c_str();
    add_time();
    logger->info("begin with config: %s", config);

    prepare_ref_sketch();
    add_time();
    logger->info("sketch loading/preparing time: %d ms", last_time());

    tie(reads, reads_should_be_deleted) = read_sequences_from_file(reads_file_name, FASTQ);
    if (to_read == 0)
        to_read = static_cast<int>(reads.size());
    add_time();
    logger->info("reads loading time: %d ms", last_time());

    reads_chunks.resize(reads.size());
    multiproc(threads_count, find_read_chunks, to_read, from_read);
    remove_vector<vector<tuple<int, int, int>>>(chunks_sketchs_CT);
    remove_vector<vector<tuple<int, int, int>>>(chunks_sketchs_GA);
    add_time();
    logger->info("find read chunks time: %d ms", last_time());

    chunks_reads.resize(chunks_size);
    for (int i = 0; i < reads_chunks.size(); ++i)
        for (auto &read_chunk: reads_chunks[i]) {
            bool is_ct = true;
            if (read_chunk >= chunks_size) {
                read_chunk -= chunks_size;
                is_ct = false;
            }
            chunks_reads[read_chunk].emplace_back(i, is_ct);
        }
    remove_vector<vector<int>>(reads_chunks);
    add_time();
    logger->info("assign reads to chunks: %d ms", last_time());

    if (chunks.empty()) {
        read_chunks();
        add_time();
        logger->info("read_chunks: %d ms", last_time());
    }
    output_map_least_penalty.resize(reads.size());
    output_map.resize(reads.size());
    multiproc(threads_count, align_chunk_reads_phase1, chunks_size);
    add_time();
    logger->info("align reads to chunks: %d ms", last_time());

    output_file = fopen(output_file_name, "w");
    for (int i = 0; i < output_map.size(); ++i) {
        int min_penalty = output_map_least_penalty[i];
        for (auto &p: output_map[i]) {
            if (p.first < ALT_RATIO_L2 * min_penalty)
                fprintf(output_file, "%s", p.second);
            free(p.second);
        }
    }
    fclose(output_file);
    add_time();
    logger->info("outputing: %d ms", last_time());
    free(ref_genome_should_be_deleted);
    free(reads_should_be_deleted);
    return 0;
}
