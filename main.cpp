#include <algorithm>
#include <map>
#include <cmath>
#include <climits>
#include <set>
#include <cstring>
#include "configs.h"
#include "utils/multiproc.h"
#include "utils/logger.h"
#include "utils/time_profile.h"
#include "bio_utils/sequence_reader.h"
#include "bio_utils/sam_writer.h"
#include "hashing/family_min_hash.h"
#include "indexing/object_writer.h"
#include "smith.h"
#include "params.h"

using namespace std;

extern Logger *logger;

FamilyMinHash *family_min_hash;

vector<Sequence> ref_genome, chunks, reads;
vector<pair<map<int, vector<int>>, map<int, vector<int>>>> chunks_sketchs_CT_tree_vector, chunks_sketchs_GA_tree_vector;
vector<pair<vector<int *>, vector<int *>>> chunks_sketchs_CT_tree, chunks_sketchs_GA_tree;
vector<vector<int>> reads_chunks;
vector<vector<pair<int, bool>>> chunks_reads;
vector<vector<pair<int, SamLine *>>> output_map;
vector<int> output_map_least_penalty;

FILE *output_file;
string config_str;

BisHashArgs args;
int log_level_power;
int chunks_count;
char *ref_genome_should_be_deleted, *reads_should_be_deleted;
int *sketchs_should_be_deleted;

template<typename T>
void remove_vector(vector<T> &vec) {
    vector<T>().swap(vec);
}

pthread_mutex_t chunks_sketchs_CT_tree_vector_mutex, chunks_sketchs_GA_tree_vector_mutex;

// Calculate the chunk sketch and add the chunk id to hashes that exist in chunk sketch
int make_chunk_sketch(const int chunk_i) {
    auto chunk_sketch_CT = family_min_hash->get_sketch(chunks[chunk_i].seq_str,
                                                       static_cast<int>(chunks[chunk_i].size),
                                                       SKETCH_MODE_WITH_C_T_CONVERSION);
    for (const auto &part: chunk_sketch_CT) {
        pthread_mutex_lock(&chunks_sketchs_CT_tree_vector_mutex);
        chunks_sketchs_CT_tree_vector[get<0>(part)].first[get<1>(part)].push_back(chunk_i);
        chunks_sketchs_CT_tree_vector[get<0>(part)].second[get<2>(part)].push_back(chunk_i);
        pthread_mutex_unlock(&chunks_sketchs_CT_tree_vector_mutex);
    }
    auto chunk_sketch_GA = family_min_hash->get_sketch(chunks[chunk_i].seq_str,
                                                       static_cast<int>(chunks[chunk_i].size),
                                                       SKETCH_MODE_WITH_G_A_CONVERSION);
    for (const auto &part: chunk_sketch_GA) {
        pthread_mutex_lock(&chunks_sketchs_GA_tree_vector_mutex);
        chunks_sketchs_GA_tree_vector[get<0>(part)].first[get<1>(part)].push_back(chunk_i);
        chunks_sketchs_GA_tree_vector[get<0>(part)].second[get<2>(part)].push_back(chunk_i);
        pthread_mutex_unlock(&chunks_sketchs_GA_tree_vector_mutex);
    }
    return 1;
}

void read_chunks() {
    tie(ref_genome, ref_genome_should_be_deleted) = read_sequences_from_file(args.ref_file_name, FASTA, true);
    chunks = Sequence::chunkenize_big_sequence(ref_genome, args.chunk_size, true, args.chunk_size - args.chunk_overlap);
}

void prepare_ref_sketch() {
    auto index_file_name = config_str + ".bsh";
    add_time();
    if (args.read_index)
        try {
            auto indexes = read_data(index_file_name.c_str(), sketchs_should_be_deleted, chunks_count);
            chunks_sketchs_CT_tree = std::move(indexes[0]);
            chunks_sketchs_GA_tree = std::move(indexes[1]);
            delete[] indexes;
            logger->info("index read completed");
        } catch (...) {
            args.read_index = false;
        }
    if (!args.read_index) {
        read_chunks();
        chunks_count = static_cast<int>(chunks.size());
        chunks_sketchs_CT_tree_vector.resize((int) pow(MAX_BASENUMBER, args.family_decompose_letters));
        chunks_sketchs_GA_tree_vector.resize((int) pow(MAX_BASENUMBER, args.family_decompose_letters));
        multiproc(args.threads_count, make_chunk_sketch, chunks_count);
        vector<pair<map<int, vector<int>>, map<int, vector<int>>>> *indexes[] = {&chunks_sketchs_CT_tree_vector,
                                                                                 &chunks_sketchs_GA_tree_vector};
        write_data(index_file_name.c_str(), indexes, 2, chunks_count);
        remove_vector(chunks_sketchs_CT_tree_vector);
        remove_vector(chunks_sketchs_GA_tree_vector);
        auto indexes2 = read_data(index_file_name.c_str(), sketchs_should_be_deleted, chunks_count);
        chunks_sketchs_CT_tree = std::move(indexes2[0]);
        chunks_sketchs_GA_tree = std::move(indexes2[1]);
        delete[] indexes2;
    }
    add_time();
    logger->info("ref sketchs prepared: %sms: %d records", get_times_str(true), chunks_count);
}

inline int *find_in_vector(const vector<int *> &data, int find_what) {
    int low = 0, high = static_cast<int>(data.size());
    while (low < high) {
        int mid = (high + low) / 2;
        if (find_what > data[mid][0])
            low = mid + 1;
        if (find_what < data[mid][0])
            high = mid;
        if (find_what == data[mid][0])
            return data[mid];
    }
    return nullptr;
}

int find_read_chunks(const int read_i) {
    auto read = &reads[read_i];

    auto read_sketch_ct = family_min_hash->get_sketch(read->seq_str, static_cast<int>(read->size),
                                                      SKETCH_MODE_WITH_C_T_CONVERSION);
    auto read_sketch_ga = family_min_hash->get_sketch(read->seq_str, static_cast<int>(read->size),
                                                      SKETCH_MODE_WITH_G_A_CONVERSION);
    int similarit_length = chunks_count * 2 + 10;
    int *similarity = new int[similarit_length];
    memset(similarity, 0, sizeof(int) * similarit_length);
    int *similarity_l = similarity;
    auto read_sketch = &read_sketch_ct;
    auto chunks_sketchs_tree = &chunks_sketchs_CT_tree;
    for (int sketch_type = 0; sketch_type < 2; ++sketch_type) {
        if (sketch_type == 1) {
            read_sketch = &read_sketch_ga;
            chunks_sketchs_tree = &chunks_sketchs_GA_tree;
            similarity_l = similarity + chunks_count;
        }
        for (const auto &read_sketch_i: *read_sketch) {
            set<int> matches;
            auto candid_chunks = find_in_vector((*chunks_sketchs_tree)[get<0>(read_sketch_i)].first, get<1>(read_sketch_i));
            if (candid_chunks != nullptr) {
                candid_chunks += 1;
                for (int i = 1; i <= candid_chunks[0]; ++i)
                    matches.insert(candid_chunks[i]);
            }

            candid_chunks = find_in_vector((*chunks_sketchs_tree)[get<0>(read_sketch_i)].second, get<2>(read_sketch_i));
            if (candid_chunks != nullptr) {
                candid_chunks += 1;
                for (int i = 1; i <= candid_chunks[0]; ++i)
                    matches.insert(candid_chunks[i]);
            }

            for (const auto &match: matches)
                similarity_l[match] += 1;
        }
    }
    int last_inserted_j = 0;
    vector<int> max_sim_i(ALT_MATCHS);
    vector<int> max_sim(ALT_MATCHS);
    for (int i = 0; i < chunks_count * 2; ++i) {
        if (max_sim.size() > 3 * ALT_MATCHS) {
            max_sim.resize(2 * ALT_MATCHS);
            max_sim.reserve(4 * ALT_MATCHS);
            max_sim_i.resize(2 * ALT_MATCHS);
            max_sim_i.reserve(4 * ALT_MATCHS);
        }
        int sketch_sim = similarity[i];
        if (sketch_sim <= max_sim[ALT_MATCHS - 1])
            continue;
        if (max_sim_i[last_inserted_j] == i - 1) {
            if (max_sim[last_inserted_j] < sketch_sim) {
                max_sim_i.erase(max_sim_i.begin() + last_inserted_j);
                max_sim.erase(max_sim.begin() + last_inserted_j);
            } else {
                if (max_sim[last_inserted_j] == sketch_sim)
                    max_sim_i[last_inserted_j] = i;
                continue;
            }
        }
        for (int j = 0; j < ALT_MATCHS; ++j)
            if (sketch_sim > max_sim[j]) {
                max_sim_i.insert(max_sim_i.begin() + j, i);
                max_sim.insert(max_sim.begin() + j, sketch_sim);
                last_inserted_j = j;
                break;
            }
    }
    int idx = 0;
    if (max_sim[0] != 0)
        for (; idx < max_sim.size(); ++idx)
            if (args.alt_matchs_ratio * max_sim[0] > max_sim[idx])
                break;
    max_sim_i.resize(min(idx, ALT_MATCHS));
    reads_chunks[read_i] = std::move(max_sim_i);
    delete[] similarity;
    if (max_sim_i.empty())
        return 0;
    return 1;
}

inline int align_chunk_reads_phase1(int chunk_i) {
    if (chunk_i % log_level_power == 0)
        logger->info("aligning chunk %d", chunk_i);
    if (chunks_reads[chunk_i].empty())
        return 0;
    add_time();
    unsigned long chr_size = ref_genome[abs(chunks[chunk_i].chr_num) - 1].size;
    align_chunk_reads_phase2(chunks[chunk_i], chunks_reads[chunk_i], reads, chr_size, args.penalty_config, output_map,
                             output_map_least_penalty);
    remove_vector(chunks_reads[chunk_i]);
    add_time();
    logger->debug("total align_chunk_reads_phase2 time: %d ms", last_time());
    return 1;
}

int main(int argc, char *argv[]) {
    try {
        args = read_args(argc, argv);
    } catch (runtime_error &e) {
        return EXIT_FAILURE;
    }
    family_min_hash = new FamilyMinHash(args.family_decompose_letters, args.kmer_length, MAX_BASENUMBER, BIG_PRIME_NUMBER, MIN_HASH_COUNT);
    log_level_power = max(static_cast<int>(pow(10, 8 - args.log_level)), 1);
    logger = new Logger(args.log_level);
    config_str = Logger::formatString("%s_%d_%d_%d_%d_%d", args.ref_file_name, args.family_decompose_letters,
                                      args.kmer_length, MIN_HASH_COUNT, args.chunk_size, args.chunk_overlap);
    const char *config = config_str.c_str();
    add_time();
    logger->info("begin with config: %s", config);

    prepare_ref_sketch();
    if (args.only_create_index)
        return EXIT_SUCCESS;
    add_time();
    logger->info("sketch loading/preparing time: %d ms", last_time());

    tie(reads, reads_should_be_deleted) = read_sequences_from_file(args.reads_file_name, FASTQ);
    args.to_read = min(static_cast<int>(reads.size()), args.to_read);
    add_time();
    logger->info("reads loading time: %d ms", last_time());

    reads_chunks.resize(reads.size());
    multiproc(args.threads_count, find_read_chunks, args.to_read, args.from_read);
    remove_vector(chunks_sketchs_CT_tree);
    remove_vector(chunks_sketchs_GA_tree);
    delete[] sketchs_should_be_deleted;
    add_time();
    logger->info("find read chunks time: %d ms", last_time());

    chunks_reads.resize(chunks_count);
    for (int i = 0; i < reads_chunks.size(); ++i)
        for (auto &read_chunk: reads_chunks[i]) {
            bool is_ct = true;
            if (read_chunk >= chunks_count) {
                read_chunk -= chunks_count;
                is_ct = false;
            }
            chunks_reads[read_chunk].emplace_back(i, is_ct);
        }
    remove_vector(reads_chunks);
    add_time();
    logger->info("assign reads to chunks: %d ms", last_time());

    if (chunks.empty()) {
        read_chunks();
        add_time();
        logger->info("read_chunks: %d ms", last_time());
    }
    output_map_least_penalty.resize(reads.size(), INT_MAX);
    output_map.resize(reads.size());
    multiproc(args.threads_count, align_chunk_reads_phase1, chunks_count);
    add_time();
    logger->info("align reads to chunks: %d ms", last_time());

    output_file = fopen(args.output_file_name, "w");
    for (int i = args.from_read; i < args.to_read; ++i) {
        if (output_map[i].empty()) {
            SamLine::create_unmapped_sam_line(&reads[i]).print_to_file(output_file);
            continue;
        }
        int min_penalty = output_map_least_penalty[i];
        sort(output_map[i].begin(), output_map[i].end());
        for (int j = 0; j < output_map[i].size(); ++j) {
            auto &p = output_map[i][j];
            if (j < ALT_MATCHS && p.first < ALT_RATIO_L2 * min_penalty) {
                if (j != 0)
                    p.second->set_as_secondary();
                if (p.second->is_reversed() && reads[i].reverse_seq_str == nullptr)
                    reads[i].create_reverse();
                p.second->print_to_file(output_file);
            }
            free((char *) p.second->cigar);
            delete p.second;
        }
    }
    fclose(output_file);
    add_time();
    logger->info("outputing: %d ms", last_time());
    delete[] ref_genome_should_be_deleted;
    delete[] reads_should_be_deleted;
    return 0;
}
