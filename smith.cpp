#include "smith.h"
#include "configs.h"
#include "utils/logger.h"
#include "utils/time_profile.h"
#include <vector>
#include <cstring>

using namespace std;

extern Logger *logger;
pthread_mutex_t output_map_inserting;

inline int write_cigar_block(int last_cigar_char_len, char last_cigar_char, char *cigar) {
    if (last_cigar_char_len == 1)
        return sprintf(cigar, "%c", last_cigar_char);
    else
        return sprintf(cigar, "%d%c", last_cigar_char_len, last_cigar_char);
}

void align_chunk_reads_phase2(const Sequence &chunk, const vector<pair<int, bool>> &chunk_reads, const vector<Sequence> &reads,
                              const unsigned long chr_size, const PenaltyConfig penalty_config,
                              vector<vector<pair<int, SamLine *>>> &output_map, vector<int> &output_map_least_penalty) {
    char reverse_cigar[2 * MAX_READ_LEN], cigar[2 * MAX_READ_LEN];
    vector<vector<int>> penalty(MAX_READ_LEN + 10), operation(MAX_READ_LEN + 10), start(MAX_READ_LEN + 10);

    add_time();
    const unsigned long max_extension = (unsigned long) MAX_READ_LEN / 4;
    const auto chunk_extension_amount_before = max(
            min(chunk.chr_num > 0 ? chunk.chr_pos : chr_size - chunk.chr_pos, max_extension), 0UL);
    const auto chunk_extension_amount_after = max(
            min((chunk.chr_num > 0 ? chr_size - chunk.chr_pos : chunk.chr_pos) - chunk.size, max_extension), 0UL);
    const int echunk_size = static_cast<int>(chunk_extension_amount_before + chunk.size + chunk_extension_amount_after);
    const auto echunk_seq_str = chunk.seq_str - chunk_extension_amount_before;
    const auto echunk_chr_pos =
            chunk.chr_num > 0 ? chunk.chr_pos - chunk_extension_amount_before : chunk.chr_pos + chunk_extension_amount_before;
    int operations[3] = {1, -1, 0};

    penalty[0].resize(echunk_size + 10);
    operation[0].resize(echunk_size + 10);
    start[0].resize(echunk_size + 10);
    for (int i = 0; i < echunk_size + 1; ++i) {
        penalty[0][i] = 0;
        operation[0][i] = 0;
        start[0][i] = i;
    }

    for (const auto &chunk_read: chunk_reads) {
        bool is_ct = chunk_read.second;
        int read_i = chunk_read.first;
        const auto &read = reads[read_i];
        for (int i = 1; i < read.size + 1; ++i) {
            if (penalty[i].empty()) {
                penalty[i].resize(echunk_size + 10);
                operation[i].resize(echunk_size + 10);
                start[i].resize(echunk_size + 10);
            }
            penalty[i][0] = penalty_config.gap_open_penalty + penalty_config.gap_extend_penalty * (i - 1);
            operation[i][0] = 1;
            start[i][0] = 0;
        }
        int min_penalty_j = 0, min_penalty = 0, chunk_from_index = 1, chunk_to_index = echunk_size + 1;
        int off = static_cast<int>(read.size) / 20;
        int delta_min = max(off / 10, 2) * penalty_config.gap_open_penalty;
        for (int i = 1; i < read.size + 1; ++i) {
            min_penalty_j = 0;
            min_penalty = penalty[i][0];
            bool narrowed = chunk_from_index != 1;
            for (int j = chunk_from_index; j < chunk_to_index; ++j) {
                bool narrowed_first_col = narrowed && j == chunk_from_index;
                int m_penalty;
                if (
                        read.seq_str[i - 1] == echunk_seq_str[j - 1] // XXX check upper lower case
                        ||
                        (is_ct and (read.seq_str[i - 1] == 'T'
//                        || read.seq_str[i - 1] == 't'
                        ) and
                         (echunk_seq_str[j - 1] == 'C'
//                         || echunk_seq_str[j - 1] == 'c'
                        ))
                        ||
                        (not is_ct and (read.seq_str[i - 1] == 'A'
//                        || read.seq_str[i - 1] == 'a'
                        ) and
                         (echunk_seq_str[j - 1] == 'G'
//                         || echunk_seq_str[j - 1] == 'g'
                        ))
                        )
                    m_penalty = -penalty_config.match_score;
                else
                    m_penalty = penalty_config.mismath_penalty;
                int penalties[3] = {
                        penalty[i - 1][j] +
                        (operation[i - 1][j] == 0 ? penalty_config.gap_open_penalty : penalty_config.gap_extend_penalty),
                        penalty[i][narrowed_first_col ? 0 : (j - 1)] +
                        (operation[i][j - 1] == 0 ? penalty_config.gap_open_penalty : penalty_config.gap_extend_penalty),
                        penalty[i - 1][narrowed_first_col ? 0 : (j - 1)] + m_penalty};
                int starts[3] = {start[i - 1][j], start[i][narrowed_first_col ? 0 : (j - 1)],
                                 start[i - 1][narrowed_first_col ? 0 : (j - 1)]};
                auto min_penalties_i = penalties[2] <= penalties[1] ?
                                       (penalties[2] <= penalties[0] ? 2 : 0) :
                                       (penalties[1] <= penalties[0] ? 1 : 0);
                penalty[i][j] = penalties[min_penalties_i];
                operation[i][j] = narrowed_first_col ? (operations[min_penalties_i] + 4) : operations[min_penalties_i];
                start[i][j] = starts[min_penalties_i];
                if (i > read.size / 10) {
                    if (min_penalty > penalty[i][j]) {
                        min_penalty = penalty[i][j];
                        min_penalty_j = j;
                    }
                }
            }
            if (min_penalty_j > 0) {
#ifdef _DEBUG
                logger->debugl2_noheader("%d\t%d\t%d\t", min_penalty_j, start[i][min_penalty_j], min_penalty);
                for (int j = 1; j < chunk_from_index; j++)
                    logger->debugl2_noheader("\t");
#endif
                int new_chunk_from_index = -1, new_chunk_to_index = 0;
                for (int j = chunk_from_index; j < chunk_to_index; ++j) {
                    if (penalty[i][j] - min_penalty < delta_min) {
                        if (new_chunk_from_index == -1)
                            new_chunk_from_index = j - off;
                        new_chunk_to_index = j + off;
                        logger->debugl2_noheader("->");
                    }
                    logger->debugl2_noheader("%d\t", penalty[i][j]);
                }
                chunk_from_index = max(new_chunk_from_index, 1);
                chunk_to_index = min(new_chunk_to_index, echunk_size + 1);
                logger->debugl2_noheader("\n");
            }
        }
        if (output_map_least_penalty[read_i] * ALT_RATIO_L2 < min_penalty)
            continue;
        auto read_pos = echunk_chr_pos + (chunk.chr_num > 0 ? start[read.size][min_penalty_j] : -min_penalty_j);
        add_time();

        int read_len = static_cast<int>(read.size), reverse_cigar_pos = 0, cigar_pos = 0;
        while (read_len > 0) {
            switch (operation[read_len][min_penalty_j]) {
                case -1:
                    reverse_cigar[reverse_cigar_pos++] = 'D';
                    min_penalty_j -= 1;
                    break;
                case 3:
                    reverse_cigar[reverse_cigar_pos++] = 'D';
                    min_penalty_j = 0;
                    break;

                case 5:
                case 1:
                    reverse_cigar[reverse_cigar_pos++] = 'I';
                    read_len -= 1;
                    break;

                case 0:
                    reverse_cigar[reverse_cigar_pos++] = 'M';
                    min_penalty_j -= 1;
                    read_len -= 1;
                    break;
                case 4:
                    reverse_cigar[reverse_cigar_pos++] = 'M';
                    min_penalty_j = 0;
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
        pthread_mutex_lock(&output_map_inserting);
        if (output_map_least_penalty[read_i] * ALT_RATIO_L2 > min_penalty) {
            SamLine *sam_line = SamLine::create_minimal_mapped_sam_line(&read, chunk.name, read_pos, strdup(cigar), chunk.chr_num < 0);
            output_map[read_i].emplace_back(min_penalty, sam_line);
            if (output_map_least_penalty[read_i] > min_penalty)
                output_map_least_penalty[read_i] = min_penalty;
        }
        pthread_mutex_unlock(&output_map_inserting);
    }
    add_time();
    logger->debug("align_chunk_reads_phase2 %s ms", get_times_str(true));
}
