#include <getopt.h>
#include <cstring>
#include <cmath>
#include "params.h"
#include "utils/logger.h"


BisHashArgs read_args(int argc, char *argv[]) {
    BisHashArgs bis_hash_args{.log_level = Logger::INFO};

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
                bis_hash_args.threads_count = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'r':
                bis_hash_args.ref_file_name = strdup(optarg);
                break;
            case 'q':
                bis_hash_args.reads_file_name = strdup(optarg);
                break;
            case 'o':
                bis_hash_args.output_file_name = strdup(optarg);
                break;
            case 'm':
                bis_hash_args.alt_matchs_ratio = strtod(optarg, nullptr);
                break;
            case 'l':
                bis_hash_args.log_level = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'n':
                bis_hash_args.read_index = false;
                break;
            case 'w':
                bis_hash_args.write_index = false;
                break;
            case 'A':
                bis_hash_args.penalty_config.match_score = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'B':
                bis_hash_args.penalty_config.mismath_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'O':
                bis_hash_args.penalty_config.gap_open_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'E':
                bis_hash_args.penalty_config.gap_extend_penalty = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'F':
                bis_hash_args.from_read = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'T':
                bis_hash_args.to_read = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'D':
                bis_hash_args.family_decompose_letters = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'K':
                bis_hash_args.kmer_length = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'S':
                bis_hash_args.chunk_size = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'V':
                bis_hash_args.chunk_overlap = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            default:
                break;
        }
    optind = 1;
    return bis_hash_args;
}

