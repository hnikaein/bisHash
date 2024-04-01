#include <getopt.h>
#include <cstring>
#include <cmath>
#include "params.h"
#include "utils/logger.h"

using namespace std;

struct CompleteOption {
    const char short_name;
    int has_arg;
    const char *long_name;
    const char *comment;
};
CompleteOption complete_options[] =
        {
                {'r', required_argument, "ref",                      "Reference file path (in FASTA format)."},
                {'q', required_argument, "query",                    "Query file path (in FASTQ format)."},
                {'o', required_argument, "output",                   "Output file path (in SAM format)."},
                {'t', required_argument, "threads",                  "Number of threads to use (default: 1)."},
                {'l', required_argument, "log",                      "Logging level (default: 4)."},
                {'f', required_argument, "from-read",                "Starting read for matching reads, inclusive (default: 0)."},
                {'e', required_argument, "to-read",                  "Ending read for matching reads, exclusive (default: INT_MAX)."},
                {'n', no_argument,       "no-read-index",            "Skip reading index (default: not set)."},
                {'w', no_argument,       "no-write-index",           "Skip writing index (default: not set)."},
                {'i', no_argument,       "only-create-index",        "Create index mode (default: not set)."},
                {'M', required_argument, "alt-match-ratio",          "Penalty ratio for alternative matches (default: 0.98)."},
                {'A', required_argument, "match-score",              "Score for a match (default: 5)."},
                {'B', required_argument, "mismatch-penalty",         "Penalty for a mismatch (default: 5)."},
                {'O', required_argument, "gap-open-penalty",         "Penalty for gap opening (default: 5)."},
                {'E', required_argument, "gap-extend-penalty",       "Penalty for gap extension (default: 3)."},
                {'D', required_argument, "family-decompose-letters", "Number of basepairs for decomposing family minhash (default: 5)."},
                {'K', required_argument, "kmer-length",              "Length of k-mer for family minhash evaluation (default: 20)."},
                {'S', required_argument, "chunk-size",               "Size of reference chunks (default: 7000)."},
                {'V', required_argument, "chunk-overlap",            "Overlap between chunks to ensure reads are contained (default: 2000)."},
                {'h', no_argument,       "help",                     "Print this help message."},
        };

void print_usage() {
    printf("bisHash 0.6.0\n");
    printf("Usage: bisHash -r reference_file -q query_file -o output_file [OPTIONS]\n");
    printf("Usage: bisHash -i -r reference_file [OPTIONS]\n");
    printf("\n");
    for (auto &complete_option: complete_options) {
        string param_naming;
        if (complete_option.short_name > '1' && complete_option.short_name < 'z')
            param_naming += string("-") + complete_option.short_name;
        if (complete_option.long_name != nullptr) {
            if (!param_naming.empty())
                param_naming += ", ";
            else
                param_naming += "    ";
            param_naming += string("--") + complete_option.long_name;
        }
        printf("  %-30s\t\t\t%s\n", param_naming.c_str(), complete_option.comment);
    }
}

BisHashArgs read_args(int argc, char *argv[]) {
    BisHashArgs bis_hash_args{.log_level = Logger::INFO};

    string shortopts;
    for (auto &complete_option: complete_options) {
        shortopts += complete_option.short_name;
        shortopts += (complete_option.has_arg == required_argument) ? ":" : (complete_option.has_arg == optional_argument) ? "::" : "";
    }
    auto complete_options_length = size(complete_options);
    auto longopts = new option[complete_options_length];
    for (int i = 0; i < complete_options_length; ++i)
        longopts[i] = {complete_options[i].long_name, complete_options[i].has_arg, nullptr, complete_options[i].short_name};
    int option_index = 0, c;
    while ((c = getopt_long(argc, argv, shortopts.c_str(), longopts, &option_index)) >= 0)
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
            case 'M':
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
            case 'f':
                bis_hash_args.from_read = static_cast<int>(strtol(optarg, nullptr, 10));
                break;
            case 'e':
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
            case 'h':
                print_usage();
                exit(EXIT_SUCCESS);
            case 'i':
                bis_hash_args.only_create_index = true;
                break;
            default:
                print_usage();
                throw runtime_error(string("unmatched argument"));
        }
    delete[] longopts;
    if (bis_hash_args.ref_file_name == nullptr ||
        (!bis_hash_args.only_create_index && (bis_hash_args.reads_file_name == nullptr || bis_hash_args.output_file_name == nullptr))) {
        fprintf(stderr, "you should specify reference file, reads file and output file.\n");
        print_usage();
        throw runtime_error(string("you should specify reference file, reads file and output file."));
    }
    return bis_hash_args;
}

