# bisHash

## Introduction

bisHash is an advanced tool for aligning bisulfite-treated DNA sequences against a reference genome. This tool is developed in C++ to ensure accuracy and reliability in handling bisulfite sequencing data. While bisHash may not be the fastest aligner available, it is designed to provide highly accurate alignments, making it particularly suitable for studies where precision is paramount.

## Installation

### Prerequisites

Before installing bisHash, ensure you have the following prerequisites:

- C++20 or later
- CMake 3.16 or later
- GCC or Clang compiler

### Build Instructions

1. Clone the repository:

   ```bash
   git clone https://github.com/hnikaein/bisHash.git
   cd bisHash
   ```

2. Create a build directory and navigate into it:

   ```bash
   mkdir build
   cd build
   ```

3. Run CMake to configure the build:

   ```bash
   cmake ..
   ```

4. Compile the source code:

   ```bash
   make
   ```

## Usage

Once maked, you can run bisHash with the following command:

```bash
./bisHash -r reference_file -q query_file -o output_file [OPTIONS]
```


### Command Line Options

- `-r` or `--ref`: Reference file path (in FASTA format).
- `-q` or `--query`: Query file path (in FASTQ format).
- `-o` or `--output`: Output file path (in SAM format).
- `-t` or `--threads`: Number of threads to use (default: 1).
- 
- `-l` or `--log`: Logging level (default: 4).
- `-f` or `--from-read`: Starting read for matching reads, inclusive (default: 0).
- `-e` or `--to-read`: Ending read for matching reads, exclusive (default: INT_MAX).
- 
- `-n` or `--no-read-index`: Skip reading index (default: not set).
- `-i` or `--only-create-index`: Create index mode (default: not set).
- 
- `-M` or `--alt-match-ratio`: Penalty ratio for alternative matches (default: 0.98).
- `-A` or `--match-score`: Score for a match (default: 5).
- `-B` or `--mismatch-penalty`: Penalty for a mismatch (default: 5).
- `-O` or `--gap-open-penalty`: Penalty for gap opening (default: 5).
- `-E` or `--gap-extend-penalty`: Penalty for gap extension (default: 3).
- 
- `-D` or `--family-decompose-letters`: Number of basepairs for decomposing family minhash (default: 5).
- `-K` or `--kmer-length`: Length of k-mer for family minhash evaluation (default: 16).
- `-S` or `--chunk-size`: Size of reference chunks (default: 1000).
- `-V` or `--chunk-overlap`: Overlap between chunks to ensure reads are contained (default: 300).
- 
- `-h` or `--help`: Print this help message.

### Example

```bash
./bisHash -r hg19.fa -q sample_reads.fq -o aligned.sam -t 4
```

This command aligns the reads in `sample_reads.fq` to the `hg19.fa` reference genome using 4 threads and outputs the results to `aligned.sam`.

## Simulation Scripts

This project includes simulation scripts located in the simulation folder, which are useful for generating test inputs and comparing results.

### Running Simulations
1. `runAlgos.sh`: This script simulates inputs and runs various tools alongside bisHash, processing the results generated. You can configure the parameters at the beginning of this script.
2. `run_sam_read_compare.sh`: This script processes the output from runAlgos.sh to create comparison results, allowing for an evaluation of bisHash’s performance against other tools. Configuration options can also be set at the beginning of this script.

**Make sure to review and adjust the configuration parameters in both scripts to fit your testing requirements.**

## Contributing

Contributions to the project are welcome! Please follow these steps to contribute:

1. Fork the repository.
2. Create a new branch with a descriptive name (`git checkout -b feature-branch-name`).
3. Make your changes and commit them (`git commit -m "Description of changes"`).
4. Push to your branch (`git push origin feature-branch-name`).
5. Open a Pull Request on GitHub.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

If you have any questions or issues, please open an issue on the repository or contact the project maintainer at [nikaein@gmail.com](mailto:nikaein@gmail.com).
