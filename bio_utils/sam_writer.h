#ifndef BISHASH_SAM_WRITER_H
#define BISHASH_SAM_WRITER_H

#include <string>
#include "sequence.h"

struct SamLine {
    const Sequence *read; // instead of const char *qname; and const char *seq;
    int flag;
    const char *rname;
    const unsigned long pos;
    const int mapq;
    const char *cigar;
    const char *rnext;
    const unsigned long pnext;
    const unsigned long tlen;
    const char *qual;

    static SamLine create_unmapped_sam_line(const Sequence *read);

    static SamLine *create_minimal_mapped_sam_line(const Sequence *read, const char *rname, unsigned long pos, const char *cigar,
                                                   bool is_reversed, bool is_secondary = false);

    void print_to_file(FILE *output_file) const;

    void set_as_secondary();

    bool is_reversed();

    bool operator<(const SamLine &rhs) const;

};

#endif //BISHASH_SAM_WRITER_H
