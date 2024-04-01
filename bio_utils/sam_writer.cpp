#include "sam_writer.h"

const char STAR[2] = {'*', '\0'};

SamLine SamLine::create_unmapped_sam_line(const Sequence *read) {
    return SamLine{.read = read, .flag = 4, .rname = STAR, .pos = 0, .mapq = 0, .cigar = STAR, .rnext = STAR, .pnext = 0, .tlen = 0,
            .qual = STAR};
}

SamLine *SamLine::create_minimal_mapped_sam_line(const Sequence *read, const char *rname, const unsigned long pos, const char *cigar,
                                                 bool is_reversed, bool is_secondary) {
    return new SamLine{.read = read, .flag = (is_reversed ? 0x10 : 0) + (is_secondary ? 0x100 : 0), .rname = rname, .pos = pos, .mapq = 0,
            .cigar = cigar, .rnext = STAR, .pnext = 0, .tlen = 0, .qual = STAR};
}

void SamLine::print_to_file(FILE *output_file) const {
    const char *seq = (flag & 0x4) ? STAR : (flag & 0x10) ? read->reverse_seq_str : read->seq_str;
    fprintf(output_file, "%s\t%d\t%s\t%lu\t%d\t%s\t%s\t%lu\t%lu\t%s\t%s\n", read->name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen,
            seq, qual);
}

void SamLine::set_as_secondary() {
    flag += 0x100;
}

bool SamLine::is_reversed() {
    return flag & 0x10;
}

bool SamLine::operator<(const SamLine &rhs) const {
    return flag < rhs.flag;
}
