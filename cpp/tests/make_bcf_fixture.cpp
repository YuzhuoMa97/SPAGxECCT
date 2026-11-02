// SPDX-License-Identifier: GPL-3.0-or-later
// Test fixture utility only: convert a VCF to BCF using the vendored htslib.
// This does not add a conversion mode to the SPAGxE analysis CLI.
#include <htslib/vcf.h>
#include <cstdio>

int main(int argc, char **argv) {
    if (argc != 3) {
        std::fprintf(stderr, "Usage: make_bcf_fixture INPUT.vcf OUTPUT.bcf\n");
        return 1;
    }
    htsFile *input = bcf_open(argv[1], "r");
    if (!input) { std::fprintf(stderr, "Cannot open VCF input: %s\n", argv[1]); return 1; }
    bcf_hdr_t *header = bcf_hdr_read(input);
    if (!header) { std::fprintf(stderr, "Cannot read VCF header.\n"); bcf_close(input); return 1; }
    htsFile *output = bcf_open(argv[2], "wb");
    if (!output) { std::fprintf(stderr, "Cannot create BCF output: %s\n", argv[2]); bcf_hdr_destroy(header); bcf_close(input); return 1; }
    bcf1_t *record = bcf_init();
    int status = record && bcf_hdr_write(output, header) == 0 ? 0 : 1;
    int readStatus = -1;
    size_t count = 0;
    if (!status) {
        while ((readStatus = bcf_read(input, header, record)) == 0) {
            if (record->errcode || bcf_write(output, header, record) != 0) { status = 1; break; }
            ++count;
        }
        if (readStatus < -1) status = 1;
    }
    if (record) bcf_destroy(record);
    bcf_hdr_destroy(header);
    if (bcf_close(input) != 0) status = 1;
    if (bcf_close(output) != 0) status = 1;
    if (status) { std::fprintf(stderr, "VCF to BCF conversion failed.\n"); return 1; }
    std::printf("Wrote %zu BCF variants to %s\n", count, argv[2]);
    return 0;
}
