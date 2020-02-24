
#include <stdio.h>
#include <iostream>
#include <fstream>
// #include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>
#include <map>
#include <set>
#include <utility>

#include "htslib/sam.h"
#include "htslib/bgzf.h"


int main(int argc, char *argv[])
{

    samFile *in;  // open input alignment file
    bam_hdr_t *input_header; // alignment header
    std::ofstream out;
    int exit_code = 0;
    
    const char *in_name;
    const char *out_name = "-";

    int c;  // for parsing input arguments

    // while ((c = getopt(argc, argv, "DSIt:i:bCul:o:N:BZ:@:M")) >= 0) {
    while ((c = getopt(argc, argv, "i:o:")) >= 0) {
        switch(c) {
        case 'i': in_name = optarg; break;
        case 'o': out_name = optarg; break;
        default: /* '?' */
            fprintf(stderr, "Usage: sam_to_gtf -i <in.bam>|<in.bam> -o output_file\n");
            fprintf(stderr, "\n");
            fprintf(stderr, "-i: Input sam/bam file name (required).\n");
            fprintf(stderr, "-o: Output file name (required).\n");
            fprintf(stderr, "File to process.\n");
            fprintf(stderr, "\n");
            exit(EXIT_FAILURE);
        }
    }

    out.open(out_name);

    in = sam_open(in_name, "r");

    // Enable multi-threading (only effective when the library was compiled with -DBGZF_MT)
    // bgzf_mt(in, n_threads, 0);

    if (in == NULL) {
        fprintf(stderr, "Error opening \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    input_header = sam_hdr_read(in);
    if (input_header == NULL) {
        fprintf(stderr, "Couldn't read header for \"%s\"\n", argv[optind]);
        return EXIT_FAILURE;
    }

    bam1_t *aln = bam_init1(); //initialize an alignment

    float average_similarity_proportion_gap_compressed = 0.0f;
    float average_similarity_proportion_blast = 0.0f;
    float average_similarity_proportion_blast_with_soft_clipping = 0.0f;
    float average_identical_proportion_nogaps = 0.0f;
    float average_identical_proportion_gap_compressed = 0.0f;
    float average_identical_proportion_blast = 0.0f;
    float average_identical_proportion_with_soft_clipping = 0.0f;
    float average_insertion_gap_rate = 0.0f;
    float average_insertion_proportion_blast = 0.0f;
    float average_deletion_gap_rate = 0.0f;
    float average_deletion_proportion_blast = 0.0f;
    float average_proportion_soft_clip = 0.0f;
    
    int_fast64_t i = 0;
    // int_fast64_t missing_nm_i = 0;


    out <<
            "read_name" << "\t" << "total_length" << "\t" <<
            "similarity_length" << "\t" << "similarity_proportion_gap_compressed" << "\t" << "similarity_proportion_blast" << "\t" << "similarity_proportion_blast_with_soft_clipping" << "\t" <<
            "identical_length" << "\t" << "identical_proportion_nogaps" << "\t" << "identical_proportion_gap_compressed" << "\t" << "identical_proportion_blast" << "\t" << "identical_proportion_with_soft_clipping" << "\t" <<
            "n_insertions" << "\t" << "insertion_length" << "\t" << "insertion_gap_rate" << "\t" << "insertion_overall_proportion" << "\t" <<
            "n_deletions" << "\t" << "deletion_length" << "\t" << "deletion_gap_rate" << "\t" << "deletion_overall_proportion" << "\t" <<
            "soft_clip_length" << "\t" << "proportion_soft_clip" <<
            std::endl;

    while(sam_read1(in, input_header, aln) >= 0) {

        // make sure the read is mapped
        if ((aln->core.flag & BAM_FUNMAP) != 0)
            continue;

        // make sure the alignment is not a secondary alignment or a supplementary alignment
        if ((aln->core.flag & BAM_FSECONDARY) != 0 or (aln->core.flag & BAM_FSUPPLEMENTARY) != 0)
            continue;

        ++i;
        int_fast64_t similar_length = 0;
        int_fast64_t insertion_length = 0;
        int_fast64_t deletion_length = 0;
        int_fast64_t soft_clip_length = 0;
        int_fast32_t n_insertions = 0;
        int_fast32_t n_deletions = 0;

        // get cigar
        u_int32_t *cigar = bam_get_cigar(aln);

        for (u_int32_t k = 0; k < aln->core.n_cigar; ++k) {
            switch(bam_cigar_op(cigar[k])) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:
                    similar_length += bam_cigar_oplen(cigar[k]);
                    break;
                case BAM_CINS:
                    insertion_length += bam_cigar_oplen(cigar[k]);
                    ++n_insertions;
                    break;
                case BAM_CDEL:
                    deletion_length += bam_cigar_oplen(cigar[k]);
                    ++n_deletions;
                    break;
                case BAM_CSOFT_CLIP:
                    soft_clip_length += bam_cigar_oplen(cigar[k]);
                    break;
                default:
                    break;
            }
        }

        int64_t nm_length = 0;
        uint8_t *nm_tag = bam_aux_get(aln, "NM");
        if (nm_tag != NULL) {
            nm_length = bam_aux2i(nm_tag);
        }
        // else {
        //     ++missing_nm_i;
        // }

        // float df_prop = -1.0f;
        // uint8_t *df_tag = bam_aux_get(aln, "df");
        // if (df_tag != NULL) {
        //     df_prop = bam_aux2f(df_tag);
        // }

        int_fast64_t total_length = similar_length + insertion_length + deletion_length + soft_clip_length;
        int_fast64_t clipless_total_length = similar_length + insertion_length + deletion_length;
        int_fast64_t identical_length = similar_length - nm_length + insertion_length + deletion_length;
        int_fast64_t gap_compressed_length = similar_length + n_insertions + n_deletions;


        float similarity_proportion_gap_compressed = (float)similar_length / gap_compressed_length;
        average_similarity_proportion_gap_compressed += similarity_proportion_gap_compressed;
        float similarity_proportion_blast = (float)similar_length / clipless_total_length;
        average_similarity_proportion_blast += similarity_proportion_blast;
        float similarity_proportion_blast_with_soft_clipping = ((float)similar_length) / total_length;
        average_similarity_proportion_blast_with_soft_clipping += similarity_proportion_blast_with_soft_clipping;
        float identical_proportion_nogaps = (float)identical_length / similar_length;
        average_identical_proportion_nogaps += identical_proportion_nogaps;
        float identical_proportion_gap_compressed = (float)identical_length / gap_compressed_length;
        average_identical_proportion_gap_compressed += identical_proportion_gap_compressed;
        float identical_proportion_blast = (float)identical_length / clipless_total_length;
        average_identical_proportion_blast += identical_proportion_blast;
        float identical_proportion_with_soft_clipping = (float)(identical_length + soft_clip_length) / total_length;
        average_identical_proportion_with_soft_clipping += identical_proportion_with_soft_clipping;
        float insertion_gap_rate = (float)n_insertions / gap_compressed_length;
        average_insertion_gap_rate += insertion_gap_rate;
        float insertion_proportion_blast = (float)insertion_length / clipless_total_length;
        average_insertion_proportion_blast += insertion_proportion_blast;
        float deletion_gap_rate = (float)n_deletions / gap_compressed_length;
        average_deletion_gap_rate += deletion_gap_rate;
        float deletion_proportion_blast = (float)deletion_length / clipless_total_length;
        average_deletion_proportion_blast += deletion_proportion_blast;
        float proportion_soft_clip = (float)soft_clip_length / total_length;
        average_proportion_soft_clip += proportion_soft_clip;

        // output the alignment
        // similar_proportion is blast proportion, since gap excluded similarity would always be 100%
        out <<
                bam_get_qname(aln) << "\t" << aln->core.l_qseq << "\t" <<
                similar_length << "\t" << similarity_proportion_gap_compressed << "\t" << similarity_proportion_blast << "\t" << similarity_proportion_blast_with_soft_clipping << "\t" <<
                identical_length << "\t" << identical_proportion_nogaps << "\t" << identical_proportion_gap_compressed << "\t" << identical_proportion_blast << "\t" << identical_proportion_with_soft_clipping << "\t" <<
                n_insertions << "\t" << insertion_length << "\t" << insertion_gap_rate << "\t" << insertion_proportion_blast << "\t" <<
                n_deletions << "\t" << deletion_length << "\t" << deletion_gap_rate << "\t" << deletion_proportion_blast << "\t" <<
                soft_clip_length << "\t" << proportion_soft_clip << "\t" <<
                std::endl; 
    }

    average_similarity_proportion_gap_compressed /= i;
    average_similarity_proportion_blast /= i;
    average_similarity_proportion_blast_with_soft_clipping /= i;
    average_identical_proportion_nogaps /= i;
    average_identical_proportion_gap_compressed /= i;
    average_identical_proportion_blast /= i;
    average_identical_proportion_with_soft_clipping /= i;
    average_insertion_gap_rate /= i;
    average_insertion_proportion_blast /= i;
    average_deletion_gap_rate /= i;
    average_deletion_proportion_blast /= i;
    average_proportion_soft_clip /= i;

    out <<
            "total" << "\t" << "total_length" << "\t" <<
            "similarity_length" << "\t" << average_similarity_proportion_gap_compressed << "\t" << average_similarity_proportion_blast << "\t" << average_similarity_proportion_blast_with_soft_clipping << "\t" <<
            "identical_length" << "\t" << average_identical_proportion_nogaps << "\t" << average_identical_proportion_gap_compressed << "\t" << average_identical_proportion_blast << "\t" << average_identical_proportion_with_soft_clipping << "\t" <<
            "n_insertions" << "\t" << "insertion_length" << "\t" << average_insertion_gap_rate << "\t" << average_insertion_proportion_blast << "\t" <<
            "n_deletions" << "\t" << "deletion_length" << "\t" << average_deletion_gap_rate << "\t" << average_deletion_proportion_blast << "\t" <<
            "soft_clip_length" << "\t" << average_proportion_soft_clip <<
            std::endl;

    out.close();
    
    if (hts_close(in) < 0) {
        fprintf(stderr, "Error closing input.\n");
        exit_code = EXIT_FAILURE;
    }

    return exit_code;
}