/** \file pop_main.cpp
 *  \brief Main entry point for evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pop_main.hpp"
#include "pop_snp.hpp"
#include "pop_diverge.hpp"
#include "pop_haplo.hpp"
#include "pop_tree.hpp"
#include "pop_nucdiv.hpp"
#include "pop_ld.hpp"
#include "pop_sfs.hpp"
#include "pop_fasta.hpp"

int
main(int argc, char *argv[])
{
    if (argc < 2)
        {
            popbam_usage();
            return 1;
        }
    char *user_func = NULL;
    user_func = strdup(argv[1]);
    argv++;
    argc--;
    if (strcmp(user_func, "snp") == 0)
        {
            main_snp(argc, argv);
        }
    else if (strcmp(user_func, "haplo") == 0)
        {
            main_haplo(argc, argv);
        }
    else if (strcmp(user_func, "diverge") == 0)
        {
            main_diverge(argc, argv);
        }
    else if (strcmp(user_func, "tree") == 0)
        {
            main_tree(argc, argv);
        }
    else if (strcmp(user_func, "nucdiv") == 0)
        {
            main_nucdiv(argc, argv);
        }
    else if (strcmp(user_func, "ld") == 0)
        {
            main_ld(argc, argv);
        }
    else if (strcmp(user_func, "sfs") == 0)
        {
            main_sfs(argc, argv);
        }
    else if (strcmp(user_func, "fasta") == 0)
        {
            main_fasta(argc, argv);
        }
    else
        {
            fprintf (stderr, "Error: unrecognized command: %s\n", user_func);
            return 1;
        }
    return 0;
}

popbamData::popbamData(void)
{
    num_sites = 0;
    tid = -1;
    beg = 0;
    end = 0x7fffffff;
}

int
popbamData::assignPops(const popbamOptions *p)
{
    int i = 0;
    int si = -1;
    kstring_t buf;

    memset(&buf, 0, sizeof(kstring_t));
    for (i = 0; i < sm->n; i++)
        {
            if (sm->smpl[i])
                {
                    si = bam_smpl_sm2popid(sm, p->bamfile.c_str(), sm->smpl[i], &buf);
                }
            if (si < 0)
                {
                    si = bam_smpl_sm2popid(sm, p->bamfile.c_str(), 0, &buf);
                }
            if (si < 0)
                {
                    std::string msg;
                    std::string missing_sample(sm->smpl[i]);
                    msg = "Sample " + missing_sample +
                          " not assigned to a population.\nPlease check BAM header file definitions";
                    fatalError (msg);
                }
            pop_mask[si] |= 0x1ULL << i;
            pop_nsmpl[si]++;
        }
    return 0;
}

int
popbam_usage(void)
{
    fputchar('\n', stderr);
    fputs("Program: popbam\n", stderr);
    fputc("Tools to perform evolutionary analysis from BAM files\n", stderr);
    fprintf(stderr, "Version: %s\n", version);
    fputs("Usage: popbam <command> [options] <in.bam> [region]\n", stderr);
    fputchar('\n', stderr);
    fputs("Commands:  snp       output consensus base calls\n", stderr);
    fputs("           fasta     output alignment as multi-fasta file\n", stderr);
    fputs("           haplo     output haplotype-based analyses\n", stderr);
    fputs("           diverge   output divergence from reference\n", stderr);
    fputs("           tree      output neighbor-joining trees\n", stderr);
    fputs("           nucdiv    output nucleotide diversity statistics\n", stderr);
    fputs("           ld        output linkage disequilibrium analysis\n", stderr);
    fputs("           sfs       output site frequency spectrum analysis\n", stderr);
    fputchar('\n', stderr);
    return 0;
}
