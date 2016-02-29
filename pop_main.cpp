/** \file pop_main.cpp
 *  \brief Main entry point for evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.5
*/


#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "pop_main.hpp"
#include "pop_diverge.hpp"

int
main (int argc, char *argv[])
{
    if (argc < 2)
        {
            popbam_usage ();
            return 1;
        }
    char *user_func = NULL;
    user_func = strdup (argv[1]);
    argv++;
    argc--;
    if (strcmp (user_func, "snp") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "haplo") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "diverge") == 0)
        {
            main_diverge (argc, argv);
        }
    else if (strcmp (user_func, "tree") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "nucdiv") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "ld") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "sfs") == 0)
        {
            return 1;
        }
    else if (strcmp (user_func, "fasta") == 0)
        {
            return 1;
        }
    else
        {
            fprintf (stderr, "Error: unrecognized command: %s\n", user_func);
            popbam_usage ();
            return 1;
        }
    return 0;
}

int
popbam_usage (void)
{
    fputc ('\n', stderr);
    fputs ("Program: popbam\n", stderr);
    fputs ("Tools to perform evolutionary analysis from BAM files\n", stderr);
    fprintf (stderr, "Version: %s\n", version);
    fputs ("Usage: popbam <command> [options] <in.bam> [region]\n", stderr);
    fputc ('\n', stderr);
    fputs ("Commands:  snp       output consensus base calls\n", stderr);
    fputs ("           fasta     output alignment as multi-fasta file\n", stderr);
    fputs ("           haplo     output haplotype-based analyses\n", stderr);
    fputs ("           diverge   output divergence from reference\n", stderr);
    fputs ("           tree      output neighbor-joining trees\n", stderr);
    fputs ("           nucdiv    output nucleotide diversity statistics\n", stderr);
    fputs ("           ld        output linkage disequilibrium analysis\n", stderr);
    fputs ("           sfs       output site frequency spectrum analysis\n", stderr);
    fputc ('\n', stderr);
    return 0;
}
