/** \file pop_diverge.cpp
 *  \brief Functions for calculating divergence from reference genome
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <new>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"

#include "pop_global.hpp"
#include "pop_options.hpp"
#include "pop_sample.hpp"
#include "pop_utils.hpp"
#include "pop_main.hpp"
#include "pop_diverge_parser.hpp"
#include "pop_diverge.hpp"

template uint64_t* callBase<divergeData>(divergeData *t, int n, const bam_pileup1_t *pl);

int
main_diverge (int argc, char *argv[])
{
    char *description = NULL;              //! description of input file format
    pop_diverge_parser param(argc, argv);  //! parse command line options
    htsFormat fmt;                         //! input file format
    hFILE *fp;                             //! input file handle

    fp = hopen (param.input_arg, "r");
    if (fp == NULL)
    {
        fprintf (stderr, "read_file: can't open \"%s\"\n", param.input_arg);
        exit (EXIT_FAILURE);
    }
    if (hts_detect_format (fp, &fmt) < 0)
    {
        fprintf(stderr, "read_file: detecting \"%s\" format failed.\n", param.input_arg);
        hclose_abruptly (fp);
        exit (EXIT_FAILURE);
    }
    description = hts_format_description (&fmt);
    fprintf (stdout, "%s:\t%s\n", param.input_arg, description);
    switch (fmt.category)
    {
        case sequence_data:
            main_diverge_bam (&param);
            break;
        case variant_data:
            main_diverge_vcf (&param);
            break;
        default:
            fprintf (stderr, "htsfile: can't open %s: unknown format\n", param.input_arg);
            exit (EXIT_FAILURE);
    }
    if (fp && (hclose(fp) < 0))
    {
        fprintf(stderr, "htsfile: closing %s failed\n", param.input_arg);
        exit(EXIT_FAILURE);
    }
    free (description);
    return 0;
}

int
main_diverge_bam(pop_diverge_parser *param)
{
    bool found = false;                  //! is the outgroup sequence found?
    int bp_indict = 0;                   //! bam parse indicator
    int i = 0;
    int chr = 0;                         //! chromosome identifier
    int beg = 0;                         //! beginning coordinate for analysis
    int end = 0;                         //! end coordinate for analysis
    int ref = 0;                         //! ref
    int win_size;                        //! calculated size of sliding window
    long num_windows = 0;                //! number of windows
    std::string msg;                     //! string for error message
    bam_sample_t *sm = NULL;             //! pointer to the sample information for the input BAM file
    bam_plbuf_t *buf = NULL;             //! pileup buffer
    samfile_t *bam_in;                   //! handle for input BAM file
    bam_header_t *h;                     //! pointer to BAM header
    bam_index_t *idx;                    //! pointer to BAM index
    faidx_t *fai_file;                   //! handle for indexed reference fastA file
    errmod_t *em;                        //! pointer to error model structure

    // check input BAM file for errors
    param->checkBAM (bam_in, h, idx, fai_file);

    // initialize the sample data structure
    sm = bam_smpl_init ();

    // add samples
    bam_smpl_add (sm, h, param->input_arg);

    // initialize the diverge data structre
    divergeData t(arg);
    t.npops = sm->npops;

    // initialize error model
    em = errmod_init (0.17);

    // if outgroup option is used check to make sure it exists
    if (param->outgroup_given)
        {
            for (i = 0; i < sm->n; i++)
                {
                    if (strcmp (sm->smpl[i], param->outgroup_arg) == 0)
                        {
                            t.outidx = i;
                            found = true;
                        }
                }
            if (!found)
                {
                    fprintf (stderr, "Specified outgroup %s not found\n", param->outgroup_arg);
                    exit (EXIT_FAILURE);
                }
        }

    // parse genomic region
    bp_indict = bam_parse_region (h, param->region_arg, &chr, &beg, &end);
    if (bp_indict < 0)
        {
            fprintf (stderr, "Bad genome coordinates: %s\n", param->region_arg);
            exit (EXIT_FAILURE);
        }

    // fetch reference sequence
    t.ref_base = faidx_fetch_seq (param->ref_arg, h->target_name[chr], 0, 0x7fffffff, &(t.len));

    // calculate the number of windows
    if (param->win_size_given)
        {
            win_size = (int)(param->win_size_arg * 1000)
            num_windows = ((end - beg) - 1) / win_size;
        }
    else
        {
            win_size = end - beg;
            num_windows = 1;
        }

    // iterate through all windows along specified genomic region
    for (i = 0; i < num_windows; i++)
        {
            // construct genome coordinate string
            std::string scaffold_name(h->target_name[chr]);
            std::ostringstream winc(scaffold_name);
            winc.seekp(0, std::ios::end);
            winc << ':' << beg + (i * win_size) + 1 << '-' << ((i + 1) * win_size) + (beg - 1);
            std::string win_coord = winc.str();

            // initialize number of sites to zero
            t.num_sites = 0;

            // parse the BAM file and check if region is retrieved from the reference
            if (param->win_size_given)
                {
                    k = bam_parse_region (h, win_coord, &ref, &(t.beg), &(t.end));
                    if (k < 0)
                        {
                            fprintf (stderr, "Bad window coordinates: %s\n", win_coord);
                            exit (EXIT_FAILURE);
                        }
                }
            else
                {
                    ref = chr;
                    t.beg = beg;
                    t.end = end;
                    if (ref < 0)
                        {
                            fprintf (stderr, "Bad scaffold name: %s\n", param->region_arg);
                            exit (EXIT_FAILURE);
                        }
                }

            // initialize diverge specific variables
            t.alloc_diverge();

            // create population assignments
            t.assign_pops(arg.input_arg);

            // set default minimum sample size as
            // the number of samples in the population
            t.set_min_pop_n();

            // initialize pileup
            buf = bam_plbuf_init (make_diverge, &t);

            // fetch region from bam file
            if ((bam_fetch (bam_in->x.bam, idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
                {
                    fprintf (stderr, "Failed to retrieve region %s "
                             "due to corrupted BAM index file", param->region_arg);
                    exit (EXIT_FAILURE);
                }

            // finalize pileup
            bam_plbuf_push (0, buf);

            // print results to stdout
            t.calc_diverge (param);
            t.print_diverge (param, h->target_name[chr], win_size);

            // take out the garbage
            bam_plbuf_destroy (buf);
        }  // end of window interation

    errmod_destroy (em);
    samclose (bam_in);
    bam_index_destroy (idx);
    bam_smpl_destroy (sm);
    free (t.ref_base);
    return 0;
}

int
main_diverge_vcf(pop_diverge_parser *param)
{

}

int
make_diverge(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i = 0;
    int fq = 0;
    uint64_t sample_cov = 0;
    uint64_t *cb = NULL;
    divergeData *t = NULL;

    // get control data structure
    t = (divergeData*)data;

    // only consider sites located in designated region
    if ((t->beg <= (int)pos) && (t->end > (int)pos))
        {
            // call bases
            cb = callBase(t, n, pl);

            // resolve heterozygous sites
            if (!(t->flag & BAM_HETEROZYGOTE))
                {
                    cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t.minSNPQ);
                }

            // determine if site is segregating
            fq = segBase(t->sm->n, cb, t->ref_base[pos], t->minSNPQ);

            // determine how many samples pass the quality filters
            sample_cov = qualFilter(t->sm->n, cb, t->minRMSQ, t->minDepth, t->maxDepth);

            for (i = 0; i < t->npops; i++)
                {
                    t->pop_sample_mask[i] = sample_cov & t->pop_mask[i];
                }

            if (bitcount64(sample_cov) == t->sm->n)
                {
                    // calculate the site type
                    t->types[t->num_sites] = calculateSiteType(t->sm->n, cb);

                    if (fq > 0)
                        {
                            t->hap.pos[t->segsites] = pos;
                            t->hap.ref[t->segsites] = (uint8_t)bam_nt16_table[(int)t->ref_base[pos]];
                            for (i = 0; i < t->sm->n; i++)
                                {
                                    t->hap.rms[i][t->segsites] = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
                                    t->hap.snpq[i][t->segsites] = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;
                                    t->hap.num_reads[i][t->segsites] = (cb[i] >> (CHAR_BIT * 2)) & 0xffff;
                                    t->hap.base[i][t->segsites] = bam_nt16_table[(int)iupac[(cb[i] >> CHAR_BIT) & 0xff]];
                                    if (cb[i] & 0x2ULL)
                                        {
                                            t->hap.seq[i][t->segsites/64] |= 0x1ULL << t->segsites % 64;
                                        }
                                }
                            t->hap.idx[t->segsites] = t->num_sites;
                            t->segsites++;
                        }
                    t->num_sites++;
                }

            // take out the garbage
            delete [] cb;
        }
    return 0;
}

int
divergeData::calc_diverge (const pop_diverge_parser *param)
{
    int i = 0;
    int j = 0;
    uint16_t freq = 0;
    uint64_t pop_type = 0;

    // calculate number of differences with reference sequence
    switch (out_format)
        {
        case 0:
            for (i = 0; i < sm->n; i++)
                {
                    for (j = 0; j <= SEG_IDX(segsites); j++)
                        {
                            ind_div[i] += bitcount64(hap.seq[i][j]);
                        }
                }
            break;
        case 1:
            for (i = 0; i < npops; i++)
                {
                    num_snps[i] = 0;
                    for (j = 0; j < segsites; j++)
                        {
                            pop_type = types[hap.idx[j]] & pop_mask[i];

                            // check if outgroup is different from reference
                            if (pdp->outgroup_given && CHECK_BIT(types[hap.idx[j]], outidx))
                                {
                                    freq = pop_nsmpl[i] - bitcount64(pop_type);
                                }
                            else
                                {
                                    freq = bitcount64(pop_type);
                                }
                            if ((freq > 0) && (freq < pop_nsmpl[i]) && !pdp->single_flag)
                                {
                                    ++num_snps[i];
                                }
                            else if ((freq > 1) && (freq < pop_nsmpl[i]) && pdp->single_flag)
                                {
                                    ++num_snps[i];
                                }
                            else if (freq == pop_nsmpl[i])
                                {
                                    ++pop_div[i];
                                }
                        }
                }
            break;
        default:
            break;
        }
    return 0;
}

divergeData::divergeData(void)
{
    derived_type = DIVERGE;
}

int
divergeData::alloc_diverge(void)
{
    int i = 0;
    int length = end - beg;
    int n = sm->n;

    segsites = 0;
    try
        {
            types = new uint64_t [length]();
            pop_mask = new uint64_t [npops]();
            pop_nsmpl = new uint8_t [npops]();
            pop_sample_mask = new uint64_t [npops]();
            min_pop_n = new uint16_t [npops]();
            num_snps = new int [npops]();
            hap.pos = new uint32_t [length]();
            hap.idx = new uint32_t [length]();
            hap.ref = new uint8_t [length]();
            hap.seq = new uint64_t* [n];
            hap.base = new uint8_t* [n];
            hap.rms = new uint16_t* [n];
            hap.snpq = new uint16_t* [n];
            hap.num_reads = new uint16_t* [n];
            switch (out_format)
                {
                case 0:
                    ind_div = new uint16_t [n]();
                    break;
                case 1:
                    pop_div = new uint16_t [npops]();
                default:
                    break;
                }
            for (i = 0; i < n; i++)
                {
                    hap.seq[i] = new uint64_t [length]();
                    hap.base[i] = new uint8_t [length]();
                    hap.rms[i] = new uint16_t [length]();
                    hap.snpq[i] = new uint16_t [length]();
                    hap.num_reads[i] = new uint16_t [length]();
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    return 0;
}

divergeData::~divergeData(void)
{
    int i = 0;
    int n = sm->n;

    delete [] pop_mask;
    delete [] types;
    delete [] pop_nsmpl;
    delete [] pop_sample_mask;
    delete [] num_snps;
    delete [] min_pop_n;
    delete [] hap.pos;
    delete [] hap.idx;
    delete [] hap.ref;
    switch (out_format)
        {
        case 0:
            delete [] ind_div;
            break;
        case 1:
            delete [] pop_div;
            break;
        default:
            break;
        }
    for (i = 0; i < n; i++)
        {
            delete [] hap.seq[i];
            delete [] hap.base[i];
            delete [] hap.num_reads[i];
            delete [] hap.snpq[i];
            delete [] hap.rms[i];
        }
    delete [] hap.seq;
    delete [] hap.base;
    delete [] hap.snpq;
    delete [] hap.rms;
    delete [] hap.num_reads;
}

int
divergeData::print_diverge(const pop_diverge_parser *param, const char *scaffold, int win_size)
{
    int i = 0;
    double pdist = 0.0;
    double jc = 0.0;
    std::stringstream out;

    out << scaffold << '\t' << beg + 1 << '\t' << end + 1 << '\t' << num_sites;
    switch (param->format_arg)
        {
        case 0:
            for (i = 0; i < sm->n; i++)
                {
                    if (num_sites >= (int)(param->min_sites_arg * win_size))
                        {
                            if (strcmp(param->dist_arg, "pdist") == 0)
                                {
                                    out << "\td[" << sm->smpl[i] << "]:";
                                    out << '\t' << std::fixed << std::setprecision(5) << 
                                            (double)(ind_div[i]) / num_sites;
                                }
                            else if (strcmp(param->dist_arg, "jc") == 0)
                                {
                                    pdist = (double)(ind_div[i]) / num_sites;
                                    jc = -0.75 * log(1.0 - pdist * (4.0 / 3.0));
                                    out << "\td[" << sm->smpl[i] << "]:";
                                    out << '\t' << std::fixed << std::setprecision(5) << jc;
                                }
                            else
                                {
                                    out << "\td[" << sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
                                }
                        }
                    else
                        {
                            out << "\td[" << sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
                        }
                }
            break;
        case 1:
            for (i = 0; i < npops; i++)
                {
                    if (num_sites >= minSites)
                        {
                            out << "\tFixed[" << sm->popul[i] << "]:\t" << pop_div[i];
                            out << "\tSeg[" << sm->popul[i] << "]:\t" << num_snps[i];
                            out << "\td[" << sm->popul[i] << "]:";
                            if (strcmp(param->dist_arg, "pdist") == 0)
                                {
                                    if (param->subst_flag)
                                        {
                                            out << '\t' << std::fixed << std::setprecision(5) << 
                                                    (double)(pop_div[i]) / num_sites;
                                        }
                                    else
                                        {
                                            out << '\t' << std::fixed << std::setprecision(5) << 
                                                    (double)(pop_div[i] + num_snps[i]) / num_sites;
                                        }
                                }
                            else if (strcmp(param->dist_arg, "jc") == 0)
                                {
                                    if (param->subst_flag)
                                        {
                                            pdist = (double)(pop_div[i]) / num_sites;
                                        }
                                    else
                                        {
                                            pdist = (double)(pop_div[i] + num_snps[i]) / num_sites;
                                        }
                                    jc = -0.75 * log(1.0 - pdist * (4.0 / 3.0));
                                    out << '\t' << std::fixed << std::setprecision(5) << jc;
                                }
                            else
                                {
                                    out << "\tFixed[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                    out << "\tSeg[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                    out << "\td[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                }
                        }
                    else
                        {
                            out << "\tFixed[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                            out << "\tSeg[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                            out << "\td[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                        }
                }
            break;
        default:
            break;
        }
    std::cout << out.str() << std::endl;
    return 0;
}

// TODO: The set_min_pop_n function cannot currently take input from user

int
divergeData::set_min_pop_n(void)
{
    int i = 0;

    for (i = 0; i < npops; i++)
        {
            min_pop_n[i] = (unsigned short)pop_nsmpl[i];
        }
    return 0;
}

