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

#include "htslib/hts.h"
#include "htslib/hfile.h"

#include "pop_diverge.hpp"
#include "tables.h"

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
        fprintf (stderr, "read_file: detecting \"%s\" format failed.\n",
                 param.input_arg);
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
            fprintf (stderr, "htsfile: can't open %s: unknown format\n",
                     param.input_arg);
            exit (EXIT_FAILURE);
    }
    if (fp && (hclose(fp) < 0))
    {
        fprintf (stderr, "htsfile: closing %s failed\n", param.input_arg);
        exit (EXIT_FAILURE);
    }
    free (description);
    return 0;
}

int
main_diverge_bam (pop_diverge_parser *param)
{
    bool found = false;                  //! is the outgroup sequence found?
    int bp_indict = 0;                   //! bam parse indicator
    int i = 0;                           //! generic int iterator
    int chr = 0;                         //! chromosome identifier
    int beg = 0;                         //! beginning coordinate for analysis
    int end = 0;                         //! end coordinate for analysis
    int ref = 0;                         //! ref
    int win_size;                        //! calculated size of sliding window
    int num_windows = 0;                 //! number of windows
    bam_plbuf_t *buf = NULL;             //! pileup buffer
    diverge_data_bam *ddb;               //! diverge function data structure

    ddb = (diverge_data_bam*) malloc (sizeof(diverge_data_bam));

    // check input BAM file for errors
    check_BAM (ddb, param);

    // initialize the sample data structure
    ddb->sm = bam_smpl_init ();

    // add samples
    bam_smpl_add (ddb->sm, ddb->h, param->input_arg);

    // initialize the diverge data structre
    init_diverge_bam (ddb);

    // initialize error model
    ddb->em = errmod_init (0.17);

    // if outgroup option is used check to make sure it exists
    if (param->outgroup_given)
        {
            for (i = 0; i < ddb->sm->n; i++)
                {
                    if (strcmp (ddb->sm->smpl[i], param->outgroup_arg) == 0)
                        {
                            ddb->outidx = i;
                            found = true;
                        }
                }
            if (!found)
                {
                    fprintf (stderr, "Specified outgroup %s not found\n",
                             param->outgroup_arg);
                    exit (EXIT_FAILURE);
                }
        }

    // parse genomic region
    bp_indict = bam_parse_region (ddb->h, param->region_arg, &chr, &beg, &end);
    if (bp_indict < 0)
        {
            fprintf (stderr, "Bad genome coordinates: %s\n", param->region_arg);
            exit (EXIT_FAILURE);
        }

    // fetch reference sequence
    ddb->ref_base = faidx_fetch_seq (ddb->fai_file, ddb->h->target_name[chr],
                                     0, 0x7fffffff, &(ddb->len));

    // calculate the number of windows
    if (param->win_size_given)
        {
            win_size = (int)(param->win_size_arg * 1000);
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
            std::string scaffold_name(ddb->h->target_name[chr]);
            std::ostringstream winc(scaffold_name);
            winc.seekp(0, std::ios::end);
            winc << ':' << beg + (i * win_size) + 1 << '-';
            winc << ((i + 1) * win_size) + (beg - 1);
            std::string win_coord = winc.str();

            // initialize number of sites to zero
            ddb->num_sites = 0;

            // parse the BAM file and check if
            // region is retrieved from the reference
            if (param->win_size_given)
                {
                    int k = bam_parse_region (ddb->h, win_coord.c_str(), &ref, &(ddb->beg),
                                          &(ddb->end));
                    if (k < 0)
                        {
                            fprintf (stderr, "Bad window coordinates: %s\n",
                                     win_coord.c_str());
                            exit (EXIT_FAILURE);
                        }
                }
            else
                {
                    ref = chr;
                    ddb->beg = beg;
                    ddb->end = end;
                    if (ref < 0)
                        {
                            fprintf (stderr, "Bad scaffold name: %s\n",
                                     param->region_arg);
                            exit (EXIT_FAILURE);
                        }
                }

            // initialize diverge specific variables
            alloc_diverge_bam (ddb, param);

            // create population assignments
            assign_pops (ddb, param);

            // set default minimum sample size as
            // the number of samples in the population
            set_min_pop_n (ddb);

            // initialize pileup
            buf = bam_plbuf_init (make_diverge, param, ddb);

            // fetch region from bam file
            if ((bam_fetch (ddb->bam_in->x.bam, ddb->idx, ref, ddb->beg, ddb->end, buf, fetch_func)) < 0)
                {
                    fprintf (stderr, "Failed to retrieve region %s "
                             "due to corrupted BAM index file", param->region_arg);
                    exit (EXIT_FAILURE);
                }

            // finalize pileup
            bam_plbuf_push (0, buf);

            // print results to stdout
            calc_diverge (ddb, param);
            print_diverge_bam (ddb, param, ddb->h->target_name[chr], win_size);

            // take out the garbage
            bam_plbuf_destroy (buf);
        }  // end of window interation

    errmod_destroy (ddb->em);
    samclose (ddb->bam_in);
    bam_index_destroy (ddb->idx);
    bam_smpl_destroy (ddb->sm);
    free (ddb->ref_base);
    return 0;
}

int
main_diverge_vcf(pop_diverge_parser *param)
{
    diverge_data_vcf *ddv;               //! diverge function data structure
    init_diverge_vcf (ddv);
    return 0;
}

int
make_diverge (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl,
              void *par, void *data)
{
    int i = 0;
    int fq = 0;
    uint64_t sample_cov = 0;
    uint64_t *cb = NULL;
    diverge_data_bam *ddb = NULL;
    pop_diverge_parser *param = static_cast<pop_diverge_parser*>(par);

    // get control data structure
    ddb = (diverge_data_bam*)data;

    // only consider sites located in designated region
    if ((ddb->beg <= (int)pos) && (ddb->end > (int)pos))
        {
            // call bases
            cb = call_base_diverge (ddb, param, n, pl);

            // resolve heterozygous sites
<<<<<<< HEAD
            if (param->clean_hets)
=======
            //if (!(param->clean_hets))
            if (1)
>>>>>>> 38db4a4ab02cd25d0d234a800bb1a25b4d4aebac
                {
                    clean_hets (ddb->sm->n, cb, (int)ddb->ref_base[pos], param->min_snp_arg);
                }

            // determine if site is segregating
            fq = seg_base (ddb->sm->n, cb, ddb->ref_base[pos], param->min_snp_arg);

            // determine how many samples pass the quality filters
            sample_cov = qual_filter (ddb->sm->n, cb, param->min_rms_arg, param->min_depth_arg, param->max_depth_arg);
            for (i = 0; i < ddb->npops; i++)
                {
                    ddb->pop_sample_mask[i] = sample_cov & ddb->pop_mask[i];
                }

<<<<<<< HEAD
            if (bitcount64 (sample_cov) >= (int)(0.5 * ddb->sm->n))
=======
            if (bitcount64(sample_cov) == ddb->sm->n)
>>>>>>> 38db4a4ab02cd25d0d234a800bb1a25b4d4aebac
                {
                    // calculate the site type
                    ddb->types[ddb->num_sites] = calculate_sitetype (ddb->sm->n, cb);
                    if (fq > 0)
                        {
                            ddb->hap.pos[ddb->segsites] = pos;
                            ddb->hap.ref[ddb->segsites] = (uint8_t)bam_nt16_table[(int)ddb->ref_base[pos]];
                            for (i = 0; i < ddb->sm->n; i++)
                                {
                                    ddb->hap.rms[i][ddb->segsites] = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
                                    ddb->hap.snpq[i][ddb->segsites] = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;
                                    ddb->hap.num_reads[i][ddb->segsites] = (cb[i] >> (CHAR_BIT * 2)) & 0xffff;
                                    ddb->hap.base[i][ddb->segsites] = bam_nt16_table[(int)iupac[(cb[i] >> CHAR_BIT) & 0xff]];
                                    if (cb[i] & 0x2ULL)
                                        {
                                            ddb->hap.seq[i][ddb->segsites / 64] |= 0x1ULL << ddb->segsites % 64;
                                        }
                                }
                            ddb->hap.idx[ddb->segsites] = ddb->num_sites;
                            ddb->segsites++;
                        }
                    ddb->num_sites++;
                }

            // take out the garbage
            delete [] cb;
        }
    return 0;
}

int
calc_diverge (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    int i = 0;
    int j = 0;

    // calculate number of differences with reference sequence
    if (param->pop_flag)
        {
            uint16_t freq = 0;
            uint64_t pop_type = 0;
            for (i = 0; i < ddb->npops; i++)
                {
                    ddb->num_snps[i] = 0;
                    for (j = 0; j < ddb->segsites; j++)
                        {
                            pop_type = ddb->types[ddb->hap.idx[j]] & ddb->pop_mask[i];

                            // check if outgroup is different from reference
                            if (param->outgroup_given &&
                                CHECK_BIT(ddb->types[ddb->hap.idx[j]], ddb->outidx))
                                {
                                    freq = ddb->pop_nsmpl[i] - bitcount64(pop_type);
                                }
                            else
                                {
                                    freq = bitcount64 (pop_type);
                                }
                            if ((freq > 0) && (freq < ddb->pop_nsmpl[i]) &&
                                !param->single_flag)
                                {
                                    ++ddb->num_snps[i];
                                }
                            else if ((freq > 1) && (freq < ddb->pop_nsmpl[i]) &&
                                     param->single_flag)
                                {
                                    ++ddb->num_snps[i];
                                }
                            else if (freq == ddb->pop_nsmpl[i])
                                {
                                    ++ddb->pop_div[i];
                                }
                        }
                }
        }
    else
        {
            for (i = 0; i < ddb->sm->n; i++)
                {
                    for (j = 0; j <= SEG_IDX(ddb->segsites); j++)
                        {
                            ddb->ind_div[i] += bitcount64 (ddb->hap.seq[i][j]);
                        }
                }
        }
    return 0;
}

int
init_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    ddb->npops = ddb->sm->npops;
    return 0;
}

int
alloc_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    int i = 0;
    int length = ddb->end - ddb->beg;
    int npops = ddb->npops;
    int n = ddb->sm->n;

    ddb->segsites = 0;
    try
        {
            ddb->types = new uint64_t [length]();
            ddb->pop_mask = new uint64_t [npops]();
            ddb->pop_nsmpl = new uint8_t [npops]();
            ddb->pop_sample_mask = new uint64_t [npops]();
            ddb->min_pop_n = new uint16_t [npops]();
            ddb->num_snps = new int [npops]();
            ddb->hap.pos = new uint32_t [length]();
            ddb->hap.idx = new uint32_t [length]();
            ddb->hap.ref = new uint8_t [length]();
            ddb->hap.seq = new uint64_t* [n];
            ddb->hap.base = new uint8_t* [n];
            ddb->hap.rms = new uint16_t* [n];
            ddb->hap.snpq = new uint16_t* [n];
            ddb->hap.num_reads = new uint16_t* [n];
            if (param->pop_flag)
                {
                    ddb->pop_div = new uint16_t [npops]();
                }
            else
                {
                    ddb->ind_div = new uint16_t [n]();
                }
            for (i = 0; i < n; i++)
                {
                    ddb->hap.seq[i] = new uint64_t [length]();
                    ddb->hap.base[i] = new uint8_t [length]();
                    ddb->hap.rms[i] = new uint16_t [length]();
                    ddb->hap.snpq[i] = new uint16_t [length]();
                    ddb->hap.num_reads[i] = new uint16_t [length]();
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    return 0;
}

int
init_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param)
{
    return 0;
}

int
alloc_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param)
{
    return 0;
}

void dealloc_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    int i = 0;
    int n = ddb->sm->n;

    delete [] ddb->pop_mask;
    delete [] ddb->types;
    delete [] ddb->pop_nsmpl;
    delete [] ddb->pop_sample_mask;
    delete [] ddb->num_snps;
    delete [] ddb->min_pop_n;
    delete [] ddb->hap.pos;
    delete [] ddb->hap.idx;
    delete [] ddb->hap.ref;
    if (param->pop_flag)
        {
            delete [] ddb->pop_div;
        }
    else
        {
            delete [] ddb->ind_div;
        }
    for (i = 0; i < n; i++)
        {
            delete [] ddb->hap.seq[i];
            delete [] ddb->hap.base[i];
            delete [] ddb->hap.num_reads[i];
            delete [] ddb->hap.snpq[i];
            delete [] ddb->hap.rms[i];
        }
    delete [] ddb->hap.seq;
    delete [] ddb->hap.base;
    delete [] ddb->hap.snpq;
    delete [] ddb->hap.rms;
    delete [] ddb->hap.num_reads;
}

void
dealloc_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param)
{
    return;
}

int
print_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param,
                   const char *scaffold, int win_size)
{
    int i = 0;
    double pdist = 0.0;
    double jc = 0.0;
    std::stringstream out;

    out << scaffold << '\t' << ddb->beg + 1 << '\t' << ddb->end + 1 << '\t' << ddb->num_sites;
    if (param->pop_flag)
        {
            for (i = 0; i < ddb->npops; i++)
                {
                    if (ddb->num_sites >= (int)(param->min_sites_arg * win_size))
                        {
                            out << "\tFixed[" << ddb->sm->popul[i] << "]:\t" << ddb->pop_div[i];
                            out << "\tSeg[" << ddb->sm->popul[i] << "]:\t" << ddb->num_snps[i];
                            out << "\td[" << ddb->sm->popul[i] << "]:";
                            if (strcmp (param->dist_arg, "pdist") == 0)
                                {
                                    if (param->subst_flag)
                                        {
                                            out << '\t' << std::fixed << std::setprecision(5) <<
                                            (double)(ddb->pop_div[i]) / ddb->num_sites;
                                        }
                                    else
                                        {
                                            out << '\t' << std::fixed << std::setprecision(5) <<
                                            (double)(ddb->pop_div[i] + ddb->num_snps[i]) / ddb->num_sites;
                                        }
                                }
                            else if (strcmp (param->dist_arg, "jc") == 0)
                                {
                                    if (param->subst_flag)
                                        {
                                            pdist = (double)(ddb->pop_div[i]) / ddb->num_sites;
                                        }
                                    else
                                        {
                                            pdist = (double)(ddb->pop_div[i] + ddb->num_snps[i]) / ddb->num_sites;
                                        }
                                    jc = -0.75 * log (1.0 - pdist * (4.0 / 3.0));
                                    out << '\t' << std::fixed << std::setprecision(5) << jc;
                                }
                            else
                                {
                                    out << "\tFixed[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                    out << "\tSeg[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                    out << "\td[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                                }
                        }
                    else
                        {
                            out << "\tFixed[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                            out << "\tSeg[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                            out << "\td[" << ddb->sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                        }
                }
        }
    else
        {
            for (i = 0; i < ddb->sm->n; i++)
                {
                    if (ddb->num_sites >= (int)(param->min_sites_arg * win_size))
                        {
                            if (strcmp (param->dist_arg, "pdist") == 0)
                                {
                                    out << "\td[" << ddb->sm->smpl[i] << "]:";
                                    out << '\t' << std::fixed << std::setprecision(5) <<
                                            (double)(ddb->ind_div[i]) / ddb->num_sites;
                                }
                            else if (strcmp (param->dist_arg, "jc") == 0)
                                {
                                    pdist = (double)(ddb->ind_div[i]) / ddb->num_sites;
                                    jc = -0.75 * log (1.0 - pdist * (4.0 / 3.0));
                                    out << "\td[" << ddb->sm->smpl[i] << "]:";
                                    out << '\t' << std::fixed << std::setprecision(5) << jc;
                                }
                            else
                                {
                                    out << "\td[" << ddb->sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
                                }
                        }
                    else
                        {
                            out << "\td[" << ddb->sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
                        }
                }

        }
    std::cout << out.str() << std::endl;
    return 0;
}

int
print_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param,
                   const char *scaffold, int win_size)
{
    return 0;
}

// TODO: The set_min_pop_n function cannot currently take input from user

int
set_min_pop_n (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    int i = 0;

    for (i = 0; i < ddb->npops; i++)
        {
            ddb->min_pop_n[i] = (uint16_t)(0.5 * ddb->pop_nsmpl[i]);
        }
    return 0;
}

int
assign_pops (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    int i = 0;
    int si = -1;
    kstring_t buf;

    memset (&buf, 0, sizeof(kstring_t));
    for (i = 0; i < ddb->sm->n; i++)
        {
            if (ddb->sm->smpl[i])
                {
                    si = bam_smpl_sm2popid (ddb->sm, param->input_arg, ddb->sm->smpl[i], &buf);
                }
            if (si < 0)
                {
                    si = bam_smpl_sm2popid (ddb->sm, param->input_arg, 0, &buf);
                }
            if (si < 0)
                {
                    fprintf (stderr, "Sample %s not assigned to a population.\n"
                             "Please check BAM header file definitions.", ddb->sm->smpl[i]);
                    exit (EXIT_FAILURE);
                }
            ddb->pop_mask[si] |= 0x1ULL << i;
            ddb->pop_nsmpl[si]++;
        }
    return 0;
}

uint64_t *
call_base_diverge (diverge_data_bam *ddb, const pop_diverge_parser *param, int n,
                   const bam_pileup1_t *pl)
{
    int i = 0;
    int j = 0;
    int qq = 0;
    int base_q = 0;
    int tmp_base_q = 0;
    int b = 0;
    int si = -1;
    int rmsq = 0;
    int n_smpl = ddb->sm->n;
    uint16_t k = 0;
    uint16_t *bases = NULL;
    uint64_t rms = 0;
    uint64_t *cb;
    std::vector<int> depth(n_smpl, 0);
    uint8_t *s = NULL;
    float q[16];
    std::string msg;
    bam_pileup1_t ***p = NULL;
    kstring_t buf;

    // allocate memory pileup data
    try
        {
            cb = new uint64_t [n_smpl]();
            p = new bam_pileup1_t** [n_smpl];
            for (i = 0; i < n_smpl; i++)
                {
                    p[i] = new bam_pileup1_t* [param->max_depth_arg];
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    memset (&buf, 0, sizeof(kstring_t));

    // partition pileup according to sample
    for (i = 0; i < n; i++)
        {
            if ((pl+i)->is_del || (pl+i)->is_refskip || ((pl+i)->b->core.flag & BAM_FUNMAP))
                {
                    continue;
                }
            s = bam_aux_get ((pl+i)->b, "RG");

            // skip reads with no read group tag
            if (!s)
                {
                    continue;
                }
            else
                {
                    si = bam_smpl_rg2smid (ddb->sm, param->input_arg, (char*)(s+1), &buf);
                }
            if (si < 0)
                {
                    si = bam_smpl_rg2smid (ddb->sm, param->input_arg, 0, &buf);
                }
            if (si < 0)
                {
                    free (buf.s);
                    char *rogue_rg = bam_aux2Z(s);
                    fprintf (stderr, "Problem assigning read group %s to a sample.\n"
                             "Please check BAM header for correct SM and PO tags",
                             rogue_rg);
                    exit (EXIT_FAILURE);
                }
            if (depth[si] < param->max_depth_arg)
                {
                    p[si][depth[si]] = const_cast<bam_pileup1_t*>(pl+i);
                    depth[si]++;
                }
            else
                {
                    continue;
                }
        }

    // fill in the base array
    for (j = 0; j < n_smpl; ++j)
        {
            rmsq = 0;
            if (depth[j] > 0)
                {
                    try
                        {
                            bases = new uint16_t [depth[j]]();
                        }
                    catch (std::bad_alloc& ba)
                        {
                            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
                        }

                    for (i = k = 0; i < depth[j]; ++i)
                        {
                            tmp_base_q = bam1_qual(p[j][i]->b)[p[j][i]->qpos];

                            if (param->illumina_flag)
                                {
                                    base_q = tmp_base_q > 31 ? tmp_base_q - 31 : 0;
                                }
                            else
                                {
                                    base_q = tmp_base_q;
                                }
                            if (base_q < 0)
                                {
                                    fputs ("Scaling problem for base scores\n", stderr);
                                    exit (EXIT_FAILURE);
                                }
                            if ((base_q < param->min_base_arg) ||
                                (p[j][i]->b->core.qual < param->min_map_arg))
                                {
                                    continue;
                                }
                            b = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p[j][i]->b), p[j][i]->qpos)];
                            if (b > 3)
                                {
                                    continue;
                                }
                            qq = base_q < p[j][i]->b->core.qual ? base_q : p[j][i]->b->core.qual;
                            if (qq < 4)
                                {
                                    qq = 4;
                                }
                            if (qq > 63)
                                {
                                    qq = 63;
                                }
                            bases[k++] = qq << 5 | (uint16_t)bam1_strand(p[j][i]->b) << 4 | b;
                            rmsq += p[j][i]->b->core.qual * param->max_depth_arg;
                        }

                    // calculate genotype likelihoods
                    errmod_cal (ddb->em, k, NBASES, bases, q);

                    // finalize root mean quality score
                    rms = (uint64_t)(sqrt((float)(rmsq) / k) + 0.499);

                    // get consensus base call
                    cb[j] = gl2cns(q, k);

                    // add root-mean map quality score to cb array
                    cb[j] |= rms << (CHAR_BIT * 6);

                    // take out some garbage
                    delete [] bases;
                    bases = NULL;
                }
            else
                {
                    continue;
                }
        }

    // take out garbage
    for (i = 0; i < n_smpl; i++)
        {
            delete [] p[i];
        }
    delete [] p;
    free (buf.s);

    return cb;
}

int
check_BAM (diverge_data_bam *ddb, const pop_diverge_parser *param)
{
    ddb->bam_in = samopen (param->input_arg, "rb", 0);

    // check if BAM file is readable
    if (!ddb->bam_in)
        {
            fprintf (stderr, "Cannot read infile \"%s\"\n", param->input_arg);
            exit (EXIT_FAILURE);
        }

    // check if BAM header is returned
    if (!ddb->bam_in->header)
        {
            fprintf (stderr, "Cannot read BAM header from file \"%s\"\n", param->input_arg);
            exit (EXIT_FAILURE);
        }
    else
        {
            ddb->h = ddb->bam_in->header;
        }

    // read in new header text
    //if (flag & BAM_HEADERIN)
    //    {
    //        std::ifstream headin(headfile);
    //        headin.seekg(0, std::ios::end);
    //        h->l_text = headin.tellg();
    //        headin.seekg(0, std::ios::beg);
    //        h->text = (char*)realloc(h->text, (size_t)h->l_text);
    //        headin.read(h->text, h->l_text);
    //        headin.close();
    //    }

    // check for bam index file
    if (!(ddb->idx = bam_index_load (param->input_arg)))
        {
            fprintf (stderr, "Index file not available for BAM file \"%s\"\n", param->input_arg);
            exit (EXIT_FAILURE);
        }

    // check if fastA reference index is available
    ddb->fai_file = fai_load (param->ref_arg);
    if (!ddb->fai_file)
        {
            fprintf (stderr, "Failed to load index for fastA reference file \"%s\"\n", param->ref_arg);
            exit (EXIT_FAILURE);
        }
    return 0;
}
