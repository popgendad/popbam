/** \file pop_sfs.cpp
 *  \brief Functions for calculating frequency spectrum statistics
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include <cstdlib>
#include <cstdint>
#include <string>
#include <limits>
#include <iostream>
#include <new>
#include <vector>
#include "bam.h"
#include "faidx.h"
#include "sam.h"
#include "kstring.h"
#include "pop_sample.h"
#include "pop_utils.h"
#include "popbam.h"
#include "pop_sfs.h"

template uint64_t* callBase<sfsData>(sfsData *t, int n, const bam_pileup1_t *pl);
int makeSFS(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
void usageSFS(const std::string);

int
mainSFS(int argc, char *argv[])
{
    bool found = false;           //! is the outgroup sequence found?
    int chr = 0;                  //! chromosome identifier
    int beg = 0;                  //! beginning coordinate for analysis
    int end = 0;                  //! end coordinate for analysis
    int ref = 0;                  //! ref
    long nWindows = 0;            //! number of windows
    std::string msg;              //! string for error message
    bam_sample_t *sm = nullptr;   //! Pointer to the sample information for the input BAM file
    bam_plbuf_t *buf = nullptr;   //! pileup buffer

    // initialize user command line options
    popbamOptions p(argc, argv);

    if (p.errorCount > 0)
        {
            usageSFS(p.errorMsg);
        }

    // check input BAM file for errors
    p.checkBAM();

    // initialize the sample data structure
    sm = bam_smpl_init();

    // add samples
    bam_smpl_add(sm, &p);

    // initialize the sfs data structre
    sfsData t(p);
    t.sm = sm;

    // initialize error model
    t.em = errmod_init(0.17);

    // if outgroup option is used check to make sure it exists
    if (p.flag & BAM_OUTGROUP)
        {
            for (int i = 0; i < t.sm->n; ++i)
                {
                    if (strcmp(t.sm->smpl[i], t.outgroup.c_str()) == 0)
                        {
                            t.outidx = i;
                            found = true;
                        }
                }
            if (!found)
                {
                    msg = "Specified outgroup " + t.outgroup + " not found";
                    fatalError(msg);
                }
        }

    // calculate the constants for computation of Tajima's D
    t.calc_a1();
    t.calc_a2();
    t.calc_e1();
    t.calc_e2();
    t.calc_dw();
    t.calc_hw();

    // parse genomic region
    int k = bam_parse_region(p.h, p.region, &chr, &beg, &end);
    if (k < 0)
        {
            msg = "Bad genome coordinates: " + p.region;
            fatalError(msg);
        }

    // fetch reference sequence
    t.ref_base = faidx_fetch_seq(p.fai_file, p.h->target_name[chr], 0, 0x7fffffff, &(t.len));

    // calculate the number of windows
    if (p.flag & BAM_WINDOW)
        {
            nWindows = ((end - beg) - 1) / p.winSize;
        }
    else
        {
            p.winSize = end - beg;
            nWindows = 1;
        }

    // iterate through all windows along specified genomic region
    for (long j = 0; j < nWindows; ++j)
        {
            // construct genome coordinate string
            std::string scaffold_name(p.h->target_name[chr]);
            std::ostringstream winc(scaffold_name);
            winc.seekp(0, std::ios::end);
            winc << ':' << beg + (j * p.winSize) + 1 << '-' << ((j + 1) * p.winSize) + (beg - 1);
            std::string winCoord = winc.str();

            // initialize number of sites to zero
            t.num_sites = 0;

            // parse the BAM file and check if region is retrieved from the reference
            if (p.flag & BAM_WINDOW)
                {
                    k = bam_parse_region(p.h, winCoord, &ref, &(t.beg), &(t.end));
                    if (k < 0)
                        {
                            msg = "Bad window coordinates " + winCoord;
                            fatalError(msg);
                        }
                }
            else
                {
                    ref = chr;
                    t.beg = beg;
                    t.end = end;
                    if (ref < 0)
                        {
                            msg = "Bad scaffold name: " + p.region;
                            fatalError(msg);
                        }
                }

            // initialize nucdiv variables
            t.allocSFS();

            // create population assignments
            t.assignPops(&p);

            // assign outgroup population
            if ((p.flag & BAM_OUTGROUP) && found)
                {
                    t.assignOutpop();
                }

            // initialize pileup
            buf = bam_plbuf_init(makeSFS, &t);

            // fetch region from bam file
            if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
                {
                    msg = "Failed to retrieve region " + p.region + " due to corrupted BAM index file";
                    fatalError(msg);
                }

            // finalize pileup
            bam_plbuf_push(0, buf);

            // calculate site frequency spectrum statistics
            t.calcSFS();

            // print results to stdout
            t.printSFS(std::string(p.h->target_name[chr]));

            // take out the garbage
            bam_plbuf_destroy(buf);
        }   // end of window interation

    errmod_destroy(t.em);
    samclose(p.bam_in);
    bam_index_destroy(p.idx);
    for (int i = 0; i <= t.sm->n; ++i)
        {
            delete [] t.dw[i];
            delete [] t.hw[i];
        }
    delete [] t.dw;
    delete [] t.hw;
    delete [] t.a1;
    delete [] t.a2;
    delete [] t.e1;
    delete [] t.e2;
    free(t.ref_base);
    return 0;
}

int
makeSFS(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i = 0;
    int j = 0;
    int fq = 0;
    uint64_t sample_cov = 0;
    uint64_t *cb = nullptr;
    sfsData *t = nullptr;

    // get control data structure
    t = (sfsData*)data;

    // only consider sites located in designated region
    if ((t->beg <= (int)pos) && (t->end > (int)pos))
        {
            // call bases
            cb = callBase(t, n, pl);

            // resolve heterozygous sites
            if (!(t->flag & BAM_HETEROZYGOTE))
                {
                    cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->minSNPQ);
                }

            // determine if site is segregating
            fq = segBase(t->sm->n, cb, t->ref_base[pos], t->minSNPQ);

            // determine how many samples pass the quality filters
            sample_cov = qualFilter(t->sm->n, cb, t->minRMSQ, t->minDepth, t->maxDepth);

            uint32_t *ncov = nullptr;
            ncov = new uint32_t [t->npops]();

            // determine population coverage
            for (i = 0; i < t->npops; ++i)
                {
                    uint64_t pc = 0;
                    pc = sample_cov & t->pop_mask[i];
                    ncov[i] = bitcount64(pc);
                    uint32_t req = (uint32_t)((t->minPop * t->pop_nsmpl[i]) + 0.4999);
                    if (ncov[i] >= req)
                        {
                            t->pop_cov[t->num_sites] |= 0x1U << i;
                        }
                }

            // record site type if the site is variable
            if (t->pop_cov[t->num_sites] > 0)
                {
                    t->num_sites++;
                    if (fq > 0)
                        {
                            for (j = 0; j < t->npops; ++j)
                                {
                                    t->ncov[j][t->segsites] = ncov[j];
                                }
                            t->types[t->segsites++] = calculateSiteType(t->sm->n, cb);
                        }
                }

            // take out the garbage
            delete [] cb;
            delete [] ncov;
        }
    return 0;
}

int
sfsData::calcSFS(void)
{
    int i = 0;
    int j = 0;
    int n = 0;
    int s = 0;
    int avgn = 0;
    uint16_t freq = 0;
    uint64_t pop_type = 0;

    // count number of aligned sites in each population
    for (i = 0; i < num_sites; ++i)
        {
            for (j = 0; j < npops; ++j)
                {
                    if (CHECK_BIT(pop_cov[i],j))
                        {
                            ++ns[j];
                        }
                }
        }
    for (i = 0; i < npops; i++)
        {
            if (ns[i] >= (uint32_t)((end - beg) * minSites))
                {
                    // get site frequency spectra and number of segregating sites
                    num_snps[i] = 0;
                    avgn = 0;
                    for (j = 0; j < segsites; j++)
                        {
                            pop_type = types[j] & pop_mask[i];

                            // check if outgroup is aligned and different from reference
                            // else the reference base is assumed to be ancestral
                            if ((flag & BAM_OUTGROUP) && (ncov[outpop][j] > 0) && CHECK_BIT(types[j], outidx))
                                {
                                    freq = ncov[i][j] - bitcount64(pop_type);
                                }
                            else
                                {
                                    freq = bitcount64(pop_type);
                                }
                            if ((freq > 0) && (freq < ncov[i][j]))
                                {
                                    td[i] += dw[ncov[i][j]][freq];
                                    fwh[i] += hw[ncov[i][j]][freq];
                                    avgn += ncov[i][j];
                                    ++num_snps[i];
                                }
                        }

                    // finalize calculation of sfs statistics
                    n = (int)(((double)(avgn) / num_snps[i]) + 0.4999);
                    s = num_snps[i];
                    if ((n > 1) && (s > 0))
                        {
                            td[i] /= sqrt(e1[n] * s + e2[n] * s * (s - 1));
                            fwh[i] /= sqrt(((n - 2) * (s / a1[n]) / (6.0 * (n - 1))) + ((s * (s - 1) / (SQ(a1[n]) + a2[n])) *
                                           (18.0 * SQ(n) * (3.0 * n + 2.0) * a2[n+1] - (88.0 * n * n * n + 9.0 * SQ(n) - 13.0 * n + 6.0)) / (9.0 * n * (SQ(n-1)))));
                        }
                    else
                        {
                            td[i] = std::numeric_limits<double>::quiet_NaN();
                            fwh[i] = std::numeric_limits<double>::quiet_NaN();
                        }
                }
            else
                {
                    td[i] = std::numeric_limits<double>::quiet_NaN();
                    fwh[i] = std::numeric_limits<double>::quiet_NaN();
                }
        }
    return 0;
}

int
sfsData::printSFS(const std::string scaffold)
{
    int i = 0;
    std::stringstream out;

    out << scaffold << '\t' << beg + 1 << '\t' << end + 1;
    for (i = 0; i < npops; i++)
        {
            out << "\tns[" << sm->popul[i] << "]:\t" << ns[i];
            if (isnan(td[i]))
                {
                    out << "\tD[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                }
            else
                {
                    out << "\tD[" << sm->popul[i] << "]:";
                    out << '\t' << std::fixed << std::setprecision(5) << td[i];
                }
            if (isnan(fwh[i]))
                {
                    out << "\tH[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
                }
            else
                {
                    out << "\tH[" << sm->popul[i] << "]:";
                    out << '\t' << std::fixed << std::setprecision(5) << fwh[i];
                }
        }
    std::cout << out.str() << std::endl;
    return 0;
}

sfsData::sfsData(const popbamOptions &p)
{
    // inherit values from popbamOptions
    bamfile = p.bamfile;
    flag = p.flag;
    minDepth = p.minDepth;
    maxDepth = p.maxDepth;
    minRMSQ = p.minRMSQ;
    minSNPQ = p.minSNPQ;
    minMapQ = p.minMapQ;
    minBaseQ = p.minBaseQ;
    hetPrior = p.hetPrior;
    minSites = p.minSites;
    minPop = p.minPop;

    // initialize native variables
    derived_type = SFS;
    outidx = 0;
}

int
sfsData::allocSFS(void)
{
    int length = end - beg;

    segsites = 0;
    try
        {
            types = new uint64_t [length]();
            ns = new uint32_t [npops]();
            ncov = new uint32_t* [npops];
            pop_mask = new uint64_t [npops]();
            pop_nsmpl = new uint8_t [npops]();
            pop_cov = new uint32_t [length]();
            num_snps = new int [npops]();
            td = new double [npops]();
            fwh = new double [npops]();
            for (int i = 0; i < npops; ++i)
                {
                    ncov[i] = new uint32_t [length]();
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    return 0;
}

sfsData::~sfsData(void)
{
    delete [] pop_mask;
    delete [] ns;
    delete [] types;
    delete [] pop_nsmpl;
    delete [] pop_cov;
    delete [] num_snps;
    delete [] td;
    delete [] fwh;
    for (int i = 0; i < npops; ++i)
        {
            delete [] ncov[i];
        }
    delete [] ncov;
}

int
sfsData::calc_dw(void)
{
    int i = 0;

    dw = new double* [sm->n+1];
    for (i = 0; i <= sm->n; ++i)
        {
            dw[i] = new double [sm->n+1]();
        }
    for (int n = 2; n <= sm->n; ++n)
        {
            for(i = 1; i <= sm->n; ++i)
                {
                    dw[n][i] = (((2.0 * i * (n - i)) / (SQ(n-1))) - (1.0 / a1[n]));
                }
        }
    return 0;
}

int
sfsData::assignOutpop(void)
{
    uint64_t u = 0x1ULL << outidx;

    for (int i = 0; i < npops; ++i)
        {
            if (pop_mask[i] & u)
                {
                    outpop = i;
                }
        }
    return 0;
}

int
sfsData::calc_hw(void)
{
    int i = 0;

    hw = new double* [sm->n+1];
    for (i = 0; i <= sm->n; ++i)
        {
            hw[i] = new double [sm->n+1]();
        }
    for (int n = 2; n <= sm->n; ++n)
        {
            for (i = 1; i <= sm->n; ++i)
                {
                    hw[n][i] = ((1.0 / a1[n]) - ((double)(i) / (n - 1)));
                }
        }
    return 0;
}

int
sfsData::calc_a1(void)
{
    a1 = new double [sm->n+1];
    a1[0] = a1[1] = 1.0;

    // consider all sample sizes
    for (int i = 2; i <= sm->n; i++)
        {
            a1[i] = 0;
            for (int j = 1; j < i; j++)
                {
                    a1[i] += 1.0 / (double)(j);
                }
        }
    return 0;
}

int
sfsData::calc_a2(void)
{
    a2 = new double [sm->n+2];
    a2[0] = a2[1] = 1.0;

    // consider all sample sizes
    for (int i = 2; i <= sm->n+1; i++)
        {
            a2[i] = 0;
            for (int j = 1; j < i; j++)
                {
                    a2[i] += 1.0 / (double)SQ(j);
                }
        }
    return 0;
}

int
sfsData::calc_e1(void)
{
    e1 = new double [sm->n+1];
    e1[0] = e1[1] = 1.0;

    for (int i = 2; i <= sm->n; i++)
        {
            double b1 = (i + 1.0) / (3.0 * (i-1));
            e1[i] = (b1 - (1.0 / a1[i])) / a1[i];
        }
    return 0;
}

int
sfsData::calc_e2(void)
{
    e2 = new double [sm->n+1];
    e2[0] = e2[1] = 1.0;

    for (int i = 2; i <= sm->n; i++)
        {
            double b2 = (2.0 * (SQ(i) + i + 3.0)) / (9.0 * i * (i-1));
            e2[i] = (b2 - ((i + 2.0) / (a1[i] * i)) + (a2[i] / SQ(a1[i]))) / (SQ(a1[i]) + a2[i]);
        }
    return 0;
}

void
usageSFS(const std::string msg)
{
    std::cerr << msg << std::endl << std::endl;
    std::cerr << "Usage:   popbam sfs [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options: -i          base qualities are Illumina 1.3+               [ default: Sanger ]" << std::endl;
    std::cerr << "         -h  FILE    Input header file                              [ default: none ]" << std::endl;
    std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
    std::cerr << "         -p  STR     sample name of outgroup                        [ default: reference ]" << std::endl;
    std::cerr << "         -k  FLT     minimum proportion of sites covered in window  [ default: 0.5 ]" << std::endl;
    std::cerr << "         -n  FLT     minimum proportion of population covered       [ default: 1.0 ]" << std::endl;
    std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
    std::cerr << "         -m  INT     minimum read coverage                          [ default: 3 ]" << std::endl;
    std::cerr << "         -x  INT     maximum read coverage                          [ default: 255 ]" << std::endl;
    std::cerr << "         -q  INT     minimum rms mapping quality                    [ default: 25 ]" << std::endl;
    std::cerr << "         -s  INT     minimum snp quality                            [ default: 25 ]" << std::endl;
    std::cerr << "         -a  INT     minimum map quality                            [ default: 13 ]" << std::endl;
    std::cerr << "         -b  INT     minimum base quality                           [ default: 13 ]" << std::endl;
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
}

