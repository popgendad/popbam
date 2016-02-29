/** \file pop_utils.cpp
 *  \brief Utility functions evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.5
 *
 * Notes:
 * Byte/bit ordering of unsigned long long consensus base data type ("cb"):
 * Byte 1:     Boolean flags (cb[i]&0xff)
 *    bit 1:      Does this site pass the quality filters? (cb[i]&0x1)
 *    bit 2:      Is there a variant present at this site? (cb[i]&0x2)
 *    bit 3:      Not implemented
 *    bit 4:      Not implemented
 *    bit 5:      Not implemented
 *    bit 6:      Not implemented
 *    bit 7:      Not implemented
 *    bit 8:      Not implemented
 * Byte 2:     unsigned char-- the IUPAC consensus genotype (cb[i]>>8)&0xff
 *    bits 1-2:   unsigned char-- the IUPAC base call for allele 1 (cb[i]>>8)&0x3
 *    bits 3-4:   unsigned char-- the IUPAC base call for allele 2 (cb[i]>>10)&0x3
 *    bits 5-8:   Not implemented
 * Bytes 3-4:  unsigned short-- the number of reads mapped to that site (cb[i]>>16)&0xffff
 * Bytes 5-6:  unsigned short-- the SNP quality score (snpQ) (cb[i] >> 32)&0xffff
 * Bytes 7-8:  unsigned short-- the root-mean quality score (rmsQ) (cb[i]>>48)&0xffff
 *
**/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdint>
#include <cstddef>
#include <cstring>

#include "pop_sample.hpp"
#include "pop_utils.hpp"
#include "tables.h"

#define lfact(n) lgamma(n+1)

uint64_t
gl2cns(float q[16], uint16_t k)
{
    uint8_t i = 0;
    uint8_t j = 0;
    uint16_t min_ij = 0;
    uint64_t snp_quality = 0;
    uint64_t num_reads = 0;
    uint64_t genotype = 0;
    float min = FLT_MAX;
    float min_next = FLT_MAX;
    float likelihood = 0.0;

    for (i = 0; i < NBASES; ++i)
        {
            for (j = i; j < NBASES; ++j)
                {
                    likelihood = q[i << 2 | j];
                    if (likelihood < min)
                        {
                            min_ij = i << 2 | j;
                            min_next = min;
                            min = likelihood;
                        }
                    else if (likelihood < min_next)
                        {
                            min_next = likelihood;
                        }
                }
        }

    // return consensus base
    snp_quality = (uint64_t)((min_next - min) + 0.499) << (CHAR_BIT * 4);
    num_reads = (uint64_t)(k) << (CHAR_BIT * 2);
    genotype = (uint64_t)(min_ij) << CHAR_BIT;
    return snp_quality + num_reads + genotype;
}

uint64_t
qual_filter (int num_samples, uint64_t *cb, int min_rmsQ, int min_depth,
             int max_depth)
{
    int i = 0;
    uint16_t rms = 0;
    uint16_t num_reads = 0;
    uint64_t coverage = 0;

    for (i = 0; i < num_samples; ++i)
        {
            rms = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
            num_reads = (cb[i] >> (CHAR_BIT * 2)) & 0xffff;

            if ((rms >= min_rmsQ) && (num_reads >= min_depth) && (num_reads <= max_depth))
                {
                    cb[i] |= 0x1ULL;
                    coverage |= 0x1ULL << i;
                }
        }
    return coverage;
}

int
seg_base (int num_samples, uint64_t *cb, char ref, int min_snpq)
{
    int i = 0;
    int j = 0;
    int k = 0;
    uint8_t genotype = 0;
    uint8_t allele1 = 0;
    uint8_t allele2 = 0;
    uint16_t snp_quality = 0;
    int baseCount[NBASES] = {0, 0, 0, 0};

    for (i = 0; i < num_samples; ++i)
        {
            genotype = (cb[i] >> CHAR_BIT) & 0xff;
            allele1 = (genotype >> 2) & 0x3;
            allele2 = genotype & 0x3;
            snp_quality = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;

            // if homozygous and different from reference with high SNP quality
            if ((allele1 == allele2) && (iupac[genotype] != ref))
                {
                    if (snp_quality >= min_snpq)
                        {
                            cb[i] |= 0x2ULL;
                            ++baseCount[allele1];
                        }
                    // if SNP quality is low, revert both alleles to the reference allele
                    else
                        {
                            cb[i] -= (genotype-iupac_rev[(int)ref]) << CHAR_BIT;
                            cb[i] -= (genotype-iupac_rev[(int)ref]) << (CHAR_BIT+2);
                        }
                }
            else if ((allele1 != allele2) && (snp_quality >= min_snpq))
                {
                    cb[i] |= 0x2ULL;
                    if (iupac[allele1] != ref)
                        {
                            ++baseCount[allele1];
                        }
                    else
                        {
                            ++baseCount[allele2];
                        }
                }
        }

    // check for infinite sites model
    for (i = 0, j = 0, k = 0; i < NBASES; ++i)
        {
            if (baseCount[i] > 0)
                {
                    ++j;
                    k = i;
                }
        }

    //if (j > 1)
    //{
    //    return -1;
    //}
    //else
    //{
    return baseCount[k];
    //}
}

void
clean_hets (int num_samples, uint64_t *cb, int ref, int min_snpq)
{
    int i = 0;
    uint16_t snp_quality = 0;
    uint8_t genotype = 0;
    uint8_t allele1 = 0;
    uint8_t allele2 = 0;

    for (i = 0; i < num_samples; ++i)
        {
            genotype = (cb[i] >> CHAR_BIT) & 0xff;
            allele1 = (genotype >> 2) & 0x3;
            allele2 = genotype & 0x3;
            snp_quality = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;

            // if heterozygous and high quality SNP--make homozygous derived
            if ((allele1 != allele2) && (snp_quality >= min_snpq))
                {
                    if (allele1 == iupac_rev[ref])
                        {
                            cb[i] += (allele2 - allele1) << (CHAR_BIT+2);
                        }
                    if (allele2 == iupac_rev[ref])
                        {
                            cb[i] -= (allele2 - allele1) << CHAR_BIT;
                        }
                }
            // if heterozygous but poor quality--make homozygous ancestral
            if ((allele1 != allele2) && (snp_quality < min_snpq))
                {
                    if (allele1 != iupac_rev[ref])
                        {
                            cb[i] += (allele2 - allele1) << (CHAR_BIT + 2);
                        }
                    if (allele2 != iupac_rev[ref])
                        {
                            cb[i] -= (allele2 - allele1) << CHAR_BIT;
                        }
                }
        }
}

double *
logbinomial_table (const int n_size)
{
    int k = 0;
    int n = 0;
    double *logbinom = nullptr;

    logbinom = (double*)calloc (n_size * n_size, sizeof(double));
    for (n = 1; n < n_size; ++n)
        {
            double lfn = lfact(n);
            for (k = 1; k <= n; ++k)
                {
                    logbinom[n << 8 | k] = lfn - lfact(k) - lfact(n-k);
                }
        }
    return logbinom;
}

char *
get_refid (char *htext)
{
    int z = 0;
    const int idblock = 200;
    char *u = NULL;
    char *v = NULL;
    char *w = NULL;
    char *refid = NULL;

    u = htext;
    v = strstr (htext, "AS:");
    if (!v)
        {
            fputs("Unable to parse reference sequence name\nBe sure "
                  "the AS tag is defined in the sequence dictionary", stderr);
            exit (EXIT_FAILURE);
        }
    u = v + 3;
    for (z = 0, w = (char*)u; *w && (*w != '\t') && (*w != '\n'); ++w, ++z);
    refid = (char*) malloc (idblock * sizeof(char));
    refid[0] = '\0';
    strncpy (refid, u, z);
    refid[z] = '\0';
    return refid;
}

int
fetch_func (const bam1_t *b, void *data)
{
    bam_plbuf_t *buf;

    buf = (bam_plbuf_t*)data;
    bam_plbuf_push (b, buf);
    return 0;
}
