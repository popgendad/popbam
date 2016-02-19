/** \file pop_utils.cpp
 *  \brief Utility functions evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.4
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
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <ifstream>
#include <string>
#include <new>
#include <vector>
#include "ksort.h"
#include "khash.h"
#include "kstring.h"
#include "bam.h"
#include "sam.h"
#include "pop_utils.h"
#include "pop_sample.h"
#include "popbam.h"
#include "tables.h"

#define lfact(n) lgamma(n+1)

typedef char *str_p;

KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)
KSORT_INIT_GENERIC(uint16_t)

template <class T> uint64_t* callBase(T *t, int n, const bam_pileup1_t *pl)
{
    int i = 0;
    int j = 0;
    int qq = 0;
    int baseQ = 0;
    int tmp_baseQ = 0;
    int b = 0;
    int si = -1;
    int rmsq = 0;
    int n_smpl = t->sm->n;
    uint16_t k = 0;
    uint16_t *bases = nullptr;
    uint64_t rms = 0;
    uint64_t *cb;
    std::vector<int> depth(n_smpl, 0);
    uint8_t *s = nullptr;
    float q[16];
    std::string msg;
    bam_pileup1_t ***p = nullptr;
    kstring_t buf;

    // allocate memory pileup data
    try
        {
            cb = new uint64_t [n_smpl]();
            p = new bam_pileup1_t** [n_smpl];
            for (i = 0; i < n_smpl; i++)
                {
                    p[i] = new bam_pileup1_t* [t->maxDepth];
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    memset(&buf, 0, sizeof(kstring_t));

    // partition pileup according to sample
    for (i = 0; i < n; i++)
        {
            if ((pl+i)->is_del || (pl+i)->is_refskip || ((pl+i)->b->core.flag & BAM_FUNMAP))
                {
                    continue;
                }
            s = bam_aux_get((pl+i)->b, "RG");

            // skip reads with no read group tag
            if (!s)
                {
                    continue;
                }
            else
                {
                    si = bam_smpl_rg2smid(t->sm, t->bamfile.c_str(), (char*)(s+1), &buf);
                }
            if (si < 0)
                {
                    si = bam_smpl_rg2smid(t->sm, t->bamfile.c_str(), 0, &buf);
                }
            if (si < 0)
                {
                    free(buf.s);
                    std::string rogue_rg(bam_aux2Z(s));
                    msg = "Problem assigning read group " + rogue_rg + " to a sample.\nPlease check BAM header for correct SM and PO tags";
                    fatalError(msg);
                }
            if (depth[si] < t->maxDepth)
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
                            tmp_baseQ = bam1_qual(p[j][i]->b)[p[j][i]->qpos];

                            if (t->flag & BAM_ILLUMINA)
                                {
                                    baseQ = tmp_baseQ > 31 ? tmp_baseQ - 31 : 0;
                                }
                            else
                                {
                                    baseQ = tmp_baseQ;
                                }
                            assert(baseQ >= 0);
                            if ((baseQ < t->minBaseQ) || (p[j][i]->b->core.qual < t->minMapQ))
                                {
                                    continue;
                                }
                            b = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p[j][i]->b), p[j][i]->qpos)];
                            if (b > 3)
                                {
                                    continue;
                                }
                            qq = baseQ < p[j][i]->b->core.qual ? baseQ : p[j][i]->b->core.qual;
                            if (qq < 4)
                                {
                                    qq = 4;
                                }
                            if (qq > 63)
                                {
                                    qq = 63;
                                }
                            bases[k++] = qq << 5 | (uint16_t)bam1_strand(p[j][i]->b) << 4 | b;
                            rmsq += SQ(p[j][i]->b->core.qual);
                        }

                    // calculate genotype likelihoods
                    errmod_cal(t->em, k, NBASES, bases, q);

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
    free(buf.s);

    return cb;
}

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
qualFilter(int num_samples, uint64_t *cb, int min_rmsQ, int min_depth, int max_depth)
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
segBase(int num_samples, uint64_t *cb, char ref, int min_snpq)
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
cleanHeterozygotes(int num_samples, uint64_t *cb, int ref, int min_snpq)
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
                            cb[i] += (allele2 - allele1) << (CHAR_BIT+2);
                        }
                    if (allele2 != iupac_rev[ref])
                        {
                            cb[i] -= (allele2 - allele1) << CHAR_BIT;
                        }
                }
        }
}

double*
logbinomial_table(const int n_size)
{
    int k = 0;
    int n = 0;
    double *logbinom = nullptr;

    logbinom = (double*)calloc(SQ(n_size), sizeof(double));
    for (n = 1; n < n_size; ++n)
        {
            double lfn = lfact(n);
            for (k = 1; k <= n; ++k)
                {
                    logbinom[n<<8|k] = lfn - lfact(k) - lfact(n-k);
                }
        }
    return logbinom;
}

errmod_coef_t*
cal_coef(double depcorr, double eta)
{
    int k = 0;
    int n = 0;
    int q = 0;
    long double sum = 0.0;
    long double sum1 = 0.0;
    double *lC = nullptr;
    errmod_coef_t *ec;

    ec = (errmod_coef_t*)calloc(1, sizeof(errmod_coef_t));

    // initialize ->fk
    ec->fk = (double*)calloc(256, sizeof(double));
    ec->fk[0] = 1.0;
    for (n = 1; n < 256; ++n)
        {
            ec->fk[n] = pow(1.0 - depcorr, n) * (1.0 - eta) + eta;
        }

    // initialize ->coef
    ec->beta = (double*)calloc(SQ(256) * 64, sizeof(double));

    lC = logbinomial_table(256);
    for (q = 1; q < 64; ++q)
        {
            double e = pow(10.0, -q/10.0);
            double le = log(e);
            double le1 = log(1.0 - e);
            for (n = 1; n <= 255; ++n)
                {
                    double *beta = ec->beta + (q << 16 | n << 8);
                    sum1 = sum = 0.0;
                    for (k = n; k >= 0; --k, sum1 = sum)
                        {
                            sum = sum1 + expl(lC[n << 8 | k] + k * le + (n - k) * le1);
                            beta[k] = -10.0 / M_LN10 * logl(sum1 / sum);
                        }
                }
        }

    // initialize ->lhet
    ec->lhet = (double*)calloc(SQ(256), sizeof(double));

    for (n = 0; n < 256; ++n)
        {
            for (k = 0; k < 256; ++k)
                {
                    ec->lhet[n << 8 | k] = lC[n << 8 | k] - M_LN2 * n;
                }
        }
    free(lC);
    return ec;
}

errmod_t*
errmod_init(float depcorr)
{
    errmod_t *em;

    em = (errmod_t*)calloc(1, sizeof(errmod_t));
    em->depcorr = depcorr;
    em->coef = cal_coef(depcorr, 0.03);
    return em;
}

void
errmod_destroy(errmod_t *em)
{
    if (em == 0)
        {
            return;
        }
    free(em->coef->lhet);
    free(em->coef->fk);
    free(em->coef->beta);
    free(em->coef);
    free(em);
}

// qual:6, strand:1, base:4
int
errmod_cal(const errmod_t *em, uint16_t n, int m, uint16_t *bases, float *q)
{
    call_aux_t aux;
    int i = 0;
    int j = 0;
    int k = 0;
    int w[32];

    memset(q, 0, SQ(m) * sizeof(float));
    if (n == 0)
        {
            return 0;
        }

    // calculate aux.esum and aux.fsum
    // then sample 255 bases
    if (n > 255)
        {
            ks_shuffle(uint16_t, n, bases);
            n = 255;
        }
    ks_introsort(uint16_t, n, bases);
    memset(w, 0, 32 * sizeof(int));
    memset(&aux, 0, sizeof(call_aux_t));

    // calculate esum and fsum
    for (j = n - 1; j >= 0; --j)
        {
            uint16_t b = bases[j];
            int qlty = b >> 5 < NBASES ? NBASES : b >> 5;
            if (qlty > 63)
                {
                    qlty = 63;
                }
            int basestrand = b & 0x1f;
            int base = b & 0xf;
            aux.fsum[base] += em->coef->fk[w[basestrand]];
            aux.bsum[base] += em->coef->fk[w[basestrand]] * em->coef->beta[qlty << 16 | n << 8 | aux.c[base]];
            ++aux.c[base];
            ++w[basestrand];
        }

    // generate likelihood
    for (j = 0; j < m; ++j)
        {
            float tmp1 = 0.0;
            float tmp3 = 0.0;
            int tmp2 = 0;

            // homozygous
            for (k = 0, tmp1 = tmp3 = 0.0, tmp2 = 0; k < m; ++k)
                {
                    if (k == j)
                        {
                            continue;
                        }
                    tmp1 += aux.bsum[k];
                    tmp2 += aux.c[k];
                    tmp3 += aux.fsum[k];
                }
            if (tmp2)
                {
                    q[j*m+j] = tmp1;
                }
            // heterozygous
            for (k = j + 1; k < m; ++k)
                {
                    int cjk = aux.c[j] + aux.c[k];
                    for (i = 0, tmp2 = 0, tmp1 = tmp3 = 0.0; i < m; ++i)
                        {
                            if ((i == j) || (i == k))
                                {
                                    continue;
                                }
                            tmp1 += aux.bsum[i];
                            tmp2 += aux.c[i];
                            tmp3 += aux.fsum[i];
                        }
                    if (tmp2)
                        {
                            q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk << 8 | aux.c[k]] + tmp1;
                        }
                    // all the bases are either j or k
                    else
                        {
                            q[j*m+k] = q[k*m+j] = -4.343 * em->coef->lhet[cjk << 8 | aux.c[k]];
                        }
                }

            for (k = 0; k < m; ++k)
                {
                    if (q[j*m+k] < 0.0)
                        {
                            q[j*m+k] = 0.0;
                        }
                }
        }
    return 0;
}

void
bam_init_header_hash(bam_header_t *header)
{
    int i = 0;

    if (header->hash == NULL)
        {
            int ret = 0;
            khiter_t iter;
            khash_t(s) *h;

            header->hash = h = kh_init(s);
            for (i = 0; i < header->n_targets; ++i)
                {
                    iter = kh_put(s, h, header->target_name[i], &ret);
                    kh_value(h, iter) = i;
                }
        }
}

int
bam_parse_region(bam_header_t *header, std::string region, int *ref_id, int *beg, int *end)
{
    std::size_t l = 0;
    std::size_t name_end = 0;
    khiter_t iter;
    khash_t(s) *h;

    bam_init_header_hash(header);
    h = (khash_t(s)*)header->hash;
    *ref_id = *beg = *end = -1;
    name_end = l = region.length();

    // remove spaces and commas
    std::string::iterator end_pos = std::remove(region.begin(), region.end(), ' ');
    region.erase(end_pos, region.end());
    end_pos = std::remove(region.begin(), region.end(), ',');
    region.erase(end_pos, region.end());
    l = region.length();

    // determine the sequence name
    // look for colon from the end
    name_end = region.find(":");
    if (name_end == std::string::npos)
        {
            name_end = l;
        }

    // check if this is really the end
    if (name_end < l)
        {
            std::string coords = region.substr(name_end + 1);
            std::size_t n_hyphen = std::count(coords.begin(), coords.end(), '-');
            std::size_t n_nondigits = coords.find_first_not_of("0123456789,-");

            // malformated region string; then take str as the name
            if ((n_nondigits != std::string::npos) || (n_hyphen > 1))
                {
                    name_end = l;
                }
            std::string scaffold_name = region.substr(0, name_end);

            // get a hash iterator to the scaffold name in the header
            iter = kh_get(s, h, scaffold_name.c_str());

            // cannot find the sequence name
            if (iter == kh_end(h))
                {
                    // try the entire region string as the lookup key
                    iter = kh_get(s, h, region.c_str());
                    if (iter == kh_end(h))
                        {
                            std::cerr << "Cannot find sequence name " << region << " in header" << std::endl;
                            return -1;
                        }
                }
        }
    else
        {
            iter = kh_get(s, h, region.c_str());
        }
    if (iter == kh_end(h))
        {
            return -1;
        }
    *ref_id = kh_val(h, iter);

    // parse the interval
    if (name_end < l)
        {
            std::string coords = region.substr(name_end+1);
            std::size_t parse = coords.find("-");
            std::string first = coords.substr(0, parse);
            *beg = atoi(first.c_str());
            if (*beg > 0)
                {
                    --*beg;
                }
            std::string last = coords.substr(parse+1);
            *end = atoi(last.c_str());
        }
    else
        {
            *beg = 0;
            *end = header->target_len[*ref_id];
        }
    return *beg <= *end ? 0 : -1;
}

char*
get_refid(char *htext)
{
    char *u = nullptr;
    char *v = nullptr;
    char *w = nullptr;
    const int idblock = 200;
    int z = 0;
    char *refid = nullptr;

    u = htext;
    v = strstr(htext, "AS:");
    if (!v)
        {
            fatalError("Unable to parse reference sequence name\nBe sure the AS tag is defined in the sequence dictionary");
        }
    u = v + 3;
    for (z = 0, w = (char*)u; *w && *w != '\t' && *w != '\n'; ++w, ++z);
    try
        {
            refid = new char [idblock];
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    refid[0] = '\0';
    strncpy(refid, u, z);
    refid[z] = '\0';
    return refid;
}

bool
is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int
fetch_func(const bam1_t *b, void *data)
{
    bam_plbuf_t *buf;

    buf = (bam_plbuf_t*)data;
    bam_plbuf_push(b, buf);
    return 0;
}

void
fatalError(const std::string msg)
{
    std::cerr << "popbam runtime error:" << std::endl;
    std::cerr << msg << std::endl;
    exit(EXIT_FAILURE);
}

