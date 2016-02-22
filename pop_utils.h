/** \file pop_utils.h
 *  \brief Header for the pop_utils.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
 */
#ifndef POP_UTILS_H
#define POP_UTILS_H

/*!
 * \struct errmod_coef_t
 * \brief A structure to hold the coefficients necessary in the error model
 */
typedef struct __errmod_coef_t
{
    double *fk;                       //!< Pointer to
    double *beta;                     //!< Pointer to
    double *lhet;                     //!< Pointer to
} errmod_coef_t;

/*!
 * \struct errmod_t
 * \brief A structure to hold data for the error model
 */
typedef struct __errmod_t
{
    double depcorr;                   //!< Dependency correlation
    errmod_coef_t *coef;              //!< Pre-computed coefficients
} errmod_t;

/*!
 * \struct call_aux_t
 * \brief A structure to hold auxiliary information for use in error model
 */
typedef struct __call_aux_t
{
    double fsum[16];                  //!< Array of
    double bsum[16];                  //!< Array of
    unsigned int c[16];               //!< Array of
} call_aux_t;

///
/// Function prototypes
///

/*!
 * \fn uint64_t gl2cns(float q[16], uint16_t k)
 * \brief Calculates a consensus base call from genotype likelihoods
 * \param q  Probabilites associated with each base
 * \param k  Number of reads mapping to a position in an individual
 */

extern uint64_t gl2cns(float q[16], uint16_t k);

/*!
 * \fn uint64_t  qualFilter(int num_samples, uint64_t *cb, int min_rmsQ, int min_depth, int max_depth)
 * \brief Filters data based on quality threshholds
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param min_rmsQ  Minimum root-mean square of mapping quality for site to be considered
 * \param min_depth  Minimum read depth per individual for site to be considered
 * \param max_depth  Maximum read depth per individual for site to be considered
 */

extern uint64_t qualFilter(int num_samples, uint64_t *cb, int min_rmsQ, int min_depth, int max_depth);

/*!
 * \fn int segbase(int num_samples, uint64_t *cb, char ref, int min_snpq)
 * \brief Determines whether a base position is segregating or not
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */

extern int segBase(int num_samples, uint64_t *cb, char ref, int min_snpq);

/*!
 * \fn void cleanHeterozygotes(int num_samples, uint64_t *cb, int ref, int min_snpq)
 * \brief Reconfigures heterozygous base calls
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */

extern void cleanHeterozygotes(int num_samples, uint64_t *cb, int ref, int min_snpq);

/*!
 * \fn static double* logbinomial_table(const int n_size)
 * \brief Construct logbinomial table
 * \param n_size Size of the table
 */

extern double* logbinomial_table(const int n_size);

/*!
 * \fn static errmod_coef_t* cal_coef(double depcorr, double eta)
 * \brief Calculate coefficients in error model
 * \param depcorr The constant for the dependency correlation
 * \param eta
 */

extern errmod_coef_t* cal_coef(double depcorr, double eta);

/*!
 * \fn errmod_t *errmod_init(float depcorr)
 * \brief Initialize the error model data structure
 * \param depcorr The constant for the dependency correlation
 */

extern errmod_t* errmod_init(float depcorr);

/*!
 * \fn void errmod_destroy(errmod_t *em)
 * \brief Deallocate memory for error model data structure
 * \param em Pointer to error model data structure
 */

extern void errmod_destroy(errmod_t *em);

/*!
 * \fn int errmod_cal(const errmod_t *em, uint16_t n, int m, uint16_t *bases, float *q)
 * \brief Calculates probability for error model
 * \param em Pointer to the error model data structure
 * \param n The number of bases
 * \param m The maximum base
 * \param bases[i] qual:6, strand:1, base:4
 * \param q[i*m+j] Phred-scaled likelihood of (i,j)
 */

extern int errmod_cal(const errmod_t *em, uint16_t n, int m, uint16_t *bases, float *q);

/*!
 * \fn void bam_init_header_hash(bam_header_t *header)
 * \brief Initialize the BAM header
 * \param header Pointer to the BAM header text
 */

extern void bam_init_header_hash(bam_header_t *header);

/*!
 * \fn int bam_parse_region(bam_header_t *header, std::string region, int *ref_id, int *begin, int *end)
  \brief Parse a region in the format: "chr2:100,000-200,000".
  \param header Pointer to the header structure
  \param str String to be parsed
  \param ref_id The returned chromosome ID
  \param begin The returned start coordinate
  \param end The returned end coordinate
  \return 0 on success; -1 on failure
 */

extern int bam_parse_region(bam_header_t *header, std::string region, int *ref_id, int *beg, int *end);

/*!
 * \fn char *get_refid(char *htext)
 * \brief Function to extract reference identifier from BAM header
 * \param htext Pointer to unformatted BAM header text
 */

extern char *get_refid(char *htext);

/*!
 * \fn bool is_file_exist(const char *fileName)
 * \brief Checks whether a file exists on disk
 */

extern bool is_file_exist(const char *fileName);

/*!
 * \fn int fetch_func (const bam1_t *b, void *data)
 * \brief Assigns functions to the pileup push
 * \param b Pointer to the alignment structure
 * \param data User defined data structure
 */
extern int fetch_func(const bam1_t *b, void *data);

/*!
 * \fn void fatalError(const std::string msg)
 * \brief Prints error message and exits program
 * \param msg Pointer to string containing error message
 */

extern void fatalError(const std::string msg);


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
                            if (baseQ < 0)
                            {
                              exit (EXIT_FAILURE);
                            }
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

/*!
 * \fn inline uint32_t log2int(const uint32_t val)
 * \brief Returns integer of log-base2 of val
 * \param val The input value
 */

extern __inline uint32_t log2int(const uint32_t va)
{
    uint32_t re;

    asm ( "\tbsr  %1, %0\n"
          : "=r" (re)
          : "r"  (va)
        );
    return re;
}

/*!
 * \fn inline uint16_t bitcount64(uint64_t x)
 * \brief Function count the number of bits set in a 64-bit integer
 * \param x the 64-bit integer
 * \return unsigned integer
 * Returns the number of bits set
 */

inline uint16_t bitcount64(uint64_t x)
{
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    return (x * 0x0101010101010101ULL) >> 56;
}

/*!
 * \fn inline uint32_t hamming_distance(uint64_t x, uint64_t y)
 * \brief Function to compute the hamming distance between two 64-bit integers
 * \param x the first 64-bit integer
 * \param y the second 64-bit integer
 * \return unsigned int of hamming distance
 * Returns the number of bits set
 */

inline uint32_t hamming_distance(uint64_t x, uint64_t y)
{
    uint32_t dist = 0;
    uint64_t val = x^y;

    while (val)
        {
            ++dist;
            val &= val-1;
        }
    return dist;
}

/*!
 * \fn inline uint64_t calculateSiteType(int n, uint64_t *cb)
 * \brief Function to calculate site type representation
 * \param n
 * \param cb Consensus base data structure
 * \return encoded site type
 */

inline uint64_t calculateSiteType(int n, uint64_t *cb)
{
    uint64_t site_type = 0;

    for (int i = 0; i < n; i++)
        {
            if ((cb[i] & 0x3ULL) == 0x3ULL)
                {
                    site_type |= 0x1ULL << i;
                }
        }
    return site_type;
}

#endif
