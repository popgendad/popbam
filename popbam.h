/** \file popbam.h
 *  \brief Header for the popbam program
 *  \author Daniel Garrigan
 *  \version 0.4
 */
#ifndef POPBAM_H
#define POPBAM_H

// If we're in C++ land, define these values so we can use various
// macros from the C99 spec (i.e. macros in stdint.h and inttypes.h
// that aren't included by default when compiling with C++ compiler).
#ifdef __cplusplus
#define __STDC_LIMIT_MACROS
#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS
#endif

// GCC specific stuff here
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

///
/// Include headers
///
#include <iostream>
#include <fstream>
#include <new>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <climits>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cerrno>
#include "faidx.h"
#include "sam.h"
#include "kstring.h"

#ifdef _MSC_VER
#define isnan(x) _isnan(x)
#endif

///
/// Definitions
///

/*! \def BAM_VARIANT
 *  \brief Flag for the -v command line switch-- only output variable sites
 */
#define BAM_VARIANT 0x01

/*! \def BAM_ILLUMINA
 *  \brief Flag for the -i command line switch-- qualities are in Illumina 1.3+ format
 */
#define BAM_ILLUMINA 0x02

/*! \def BAM_WINDOW
 *  \brief Flag for the -w command line switch-- do a sliding window analysis
 */
#define BAM_WINDOW 0x04

/*! \def BAM_MINPOPSAMPLE
 *  \brief Flag for user designated minimum sample sizes per population
 */
// TODO: This is deprecated
#define BAM_MINPOPSAMPLE 0x08

/*! \def BAM_SUBSTITUTE
 *  \brief Flag for only counting fixed substitutions in diverge function
 */
#define BAM_SUBSTITUTE 0x10

/*! \def BAM_HETEROZYGOTE
 *  \brief Flag for outputting heterozygous site calls in snp function
 */
#define BAM_HETEROZYGOTE 0x20

/*! \def BAM_OUTGROUP
 *  \brief Flag for changing the outgroup from the reference
 */
#define BAM_OUTGROUP 0x40

/*! \def BAM_HEADERIN
 *  \brief Flag for presence of user BAM header file
 */
#define BAM_HEADERIN 0x80

/*! \def BAM_NOSINGLETONS
 *  \brief Flag to exclude singleton polymorphisms from the analysis
 */
#define BAM_NOSINGLETONS 0x100

/*! \def POPBAM_RELEASE
 *  \brief Version number of popbam program
 */
#define POPBAM_RELEASE "0.4b"

/*! \def NBASES
 *  \brief The number of possible bases
 */
#define NBASES 4

/*! \def IUPAC_N
 *  \brief The integer representation of the IUPAC ambiguity symbol 'N'
 */
#define IUPAC_N 0xf

/*! \def KB
 *  \brief Integer for length of a kilobase
 */
#define KB 1000

/*! \def CHECK_BIT(var,pos)
 *  \brief A macro to check if a bit is set at pos in the unsigned long long var
 */
#define CHECK_BIT(var,pos) ((var) & (0x1ULL << (pos)))

/*! \def SEG_IDX(segsite)
 *  \brief A macro access index of a segregating site
 */
#define SEG_IDX(seg) (((seg) - 1) / 64)

/*! \def UTIDX(nrows,row,col)
 *  \brief A macro to access index of 1D array from 2D data structure
 */
#define UTIDX(n,i,j) ((2*(n)-((i)+1))*(((i)+1)-1)/2-((i)+1)+((j)+1)-1)

/*! \def SQ(x)
 *  \brief A macro to calculate the square of its argument
 */
#define SQ(x) ((x) * (x))

/*! \def BINOM(x)
 *  \brief A macro to calculate binomial coefficient
 */
#define BINOM(x) ((x) * ((x) - 1) / 2)

//
// Define data structures
//

/*!
 * struct hData_t
 * \brief A structure to represent a haplotype data set
 */
typedef struct
{
	unsigned long long **seq;         //!< binary encoding of haplotype data
	unsigned int *pos;                //!< reference coordinate for each position
	unsigned int *idx;                //!< position index of each segregating site
	unsigned char *ref;               //!< reference allele at each position
	unsigned char **base;             //!< consensus base at each position in each individual
	unsigned short **rms;             //!< root mean square mapping score at each position
	unsigned short **snpq;            //!< SNP quality score at each position
	unsigned short **num_reads;       //!< number of reads at each position in each individual
} hData_t;

/*!
 * \struct bam_sample_t
 * \brief A structure to represent a sample in the BAM file
 */
typedef struct __bam_sample_t
{
	int npops;                        //!< Number of populations in the BAM file
	int b;                            //!< Counter for population configuration
	int n;                            //!< Number of samples in the BAM file
	int m;                            //!< Counter for sample configuration
	char **smpl;                      //!< Pointer to array of sample names
	char **popul;                     //!< Pointer to array of population names
	void *rg2smid;                    //!< Pointer to hash for read group to sample id lookup
	void *sm2popid;                   //!< Pointer to hash for sample to population id lookup
	void *sm2id;                      //!< Pointer to hash for sample to identifier lookup
	void *pop2sm;                     //!< Pointer to hash for population to sample lookup
} bam_sample_t;

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

//
// Define some global variables
//

/*! \def popbam_func_t
 *  \brief A enum data type that holds the popbam function identifier
 */
enum popbam_func_t {SNP, FASTA, DIVERGE, HAPLO, TREE, NUCDIV, LD, SFS};

///
/// Define classes
///

class popbamOptions
{
public:
	// constructor
	popbamOptions(int, char**);

	// destructor
	~popbamOptions() {}

	// member variables
	samfile_t *bam_in;                      //!< BAM input file stream
	faidx_t *fai_file;                      //!< Fasta reference file index
	bam_index_t *idx;                       //!< Pointer to the BAM input file index
	bam_header_t *h;                        //!< Pointer to the header text for the input BAM file
	unsigned short flag;                    //!< Bit flag to hold user options
	int output;                             //!< Analysis output option
	int errorCount;                         //!< Flag to indicate error in reading user options
	int minDepth;                           //!< User-specified minimumm read depth
	int maxDepth;                           //!< User-specified maximum read depth
	int minRMSQ;                            //!< User-specified minimum rms mapping quality
	int minSNPQ;                            //!< User-specified minimum SNP quality score
	unsigned int winSize;                   //!< User-specified window size in kilobases
	unsigned char minMapQ;                  //!< User-specified minimum individual read mapping quality
	unsigned char minBaseQ;                 //!< User-specified minimum inidividual base quality
	double minSites;                        //!< User-specified minimum number of aligned sites to perform analysis
	double minPop;                          //!< Minimum proportion of samples present
	double hetPrior;                        //!< Prior probability for calling heterozygous genotypes
	std::string dist;                       //!< Pointer to the name of the desired distance metric	(-d switch)
	std::string bamfile;                    //!< File name for the input BAM file
	std::string reffile;                    //!< File name for the input reference Fasta file
	std::string headfile;                   //!< File name for optional BAM header input file
	std::string region;                     //!< Region on which to perform the analysis
	std::string errorMsg;                   //!< String to hold any error messages
	std::string popFunc;                    //!< The popbam function being invoked

	// member functions
	int checkBAM(void);
};

/*!
 * \class popbamData
 * \brief The abstract base class for passing parameters and data
 */
class popbamData
{
	public:
		// default constructor
		popbamData();

		// destructor
		~popbamData() {}

		// member functions
		int assignPops(const popbamOptions *p);

		// virtual functions
		virtual void printUsage(std::string) = 0;

		// member variables
		std::string bamfile;                    //!< Name of bamfile used for indexing purposes
		bam_sample_t *sm;                       //!< Pointer to the sample information for the input BAM file
		char *ref_base;                         //!< Reference sequence string for specified region
		int tid;                                //!< Reference chromosome/scaffold identifier
		int beg;                                //!< Reference coordinate of the beginning of the current region
		int end;                                //!< Reference coordinate of the end of current region
		int len;                                //!< Length of the reference sequence for current region
		unsigned short flag;                    //!< Bit flag to hold user options
		int num_sites;                          //!< Total number of aligned sites
		int segsites;                           //!< Total number of segregating sites in entire sample
		unsigned char *pop_nsmpl;               //!< Sample size per population
		unsigned long long *types;              //!< The site type for each aligned site
		unsigned long long *pop_mask;           //!< Bit mask for which individuals are in which population
		int minDepth;                           //!< User-specified minimumm read depth
		int maxDepth;                           //!< User-specified maximum read depth
		int minRMSQ;                            //!< User-specified minimum rms mapping quality
		int minSNPQ;                            //!< User-specified minimum SNP quality score
		unsigned char minMapQ;                  //!< User-specified minimum individual read mapping quality
		unsigned char minBaseQ;                 //!< User-specified minimum inidividual base quality
		double hetPrior;                        //!< Prior probability of heterozygous genotype
		errmod_t *em;                           //!< Error model data structure
		popbam_func_t derived_type;             //!< Type of the derived class
};

///
/// Function prototypes
///

// Entry points to the main popbam functions
extern int mainSNP(int, char**);
extern int mainFasta(int, char**);
extern int mainHaplo(int, char**);
extern int mainDiverge(int, char**);
extern int mainTree(int, char**);
extern int mainNucdiv(int, char**);
extern int mainLD(int, char**);
extern int mainSFS(int, char**);

/*!
 * \fn inline unsigned int log2int(const unsigned int val)
 * \brief Returns integer of log-base2 of val
 * \param val The input value
 */
#ifdef _MSC_VER
static inline unsigned int log2int(const unsigned int va)
{
	unsigned int re;

  __asm {
	bsr eax, va
	mov re, eax
  }

  return re;
}
#else
extern __inline unsigned int log2int(const unsigned int va)
{
	unsigned int re;

	asm ( "\tbsr  %1, %0\n"
		 : "=r" (re)
		 : "r"  (va)
		 );

	return re;
}
#endif

/*!
 * \fn inline unsigned int bitcount64(unsigned long long x)
 * \brief Function count the number of bits set in a 64-bit integer
 * \param x the 64-bit integer
 * \return unsigned integer
 * Returns the number of bits set
 */
inline unsigned short bitcount64(unsigned long long x)
{
	x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
	return (x * 0x0101010101010101ULL) >> 56;
}

/*!
 * \fn inline unsigned int hamming_distance(unsigned long long x, unsigned long long y)
 * \brief Function to compute the hamming distance between two 64-bit integers
 * \param x the first 64-bit integer
 * \param y the second 64-bit integer
 * \return unsigned int of hamming distance
 * Returns the number of bits set
 */
inline unsigned int hamming_distance(unsigned long long x, unsigned long long y)
{
	unsigned int dist = 0;
	unsigned long long val = x^y;

	// count the number of set bits
	while (val)
	{
		++dist;
		val &= val-1;
	}

	return dist;
}

inline unsigned long long calculateSiteType(int n, unsigned long long *cb)
{
	unsigned long long site_type = 0;
	for (int i=0; i < n; i++)
		if ((cb[i] & 0x3ULL) == 0x3ULL)
			site_type |= 0x1ULL << i;
	return site_type;
}

/*!
 * \fn int popbam_usage(void)
 * \brief Prints general command usage options to stdout
 */
extern int popbam_usage(void);

/*!
 * \fn bool is_file_exist(const char *fileName)
 * \brief Checks whether a file exists on disk
 */
extern bool is_file_exist(const char *fileName);

/*!
 * \fn bam_sample_t *bam_smpl_init(void)
 * \brief Initialize the sample data structure
 */
extern bam_sample_t *bam_smpl_init(void);

/*!
 * \fn int bam_smpl_add(bam_sample_t *sm, const popbamOptions *p)
 * \brief Add a sample data structure
 * \param sm Pointer to sample data structure
 * \param abs Pointer to name of input BAM file
 * \param txt Pointer to unformatted BAM header txt
 */
extern int bam_smpl_add(bam_sample_t *sm, const popbamOptions *op);

/*!
 * \fn int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
 * \brief Get the sample id of a read group
 * \param sm Pointer to sample data structure
 * \param fn Pointer to the name of the input BAM file
 * \param rg Pointer to the name of the read group
 * \param str Pointer to the name of the sample
 */
extern int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str);

/*!
 * \fn int bam_smpl_sm2popid(const bam_sample_t *sm, const char *fn, const char *smpl, kstring_t *str)
 * \brief Get the population id of a sample
 * \param sm Pointer to sample data structure
 * \param fn Pointer to the name of the input BAM file
 * \param smpl Pointer to the name of the sample
 * \param str Pointer to the name of the population
 */
extern int bam_smpl_sm2popid(const bam_sample_t *sm, const char *fn, const char *smpl, kstring_t *str);
/*!
 * \fn void bam_smpl_destroy(bam_sample_t *sm)
 * \brief Free a sample data structure from memory
 * \param sm Pointer to sample data structure
 */
extern void bam_smpl_destroy(bam_sample_t *sm);

/*!
 * \fn char *get_refid(char *htext)
 * \brief Function to extract reference identifier from BAM header
 * \param htext Pointer to unformatted BAM header text
 */
extern char *get_refid(char *htext);

/*!
 * \fn int fetch_func (const bam1_t *b, void *data)
 * \brief Assigns functions to the pileup push
 * \param b Pointer to the alignment structure
 * \param data User defined data structure
 */
extern int fetch_func(const bam1_t *b, void *data);

/*!
 * \fn unsigned long long qualFilter(int num_samples, unsigned long long *cb, int min_rmsQ, int min_depth, int max_depth)
 * \brief Filters data based on quality threshholds
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param min_rmsQ  Minimum root-mean square of mapping quality for site to be considered
 * \param min_depth  Minimum read depth per individual for site to be considered
 * \param max_depth  Maximum read depth per individual for site to be considered
 */
extern unsigned long long qualFilter(int num_samples, unsigned long long *cb, int min_rmsQ, int min_depth, int max_depth);

/*!
 * \fn int segbase(int num_samples, unsigned long long *cb, char ref, int min_snpq)
 * \brief Determines whether a base position is segregating or not
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */
extern int segBase(int num_samples, unsigned long long *cb, char ref, int min_snpq);

/*!
 * \fn void cleanHeterozygotes(int num_samples, unsigned long long *cb, int ref, int min_snpq)
 * \brief Reconfigures heterozygous base calls
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */
extern void cleanHeterozygotes(int num_samples, unsigned long long *cb, int ref, int min_snpq);

/*!
 * \fn unsigned long long gl2cns(float q[16], unsigned short k)
 * \brief Calculates a consensus base call from genotype likelihoods
 * \param q  Probabilites associated with each base
 * \param k  Number of reads mapping to a position in an individual
 */
extern unsigned long long gl2cns(float q[16], unsigned short k);

/*!
 * \fn errmod_t *errmod_init(float depcorr)
 * \brief Initialize the error model data structure
 * \param depcorr The constant for the dependency correlation
 */
extern errmod_t *errmod_init(float depcorr);

/*!
 * \fn void errmod_destroy(errmod_t *em)
 * \brief Deallocate memory for error model data structure
 * \param em Pointer to error model data structure
 */
extern void errmod_destroy(errmod_t *em);

/*!
 * \fn int errmod_cal(const errmod_t *em, unsigned short n, int m, unsigned short *bases, float *q)
 * \brief Calculates probability for error model
 * \param em Pointer to the error model data structure
 * \param n The number of bases
 * \param m The maximum base
 * \param bases[i] qual:6, strand:1, base:4
 * \param q[i*m+j] Phred-scaled likelihood of (i,j)
 */
extern int errmod_cal(const errmod_t *em, unsigned short n, int m, unsigned short *bases, float *q);

/*!
 * \fn void fatalError(const char *msg, char* file, int line, void(*err_func)(void))
 * \brief Prints error message and exits program
 * \param msg Pointer to string containing error message
 */
extern void fatalError(const std::string msg);

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
extern int bam_parse_region(bam_header_t *header, std::string region, int *ref_id, int *begin, int *end);

/*!
 * \fn const char* svn_version(void)
 * \brief Retrieves a pointer to the subversion revision number
 */
const char* svn_version(void);

#endif
