/** \file popbam.h
 *  \brief Header for the popbam program
 *  \author Daniel Garrigan
 *  \version 0.4
 */
#ifndef POPBAM_H
#define POPBAM_H

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
	samfile_t *bam_in;                //!< BAM input file stream
	faidx_t *fai_file;                //!< Fasta reference file index
	bam_index_t *idx;                 //!< Pointer to the BAM input file index
	bam_header_t *h;                  //!< Pointer to the header text for the input BAM file
	uint16_t flag;                    //!< Bit flag to hold user options
	int output;                       //!< Analysis output option
	int errorCount;                   //!< Flag to indicate error in reading user options
	int minDepth;                     //!< User-specified minimumm read depth
	int maxDepth;                     //!< User-specified maximum read depth
	int minRMSQ;                      //!< User-specified minimum rms mapping quality
	int minSNPQ;                      //!< User-specified minimum SNP quality score
	uint32_t winSize;                 //!< User-specified window size in kilobases
	uint8_t minMapQ;                  //!< User-specified minimum individual read mapping quality
	uint8_t minBaseQ;                 //!< User-specified minimum inidividual base quality
	double minSites;                  //!< User-specified minimum number of aligned sites to perform analysis
	double minPop;                    //!< Minimum proportion of samples present
	double hetPrior;                  //!< Prior probability for calling heterozygous genotypes
	std::string dist;                 //!< Pointer to the name of the desired distance metric	(-d switch)
	std::string bamfile;              //!< File name for the input BAM file
	std::string reffile;              //!< File name for the input reference Fasta file
	std::string headfile;             //!< File name for optional BAM header input file
	std::string region;               //!< Region on which to perform the analysis
	std::string errorMsg;             //!< String to hold any error messages
	std::string popFunc;              //!< The popbam function being invoked

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

		// member variables
		std::string bamfile;           //!< Name of bamfile used for indexing purposes
		bam_sample_t *sm;              //!< Pointer to the sample information for the input BAM file
		char *ref_base;                //!< Reference sequence string for specified region
		int tid;                       //!< Reference chromosome/scaffold identifier
		int beg;                       //!< Reference coordinate of the beginning of the current region
		int end;                       //!< Reference coordinate of the end of current region
		int len;                       //!< Length of the reference sequence for current region
		uint16_t flag;                 //!< Bit flag to hold user options
		int num_sites;                 //!< Total number of aligned sites
		int segsites;                  //!< Total number of segregating sites in entire sample
		uint8_t *pop_nsmpl;            //!< Sample size per population
		uint64_t *types;               //!< The site type for each aligned site
		uint64_t *pop_mask;            //!< Bit mask for which individuals are in which population
		int minDepth;                  //!< User-specified minimumm read depth
		int maxDepth;                  //!< User-specified maximum read depth
		int minRMSQ;                   //!< User-specified minimum rms mapping quality
		int minSNPQ;                   //!< User-specified minimum SNP quality score
		uint8_t minMapQ;               //!< User-specified minimum individual read mapping quality
		uint8_t minBaseQ;              //!< User-specified minimum inidividual base quality
		double hetPrior;               //!< Prior probability of heterozygous genotype
		errmod_t *em;                  //!< Error model data structure
		popbam_func_t derived_type;    //!< Type of the derived class
};

#endif
