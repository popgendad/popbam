/** \file bam.h
 *  \brief Header for the bam.c file
 *  \author Daniel Garrigan
 *  \version 0.5
 *  Much of the code is adapted from samtools written by Heng Li
*/

#ifndef BAM_BAM_H
#define BAM_BAM_H

/*!
  @header

  BAM library provides I/O and various operations on manipulating files
  in the BAM (Binary Alignment/Mapping) or SAM (Sequence Alignment/Map)
  format. It now supports importing from or exporting to SAM, sorting,
  merging, generating pileup, and quickly retrieval of reads overlapped
  with a specified region.

 */
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#ifndef BAM_LITE
#define BAM_VIRTUAL_OFFSET16
#include "bgzf.h"
/*! @abstract BAM file handler */
typedef BGZF *bamFile;
#define bam_open(fn, mode) bgzf_open(fn, mode)
#define bam_dopen(fd, mode) bgzf_fdopen(fd, mode)
#define bam_close(fp) bgzf_close(fp)
#define bam_read(fp, buf, size) bgzf_read(fp, buf, size)
#define bam_write(fp, buf, size) bgzf_write(fp, buf, size)
#define bam_tell(fp) bgzf_tell(fp)
#define bam_seek(fp, pos, dir) bgzf_seek(fp, pos, dir)
#else
#define BAM_TRUE_OFFSET
#include <zlib.h>
typedef gzFile bamFile;
#define bam_open(fn, mode) gzopen(fn, mode)
#define bam_dopen(fd, mode) gzdopen(fd, mode)
#define bam_close(fp) gzclose(fp)
#define bam_read(fp, buf, size) gzread(fp, buf, size)
/* no bam_write/bam_tell/bam_seek() here */
#endif

/*! @typedef
  @abstract Structure for the alignment header.
  @field n_targets   number of reference sequences
  @field target_name names of the reference sequences
  @field target_len  lengths of the referene sequences
  @field dict        header dictionary
  @field hash        hash table for fast name lookup
  @field rg2lib      hash table for @RG-ID -> LB lookup
  @field l_text      length of the plain text in the header
  @field text        plain text

  @discussion Field hash points to null by default. It is a private
  member.
 */
typedef struct
{
	int n_targets;
	char **target_name;
	unsigned int *target_len;
	void *dict;
	void *hash;
	void *rg2lib;
	unsigned int l_text;
	unsigned int n_text;
	char *text;
} bam_header_t;

/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024

#define BAM_OFDEC          0
#define BAM_OFHEX          1
#define BAM_OFSTR          2

/*! @abstract defautl mask for pileup */
#define BAM_DEF_MASK (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

#define BAM_CORE_SIZE   sizeof(bam1_core_t)

/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

/*
  CIGAR operations.
 */
/*! @abstract CIGAR: M = match or mismatch*/
#define BAM_CMATCH      0
/*! @abstract CIGAR: I = insertion to the reference */
#define BAM_CINS        1
/*! @abstract CIGAR: D = deletion from the reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: N = skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: S = clip on the read with clipped sequence
  present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: H = clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: P = padding */
#define BAM_CPAD        6
/*! @abstract CIGAR: equals = match */
#define BAM_CEQUAL      7
/*! @abstract CIGAR: X = mismatch */
#define BAM_CDIFF       8
#define BAM_CBACK       9

#define BAM_CIGAR_STR  "MIDNSHP=XB"
#define BAM_CIGAR_TYPE 0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference

/*! @typedef
  @abstract Structure for core alignment information.
  @field  tid     chromosome ID, defined by bam_header_t
  @field  pos     0-based leftmost coordinate
  @field  strand  strand; 0 for forward and 1 otherwise
  @field  bin     bin calculated by bam_reg2bin()
  @field  qual    mapping quality
  @field  l_qname length of the query name
  @field  flag    bitwise flag
  @field  n_cigar number of CIGAR operations
  @field  l_qseq  length of the query sequence (read)
 */
typedef struct {
	int tid;
	int pos;
	unsigned int bin:16, qual:8, l_qname:8;
	unsigned int flag:16, n_cigar:16;
	int l_qseq;
	int mtid;
	int mpos;
	int isize;
} bam1_core_t;

/*! @typedef
  @abstract Structure for one alignment.
  @field  core       core information about the alignment
  @field  l_aux      length of auxiliary data
  @field  data_len   current length of bam1_t::data
  @field  m_data     maximum length of bam1_t::data
  @field  data       all variable-length data, concatenated; structure: cigar-qname-seq-qual-aux

  @discussion Notes:
 
   1. qname is zero tailing and core.l_qname includes the tailing '\0'.
   2. l_qseq is calculated from the total length of an alignment block
	  on reading or from CIGAR.
 */
typedef struct {
	bam1_core_t core;
	int l_aux, data_len, m_data;
	unsigned char *data;
} bam1_t;

typedef struct __bam_iter_t *bam_iter_t;

#define bam1_strand(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define bam1_mstrand(b) (((b)->core.flag&BAM_FMREVERSE) != 0)

/*! @function
  @abstract  Get the CIGAR array
  @param  b  pointer to an alignment
  @return    pointer to the CIGAR array

  @discussion In the CIGAR array, each element is a 32-bit integer. The
  lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
  length of a CIGAR.
 */
#define bam1_cigar(b) ((unsigned int*)((b)->data + (b)->core.l_qname))

/*! @function
  @abstract  Get the name of the query
  @param  b  pointer to an alignment
  @return    pointer to the name string, null terminated
 */
#define bam1_qname(b) ((char*)((b)->data))

/*! @function
  @abstract  Get query sequence
  @param  b  pointer to an alignment
  @return    pointer to sequence

  @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
  8 for T and 15 for N. Two bases are packed in one byte with the base
  at the higher 4 bits having smaller coordinate on the read. It is
  recommended to use bam1_seqi() macro to get the base.
 */
#define bam1_seq(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname)

/*! @function
  @abstract  Get query quality
  @param  b  pointer to an alignment
  @return    pointer to quality string
 */
#define bam1_qual(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))

/*! @function
  @abstract  Get a base on read
  @param  s  Query sequence returned by bam1_seq()
  @param  i  The i-th position, 0-based
  @return    4-bit integer representing the base.
 */
//#define bam1_seqi(s, i) ((s)[(i)/2] >> 4*(1-(i)%2) & 0xf)
#define bam1_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)

#define bam1_seq_seti(s, i, c) ( (s)[(i)>>1] = ((s)[(i)>>1] & 0xf<<(((i)&1)<<2)) | (c)<<((~(i)&1)<<2) )

/*! @function
  @abstract  Get query sequence and quality
  @param  b  pointer to an alignment
  @return    pointer to the concatenated auxiliary data
 */
#define bam1_aux(b) ((b)->data + (b)->core.n_cigar*4 + (b)->core.l_qname + (b)->core.l_qseq + ((b)->core.l_qseq + 1)/2)

#ifndef kroundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
 */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*!
  @abstract Whether the machine is big-endian; modified only in
  bam_header_init().
 */
extern int bam_is_be;

/*!
  @abstract Verbose level between 0 and 3; 0 is supposed to disable all
  debugging information, though this may not have been implemented.
 */
extern int bam_verbose;

extern int bam_no_B;
extern void bam_init_header_hash(bam_header_t *header);

#ifdef __cplusplus
extern "C" {
#endif

// SAM I/O functions

/*! @abstract TAM file handler */
typedef struct __tamFile_t *tamFile;

/*!
  @abstract   Open a SAM file for reading, either uncompressed or compressed by gzip/zlib.
  @param  fn  SAM file name
  @return     SAM file handler
 */
tamFile sam_open(const char *fn);

/*!
  @abstract   Close a SAM file handler
  @param  fp  SAM file handler
 */
void sam_close(tamFile fp);

/*!
  @abstract      Read one alignment from a SAM file handler
  @param  fp     SAM file handler
  @param  header header information (ordered names of chromosomes)
  @param  b      read alignment; all members in b will be updated
  @return        0 if successful; otherwise negative
 */
int sam_read1(tamFile fp, bam_header_t *header, bam1_t *b);

/*!
  @abstract       Read header information from a TAB-delimited list file.
  @param  fn_list file name for the list
  @return         a pointer to the header structure

  @discussion Each line in this file consists of chromosome name and
  the length of chromosome.
 */
bam_header_t *sam_header_read2(const char *fn_list);

/*!
  @abstract       Read header from a SAM file (if present)
  @param  fp      SAM file handler
  @return         pointer to header struct; 0 if no @SQ lines available
 */
bam_header_t *sam_header_read(tamFile fp);

/*!
  @abstract       Parse @SQ lines a update a header struct
  @param  h       pointer to the header struct to be updated
  @return         number of target sequences

  @discussion bam_header_t::{n_targets,target_len,target_name} will
  be destroyed in the first place.
 */
int sam_header_parse(bam_header_t *h);
int bam_get_tid(const bam_header_t *header, const char *seq_name);


// BAM I/O functions

/*!
  @abstract Initialize a header structure.
  @return   the pointer to the header structure

  @discussion This function also modifies the global variable
  bam_is_be.
 */
bam_header_t *bam_header_init();

/*!
  @abstract        Destroy a header structure.
  @param  header  pointer to the header
 */
void bam_header_destroy(bam_header_t *header);

/*!
  @abstract   Read a header structure from BAM.
  @param  fp  BAM file handler, opened by bam_open()
  @return     pointer to the header structure

  @discussion The file position indicator must be placed at the
  beginning of the file. Upon success, the position indicator will
  be set at the start of the first alignment.
 */
bam_header_t *bam_header_read(bamFile fp);

/*!
  @abstract      Write a header structure to BAM.
  @param  fp     BAM file handler
  @param  header pointer to the header structure
  @return        always 0 currently
 */
int bam_header_write(bamFile fp, const bam_header_t *header);

/*!
  @abstract   Read an alignment from BAM.
  @param  fp  BAM file handler
  @param  b   read alignment; all members are updated.
  @return     number of bytes read from the file

  @discussion The file position indicator must be
  placed right before an alignment. Upon success, this function
  will set the position indicator to the start of the next
  alignment. This function is not affected by the machine
  endianness.
 */
int bam_read1(bamFile fp, bam1_t *b);

int bam_remove_B(bam1_t *b);

/*! @function
	@abstract  Initiate a pointer to bam1_t struct
 */
#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))

/*! @function
  @abstract  Free the memory allocated for an alignment.
  @param  b  pointer to an alignment
 */
#define bam_destroy1(b) do {					\
		if (b) { free((b)->data); free(b); }	\
	} while (0)

/*!
  @abstract       Check whether a BAM record is plausibly valid
  @param  header  associated header structure, or NULL if unavailable
  @param  b       alignment to validate
  @return         0 if the alignment is invalid; non-zero otherwise

  @discussion  Simple consistency check of some of the fields of the
  alignment record.  If the header is provided, several additional checks
  are made.  Not all fields are checked, so a non-zero result is not a
  guarantee that the record is valid.  However it is usually good enough
  to detect when bam_seek() has been called with a virtual file offset
  that is not the offset of an alignment record.
 */
int bam_validate1(const bam_header_t *header, const bam1_t *b);

/*! @typedef
  @abstract Structure for one alignment covering the pileup position.
  @field  b      pointer to the alignment
  @field  qpos   position of the read base at the pileup site, 0-based
  @field  indel  indel length; 0 for no indel, positive for ins and negative for del
  @field  is_del 1 iff the base on the padded read is a deletion
  @field  level  the level of the read in the "viewer" mode

  @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
  difference between the two functions is that the former does not
  set bam_pileup1_t::level, while the later does. Level helps the
  implementation of alignment viewers, but calculating this has some
  overhead.
 */
typedef struct
{
	bam1_t *b;
	int qpos;
	int indel, level;
	unsigned int is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
} bam_pileup1_t;

typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);
struct __bam_plp_t;
typedef struct __bam_plp_t *bam_plp_t;

bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
int bam_plp_push(bam_plp_t iter, const bam1_t *b);
const bam_pileup1_t *bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
void bam_plp_destroy(bam_plp_t iter);

/*! @typedef
  @abstract    Type of function to be called by bam_plbuf_push().
  @param  tid  chromosome ID as is defined in the header
  @param  pos  start coordinate of the alignment, 0-based
  @param  n    number of elements in pl array
  @param  pl   array of alignments
  @param  data user provided data
  @discussion  See also bam_plbuf_push(), bam_plbuf_init() and bam_pileup1_t.
 */
typedef int (*bam_pileup_f)(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

typedef struct
{
	bam_plp_t iter;
	bam_pileup_f func;
	void *data;
} bam_plbuf_t;

void bam_plbuf_set_mask(bam_plbuf_t *buf, int mask);
void bam_plbuf_reset(bam_plbuf_t *buf);
bam_plbuf_t *bam_plbuf_init(bam_pileup_f func, void *data);
void bam_plbuf_destroy(bam_plbuf_t *buf);
int bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf);

struct __bam_lplbuf_t;
typedef struct __bam_lplbuf_t bam_lplbuf_t;

//void bam_lplbuf_reset(bam_lplbuf_t *buf);

/*! @abstract  bam_plbuf_init() equivalent with level calculated. */
bam_lplbuf_t *bam_lplbuf_init(bam_pileup_f func, void *data);

/*! @abstract  bam_plbuf_destroy() equivalent with level calculated. */
void bam_lplbuf_destroy(bam_lplbuf_t *tv);

/*! @abstract  bam_plbuf_push() equivalent with level calculated. */
int bam_lplbuf_push(const bam1_t *b, bam_lplbuf_t *buf);


// BAM indexing functions

struct __bam_index_t;
typedef struct __bam_index_t bam_index_t;

/*!
  @abstract   Build index for a BAM file.
  @discussion Index file "fn.bai" will be created.
  @param  fn  name of the BAM file
  @return     always 0 currently
 */
int bam_index_build(const char *fn);

/*!
  @abstract   Load index from file "fn.bai".
  @param  fn  name of the BAM file (NOT the index file)
  @return     pointer to the index structure
 */
bam_index_t *bam_index_load(const char *fn);

/*!
  @abstract    Destroy an index structure.
  @param  idx  pointer to the index structure
 */
void bam_index_destroy(bam_index_t *idx);

/*! @typedef
  @abstract      Type of function to be called by bam_fetch().
  @param  b     the alignment
  @param  data  user provided data
 */
typedef int (*bam_fetch_f)(const bam1_t *b, void *data);

/*!
  @abstract Retrieve the alignments that are overlapped with the
  specified region.
  @discussion A user defined function will be called for each
  retrieved alignment ordered by its start position.
  @param  fp    BAM file handler
  @param  idx   pointer to the alignment index
  @param  tid   chromosome ID as is defined in the header
  @param  beg   start coordinate, 0-based
  @param  end   end coordinate, 0-based
  @param  data  user provided data (will be transferred to func)
  @param  func  user defined function
 */
int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);

bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end);
int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b);
void bam_iter_destroy(bam_iter_t iter);

// tag handling functions

/*!
  @abstract       Retrieve data of a tag
  @param  b       pointer to an alignment struct
  @param  tag     two-character tag to be retrieved
  @return         pointer to the type and data. The first character is the
				  type that can be 'iIsScCdfAZH'.
  @discussion  Use bam_aux2?() series to convert the returned data to
			   the corresponding type.
*/
unsigned char *bam_aux_get(const bam1_t *b, const char tag[2]);

int bam_aux2i(const unsigned char *s);
float bam_aux2f(const unsigned char *s);
double bam_aux2d(const unsigned char *s);
char bam_aux2A(const unsigned char *s);
char *bam_aux2Z(const unsigned char *s);
int bam_aux_del(bam1_t *b, unsigned char *s);
void bam_aux_append(bam1_t *b, const char tag[2], char type, int len, unsigned char *data);
unsigned char *bam_aux_get_core(bam1_t *b, const char tag[2]); // an alias of bam_aux_get()


// misc functions

/*!  
  @abstract Calculate the rightmost coordinate of an alignment on the
  reference genome.

  @param  c      pointer to the bam1_core_t structure
  @param  cigar  the corresponding CIGAR array (from bam1_t::cigar)
  @return        the rightmost coordinate, 0-based
*/
unsigned int bam_calend(const bam1_core_t *c, const unsigned int *cigar);

#ifdef __cplusplus
}
#endif

#ifdef _MSC_VER
#define inline __inline
#endif

/*!
  @abstract    Calculate the minimum bin that contains a region [beg,end).
  @param  beg  start of the region, 0-based
  @param  end  end of the region, 0-based
  @return      bin
 */
static inline int bam_reg2bin(unsigned int beg, unsigned int end)
{
	--end;
	if (beg>>14 == end>>14)
		return 4681 + (beg>>14);
	if (beg>>17 == end>>17)
		return  585 + (beg>>17);
	if (beg>>20 == end>>20)
		return   73 + (beg>>20);
	if (beg>>23 == end>>23)
		return    9 + (beg>>23);
	if (beg>>26 == end>>26)
		return    1 + (beg>>26);
	return 0;
}

/*!
  @abstract     Copy an alignment
  @param  bdst  destination alignment struct
  @param  bsrc  source alignment struct
  @return       pointer to the destination alignment struct
 */
static inline bam1_t *bam_copy1(bam1_t *bdst, const bam1_t *bsrc)
{
	unsigned char *data = bdst->data;
	int m_data = bdst->m_data;

	if (m_data < bsrc->data_len)
	{
		// double the capacity
		m_data = bsrc->data_len; kroundup32(m_data);
		data = (unsigned char*)realloc(data, m_data);
	}
	
	// copy var-len data
	memcpy(data, bsrc->data, bsrc->data_len);
	
	// copy the rest
	*bdst = *bsrc; 
	
	// restore the backup
	bdst->m_data = m_data;
	bdst->data = data;

	return bdst;
}

/*!
  @abstract     Duplicate an alignment
  @param  src   source alignment struct
  @return       pointer to the destination alignment struct
 */
static inline bam1_t *bam_dup1(const bam1_t *src)
{
	bam1_t *b;

	b = bam_init1();
	*b = *src;
	b->m_data = b->data_len;
	b->data = (unsigned char*)calloc(b->data_len, 1);
	memcpy(b->data, src->data, b->data_len);

	return b;
}

static inline int bam_aux_type2size(int x)
{
	if ((x == 'C') || (x == 'c') || (x == 'A'))
		return 1;
	else if ((x == 'S') || (x == 's'))
		return 2;
	else if ((x == 'I') || (x == 'i') || (x == 'f'))
		return 4;
	else
		return 0;
}


#endif
