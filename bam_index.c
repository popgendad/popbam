#include <ctype.h>
#include <assert.h>
#include "bam.h"
#include "khash.h"
#include "ksort.h"
#include "bam_endian.h"

#ifdef _MSC_VER
#define inline __inline
#endif

#define BAM_MIN_CHUNK_GAP 32768
// 1<<14 is the size of minimum bin.
#define BAM_LIDX_SHIFT    14
#define BAM_MAX_BIN 37450 // =(8^6-1)/7+1

typedef struct
{
	uint64_t u;
	uint64_t v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)
KSORT_INIT(off, pair64_t, pair64_lt)

typedef struct
{
	uint32_t m;
	uint32_t n;
	pair64_t *list;
} bam_binlist_t;

typedef struct
{
	int32_t n;
	int32_t m;
	uint64_t *offset;
} bam_lidx_t;

KHASH_MAP_INIT_INT(i, bam_binlist_t)

struct __bam_index_t
{
	int32_t n;
	uint64_t n_no_coor;           // unmapped reads without coordinate
	khash_t(i) **index;
	bam_lidx_t *index2;
};


struct __bam_iter_t
{
	int from_first;              // read from the first record; no random access
	int tid;
	int beg;
	int end;
	int n_off;
	int i;
	int finished;
	uint64_t curr_off;
	pair64_t *off;
};

void bam_index_destroy(bam_index_t *idx)
{
	khint_t k;
	int i;

	if (idx == 0)
		return;

	for (i=0; i < idx->n; ++i)
	{
		khash_t(i) *index = idx->index[i];
		bam_lidx_t *index2 = idx->index2 + i;
		for (k=kh_begin(index); k != kh_end(index); ++k)
		{
			if (kh_exist(index, k))
				free(kh_value(index, k).list);
		}
		kh_destroy(i, index);
		free(index2->offset);
	}

	free(idx->index);
	free(idx->index2);
	free(idx);
}

static bam_index_t *bam_index_load_core(FILE *fp)
{
	int i;
	char magic[4];
	bam_index_t *idx;

	if (fp == 0)
	{
		fprintf(stderr, "[bam_index_load_core] fail to load index.\n");
		return 0;
	}

	fread(magic, 1, 4, fp);

	if (strncmp(magic, "BAI\1", 4))
	{
		fprintf(stderr, "[bam_index_load] wrong magic number.\n");
		fclose(fp);
		return 0;
	}

	idx = (bam_index_t*)calloc(1, sizeof(bam_index_t));
	fread(&idx->n, 4, 1, fp);

	if (bam_is_be)
		bam_swap_endian_4p(&idx->n);

	idx->index = (khash_t(i)**)calloc(idx->n, sizeof(void*));
	idx->index2 = (bam_lidx_t*)calloc(idx->n, sizeof(bam_lidx_t));

	for (i=0; i < idx->n; ++i)
	{
		khash_t(i) *index;
		bam_lidx_t *index2 = idx->index2 + i;
		uint32_t key;
		uint32_t size;
		khint_t k;
		int j;
		int ret;
		bam_binlist_t *p;

		index = idx->index[i] = kh_init(i);

		// load binning index
		fread(&size, 4, 1, fp);

		if (bam_is_be)
			bam_swap_endian_4p(&size);

		for (j=0; j < (int)size; ++j)
		{
			fread(&key, 4, 1, fp);

			if (bam_is_be)
				bam_swap_endian_4p(&key);

			k = kh_put(i, index, key, &ret);
			p = &kh_value(index, k);
			fread(&p->n, 4, 1, fp);

			if (bam_is_be)
				bam_swap_endian_4p(&p->n);

			p->m = p->n;
			p->list = (pair64_t*)malloc(p->m*16);
			fread(p->list, 16, p->n, fp);

			if (bam_is_be)
			{
				int x;
				for (x=0; x < p->n; ++x)
				{
					bam_swap_endian_8p(&p->list[x].u);
					bam_swap_endian_8p(&p->list[x].v);
				}
			}
		}

		// load linear index
		fread(&index2->n, 4, 1, fp);

		if (bam_is_be)
			bam_swap_endian_4p(&index2->n);

		index2->m = index2->n;
		index2->offset = (uint64_t*)calloc(index2->m, 8);
		fread(index2->offset, index2->n, 8, fp);

		if (bam_is_be)
			for (j=0; j < index2->n; ++j)
				bam_swap_endian_8p(&index2->offset[j]);
	}

	if (fread(&idx->n_no_coor, 8, 1, fp) == 0)
		idx->n_no_coor = 0;

	if (bam_is_be)
		bam_swap_endian_8p(&idx->n_no_coor);

	return idx;
}

bam_index_t *bam_index_load_local(const char *_fn)
{
	FILE *fp;
	char *fnidx;
	char *fn;

	fn = strdup(_fn);
	fnidx = (char*)calloc(strlen(fn) + 5, 1);
	strcpy(fnidx, fn);
	strcat(fnidx, ".bai");
	fp = fopen(fnidx, "rb");

	// try "{base}.bai"
	if (fp == 0)
	{
		char *s = strstr(fn, "bam");
		if (s == fn + strlen(fn) - 3)
		{
			strcpy(fnidx, fn);
			fnidx[strlen(fn)-1] = 'i';
			fp = fopen(fnidx, "rb");
		}
	}

	free(fnidx);
	free(fn);

	if (fp)
	{
		bam_index_t *idx = bam_index_load_core(fp);
		fclose(fp);
		return idx;
	}
	else
		return 0;
}

bam_index_t *bam_index_load(const char *fn)
{
	bam_index_t *idx;

	idx = bam_index_load_local(fn);
	if (idx == 0)
		fprintf(stderr, "[bam_index_load] fail to load BAM index.\n");

	return idx;
}

static inline int reg2bins(uint32_t beg, uint32_t end, uint16_t list[BAM_MAX_BIN])
{
	int i = 0;
	int k;

	if (beg >= end)
		return 0;

	if (end >= 1u << 29)
		end = 1u << 29;

	--end;
	list[i++] = 0;

	for (k=1+(beg>>26); k <= 1 + (end >> 26); ++k)
		list[i++] = k;

	for (k=9+(beg>>23); k <= 9 + (end >> 23); ++k)
		list[i++] = k;

	for (k=73+(beg>>20); k <= 73 + (end >> 20); ++k)
		list[i++] = k;

	for (k=585+(beg>>17); k <= 585 + (end >> 17); ++k)
		list[i++] = k;

	for (k=4681+(beg>>14); k <= 4681 + (end >> 14); ++k)
		list[i++] = k;

	return i;
}

static inline int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b)
{
	uint32_t rbeg = b->core.pos;
	uint32_t rend = b->core.n_cigar ? bam_calend(&b->core, bam1_cigar(b)) : b->core.pos + 1;

	return ((rend > beg) && (rbeg < end));
}

// bam_fetch helper function retrieves
bam_iter_t bam_iter_query(const bam_index_t *idx, int tid, int beg, int end)
{
	uint16_t *bins;
	int i;
	int n_bins;
	int n_off;
	pair64_t *off;
	khint_t k;
	khash_t(i) *index;
	uint64_t min_off;
	bam_iter_t iter = 0;

	if (beg < 0)
		beg = 0;

	if (end < beg)
		return 0;

	// initialize iter
	iter = (bam_iter_t)calloc(1, sizeof(struct __bam_iter_t));
	iter->tid = tid, iter->beg = beg, iter->end = end;
	iter->i = -1;
	//
	bins = (uint16_t*)calloc(BAM_MAX_BIN, 2);
	n_bins = reg2bins(beg, end, bins);
	index = idx->index[tid];

	if (idx->index2[tid].n > 0)
	{
		min_off = ((beg >> BAM_LIDX_SHIFT) >= idx->index2[tid].n) ? idx->index2[tid].offset[idx->index2[tid].n-1] : idx->index2[tid].offset[beg>>BAM_LIDX_SHIFT];

		// improvement for index files built by tabix prior to 0.1.4
		if (min_off == 0)
		{
			int n = beg >> BAM_LIDX_SHIFT;

			if (n > idx->index2[tid].n)
				n = idx->index2[tid].n;

			for (i=n-1; i >= 0; --i)
				if (idx->index2[tid].offset[i] != 0)
					break;

			if (i >= 0)
				min_off = idx->index2[tid].offset[i];
		}
	}
	// tabix 0.1.2 may produce such index files
	else
		min_off = 0;

	for (i=n_off = 0; i < n_bins; ++i)
	{
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
			n_off += kh_value(index, k).n;
	}

	if (n_off == 0)
	{
		free(bins);
		return iter;
	}

	off = (pair64_t*)calloc(n_off, 16);
	for (i=n_off = 0; i < n_bins; ++i)
	{
		if ((k = kh_get(i, index, bins[i])) != kh_end(index))
		{
			int j;
			bam_binlist_t *p = &kh_value(index, k);
			for (j=0; j < p->n; ++j)
				if (p->list[j].v > min_off)
					off[n_off++] = p->list[j];
		}
	}

	free(bins);

	if (n_off == 0)
	{
		free(off);
		return iter;
	}
	{
		bam1_t *b = (bam1_t*)calloc(1, sizeof(bam1_t));
		int l;
		ks_introsort(off, n_off, off);

		// resolve completely contained adjacent blocks
		for (i=1, l=0; i < n_off; ++i)
			if (off[l].v < off[i].v)
				off[++l] = off[i];

		n_off = l + 1;

		// resolve overlaps between adjacent blocks
		// this may happen due to the merge in indexing
		for (i = 1; i < n_off; ++i)
			if (off[i-1].v >= off[i].u)
				off[i-1].v = off[i].u;
		{
			// merge adjacent blocks
#if defined(BAM_TRUE_OFFSET) || defined(BAM_VIRTUAL_OFFSET16)
			for (i=1, l=0; i < n_off; ++i)
			{
#ifdef BAM_TRUE_OFFSET
				if (off[l].v + BAM_MIN_CHUNK_GAP > off[i].u)
					off[l].v = off[i].v;
#else
				if ((off[l].v >> 16) == (off[i].u >> 16))
					off[l].v = off[i].v;
#endif
				else
					off[++l] = off[i];
			}
			n_off = l + 1;
#endif
		}
		bam_destroy1(b);
	}

	iter->n_off = n_off;
	iter->off = off;

	return iter;
}

void bam_iter_destroy(bam_iter_t iter)
{
	if (iter)
	{
		free(iter->off);
		free(iter);
	}
}

int bam_iter_read(bamFile fp, bam_iter_t iter, bam1_t *b)
{
	int ret;

	if (iter && iter->finished)
		return -1;

	if ((iter == 0) || iter->from_first)
	{
		ret = bam_read1(fp, b);

		if (ret < 0 && iter)
			iter->finished = 1;

		return ret;
	}

	if (iter->off == 0)
		return -1;

	for (;;)
	{
		// then jump to the next chunk
		if ((iter->curr_off == 0) || (iter->curr_off >= iter->off[iter->i].v))
		{
			// no more chunks
			if (iter->i == iter->n_off - 1)
			{
				ret = -1;
				break;
			}

			// otherwise bug
			if (iter->i >= 0)
				assert(iter->curr_off == iter->off[iter->i].v);

			// not adjacent chunks; then seek
			if ((iter->i < 0) || (iter->off[iter->i].v != iter->off[iter->i+1].u))
			{
				bam_seek(fp, iter->off[iter->i+1].u, SEEK_SET);
				iter->curr_off = bam_tell(fp);
			}

			++iter->i;
		}

		if ((ret = bam_read1(fp, b)) >= 0)
		{
			iter->curr_off = bam_tell(fp);

			if ((b->core.tid != iter->tid) || (b->core.pos >= iter->end))
			{
				// no need to proceed
				// determine whether end of region or error
				ret = bam_validate1(NULL, b) ? -1 : -5;
				break;
			}
			else if (is_overlap(iter->beg, iter->end, b))
				return ret;
		}
		// end of file or error
		else
			break;
	}

	iter->finished = 1;

	return ret;
}

int bam_fetch(bamFile fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)
{
	int ret;
	bam_iter_t iter;
	bam1_t *b;

	b = bam_init1();
	iter = bam_iter_query(idx, tid, beg, end);

	while ((ret = bam_iter_read(fp, iter, b)) >= 0)
		func(b, data);

	bam_iter_destroy(iter);
	bam_destroy1(b);

	return (ret == -1) ? 0 : ret;
}
