#include <zlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "tables.h"
#include "kstring.h"
#include "bam.h"
#include "sam_header.h"
#include "kseq.h"
#include "khash.h"

typedef char *str_p;

KHASH_MAP_INIT_STR(s, int)
KHASH_MAP_INIT_STR(r2l, str_p)
KSTREAM_INIT(gzFile, gzread, 16384)
KHASH_MAP_INIT_STR(ref, unsigned long long)

void bam_init_header_hash2(bam_header_t *header);
extern void bam_destroy_header_hash(bam_header_t *header);
int bam_get_tid(const bam_header_t *header, const char *seq_name);

unsigned short bam_char2flag_table[256] =
{
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,BAM_FREAD1,BAM_FREAD2,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    BAM_FPROPER_PAIR,0,BAM_FMREVERSE,0, 0,BAM_FMUNMAP,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, BAM_FDUP,0,BAM_FQCFAIL,0, 0,0,0,0, 0,0,0,0,
    BAM_FPAIRED,0,BAM_FREVERSE,BAM_FSECONDARY, 0,BAM_FUNMAP,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

struct __tamFile_t
{
    gzFile fp;
    kstream_t *ks;
    kstring_t *str;
    unsigned long long n_lines;
    int is_first;
};

static bam_header_t *
hash2header(const kh_ref_t *hash)
{
    int i = 0;
    bam_header_t *header;
    khiter_t k;

    header = bam_header_init();
    header->n_targets = kh_size(hash);
    header->target_name = (char**)calloc(kh_size(hash), sizeof(char*));
    header->target_len = (unsigned int*)calloc(kh_size(hash), 4);
    for (k = kh_begin(hash); k != kh_end(hash); ++k)
        {
            if (kh_exist(hash, k))
                {
                    i = (int)kh_value(hash, k);
                    header->target_name[i] = (char*)kh_key(hash, k);
                    header->target_len[i] = kh_value(hash, k) >> 32;
                }
        }
    bam_init_header_hash2(header);
    return header;
}

bam_header_t *
sam_header_read2(const char *fn)
{
    int i = 0;
    int c = 0;
    int dret = 0;
    int ret = 0;
    int len = 0;
    int error = 0;
    char *s = NULL;
    gzFile fp;
    kstream_t *ks;
    kstring_t *str;
    kh_ref_t *hash;
    khiter_t k;
    bam_header_t *header;

    if (fn == 0)
        {
            return 0;
        }
    fp = (strcmp(fn, "-") == 0) ? gzdopen(fileno(stdin), "r") : gzopen(fn, "r");
    if (fp == 0)
        {
            return 0;
        }
    hash = kh_init(ref);
    ks = ks_init(fp);
    str = (kstring_t*)calloc(1, sizeof(kstring_t));
    while (ks_getuntil(ks, 0, str, &dret) > 0)
        {
            s = strdup(str->s);
            i = kh_size(hash);
            ks_getuntil(ks, 0, str, &dret);
            len = atoi(str->s);
            k = kh_put(ref, hash, s, &ret);
            if (ret == 0)
                {
                    fprintf(stderr, "[sam_header_read2] duplicated sequence name: %s\n", s);
                    error = 1;
                }
            kh_value(hash, k) = (unsigned long long)len << 32 | i;
            if (dret != '\n')
                {
                    while ((c = ks_getc(ks)) != '\n' && c != -1);
                }
        }
    ks_destroy(ks);
    gzclose(fp);
    free(str->s);
    free(str);
    fprintf(stderr, "[sam_header_read2] %d sequences loaded.\n", kh_size(hash));
    if (error)
        {
            return 0;
        }
    header = hash2header(hash);
    kh_destroy(ref, hash);
    return header;
}

static inline unsigned char *
alloc_data(bam1_t *b, int size)
{
    if (b->m_data < size)
        {
            b->m_data = size;
            kroundup32(b->m_data);
            b->data = (unsigned char*)realloc(b->data, b->m_data);
        }
    return b->data;
}

static inline void
parse_error(long long n_lines, const char * __restrict msg)
{
    fprintf(stderr, "Parse error at line %lld: %s\n", (long long)n_lines, msg);
    abort();
}

static inline void
append_text(bam_header_t *header, kstring_t *str)
{
    // 2 = 1 byte dret + 1 byte null
    size_t x = header->l_text, y = header->l_text + str->l + 2;
    kroundup32(x);
    kroundup32(y);
    if (x < y)
        {
            header->n_text = y;
            header->text = (char*)realloc(header->text, y);

            if (!header->text)
                {
                    fprintf(stderr,"realloc failed to alloc %ld bytes\n", y);
                    abort();
                }
        }

    // Sanity check
    if (header->l_text+str->l+1 >= header->n_text)
        {
            fprintf(stderr,"append_text FIXME: %ld>=%ld, x=%ld,y=%ld\n",  header->l_text+str->l+1,(long)header->n_text,x,y);
            abort();
        }

    // we cannot use strcpy() here.
    strncpy(header->text + header->l_text, str->s, str->l+1);
    header->l_text += str->l + 1;
    header->text[header->l_text] = 0;
}

int
sam_header_parse(bam_header_t *h)
{
    int i = 0;
    char **tmp = NULL;

    free(h->target_len);
    free(h->target_name);

    h->n_targets = 0;
    h->target_len = 0;
    h->target_name = 0;
    if (h->l_text < 3)
        {
            return 0;
        }
    if (h->dict == 0)
        {
            h->dict = sam_header_parse2(h->text);
        }
    tmp = sam_header2list(h->dict, "SQ", "SN", &h->n_targets);
    if (h->n_targets == 0)
        {
            return 0;
        }
    h->target_name = (char**)calloc(h->n_targets, sizeof(char*));
    for (i = 0; i < h->n_targets; ++i)
        {
            h->target_name[i] = strdup(tmp[i]);
        }
    free(tmp);
    tmp = sam_header2list(h->dict, "SQ", "LN", &h->n_targets);
    h->target_len = (unsigned int*)calloc(h->n_targets, 4);
    for (i = 0; i < h->n_targets; ++i)
        {
            h->target_len[i] = atoi(tmp[i]);
        }
    free(tmp);
    return h->n_targets;
}

void
bam_init_header_hash2(bam_header_t *header)
{
    int i = 0;

    if (header->hash == 0)
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

bam_header_t *
sam_header_read(tamFile fp)
{
    int ret = 0;
    int dret = 0;
    bam_header_t *header = bam_header_init();
    kstring_t *str = fp->str;

    while (((ret = ks_getuntil(fp->ks, KS_SEP_TAB, str, &dret)) >= 0) && (str->s[0] == '@'))
        {
            // skip header
            // note that str->s is NOT null terminated!!
            str->s[str->l] = dret;
            append_text(header, str);
            if (dret != '\n')
                {
                    ret = ks_getuntil(fp->ks, '\n', str, &dret);

                    // NOT null terminated!!
                    str->s[str->l] = '\n';
                    append_text(header, str);
                }
            ++fp->n_lines;
        }
    sam_header_parse(header);
    bam_init_header_hash2(header);
    fp->is_first = 1;
    return header;
}

tamFile
sam_open(const char *fn)
{
    tamFile fp;
    gzFile gzfp = (strcmp(fn, "-") == 0) ? gzdopen(fileno(stdin), "rb") : gzopen(fn, "rb");

    if (gzfp == 0)
        {
            return 0;
        }
    fp = (tamFile)calloc(1, sizeof(struct __tamFile_t));
    fp->str = (kstring_t*)calloc(1, sizeof(kstring_t));
    fp->fp = gzfp;
    fp->ks = ks_init(fp->fp);
    return fp;
}

void
sam_close(tamFile fp)
{
    if (fp)
        {
            ks_destroy(fp->ks);
            gzclose(fp->fp);
            free(fp->str->s);
            free(fp->str);
            free(fp);
        }
}
