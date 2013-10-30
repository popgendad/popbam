#include <string.h>
#include "faidx.h"
#include "sam.h"

#define TYPE_BAM  1
#define TYPE_READ 2

#ifdef _MSC_VER
#define access _access
#include <io.h>
const int W_OK = 2;
const int R_OK = 4;
#endif

bam_header_t *bam_header_dup(const bam_header_t *h0)
{
    bam_header_t *h;
    int i;

    h = bam_header_init();
    *h = *h0;
    h->hash = h->dict = h->rg2lib = 0;
    h->text = (char*)calloc(h->l_text + 1, 1);
    memcpy(h->text, h0->text, h->l_text);
    h->target_len = (uint32_t*)calloc(h->n_targets, 4);
    h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
    for (i=0; i < h->n_targets; ++i)
    {
        h->target_len[i] = h0->target_len[i];
        h->target_name[i] = strdup(h0->target_name[i]);
    }

    return h;
}

static void append_header_text(bam_header_t *header, char* text, int len)
{
    int x = header->l_text + 1;
    int y = header->l_text + len + 1;

    if (text == 0)
        return;
    kroundup32(x);
    kroundup32(y);
    if (x < y)
        header->text = (char*)realloc(header->text, y);
    // we cannot use strcpy() here.
    strncpy(header->text + header->l_text, text, len);
    header->l_text += len;
    header->text[header->l_text] = 0;
}

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
    samfile_t *fp;

    fp = (samfile_t*)calloc(1, sizeof(samfile_t));

    if (strchr(mode, 'r'))
    {
        fp->type |= TYPE_READ;
        if (strchr(mode, 'b'))
        {
            fp->type |= TYPE_BAM;
            fp->x.bam = strcmp(fn, "-") ? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
            if (fp->x.bam == 0)
                goto open_err_ret;
            fp->header = bam_header_read(fp->x.bam);
        }
        else
        {
            fp->x.tamr = sam_open(fn);
            if (fp->x.tamr == 0)
                goto open_err_ret;
            fp->header = sam_header_read(fp->x.tamr);
            if (fp->header->n_targets == 0)
            {
                // no @SQ fields
                if (aux)
                {
                    // check if aux is present
                    bam_header_t *textheader = fp->header;
                    fp->header = sam_header_read2((const char*)aux);
                    if (fp->header == 0)
                        goto open_err_ret;
                    append_header_text(fp->header, textheader->text, textheader->l_text);
                    bam_header_destroy(textheader);
                }
                if ((fp->header->n_targets == 0) && (bam_verbose >= 1))
                    fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
            }
            else if (bam_verbose >= 2)
                fprintf(stderr, "[samopen] SAM header is present: %d sequences.\n", fp->header->n_targets);
        }
    }
    else if (strchr(mode, 'w'))
    {
        // write
        fp->header = bam_header_dup((const bam_header_t*)aux);
        if (strchr(mode, 'b'))
        {
            // binary
            char bmode[3];
            int i;
            int compress_level = -1;

            for (i=0; mode[i]; ++i)
                if ((mode[i] >= '0') && (mode[i] <= '9'))
                    break;
            if (mode[i])
                compress_level = mode[i] - '0';
            if (strchr(mode, 'u'))
                compress_level = 0;
            bmode[0] = 'w';
            bmode[1] = compress_level < 0 ? 0 : compress_level + '0';
            bmode[2] = 0;
            fp->type |= TYPE_BAM;
            fp->x.bam = strcmp(fn, "-") ? bam_open(fn, bmode) : bam_dopen(fileno(stdout), bmode);
            if (fp->x.bam == 0)
                goto open_err_ret;
            bam_header_write(fp->x.bam, fp->header);
        }
        else
        {
            // text
            // open file
            fp->x.tamw = strcmp(fn, "-") ? fopen(fn, "w") : stdout;
            if (fp->x.tamw == 0)
                goto open_err_ret;
            if (strchr(mode, 'X'))
                fp->type |= BAM_OFSTR << 2;
            else if (strchr(mode, 'x'))
                fp->type |= BAM_OFHEX << 2;
            else
                fp->type |= BAM_OFDEC << 2;
            // write header
            if (strchr(mode, 'h'))
            {
                int i;
                bam_header_t *alt;

                // parse the header text
                alt = bam_header_init();
                alt->l_text = fp->header->l_text;
                alt->text = fp->header->text;
                sam_header_parse(alt);
                alt->l_text = 0;
                alt->text = 0;

                // check if there are @SQ lines in the header
                // FIXME: better to skip the trailing NULL
                fwrite(fp->header->text, 1, fp->header->l_text, fp->x.tamw);
                if (alt->n_targets)
                {
                    // then write the header text without dumping ->target_{name,len}
                    if ((alt->n_targets != fp->header->n_targets) && (bam_verbose >= 1))
                        fprintf(stderr, "[samopen] inconsistent number of target sequences. Output the text header.\n");
                }
                else
                {
                    // then dump ->target_{name,len}
                    for (i=0; i < fp->header->n_targets; ++i)
                        fprintf(fp->x.tamw, "@SQ\tSN:%s\tLN:%d\n", fp->header->target_name[i], fp->header->target_len[i]);
                }
                bam_header_destroy(alt);
            }
        }
    }

    return fp;

open_err_ret:
    free(fp);
    return 0;
}

void samclose(samfile_t *fp)
{
    if (fp == 0)
        return;
    if (fp->header)
        bam_header_destroy(fp->header);
    if (fp->type & TYPE_BAM)
        bam_close(fp->x.bam);
    else if (fp->type & TYPE_READ)
        sam_close(fp->x.tamr);
    else
        fclose(fp->x.tamw);
    free(fp);
}
