#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include "bam.h"
#include "bam_endian.h"
#include "kstring.h"
#include "sam_header.h"

int bam_is_be = 0;
int bam_verbose = 2;
int bam_no_B = 0;
char *bam_flag2char_table = "pPuUrR12sfd\0\0\0\0\0";
char bam_nt16_rev_table2[16] = {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

unsigned int
bam_calend(const bam1_core_t *c, const unsigned int *cigar)
{
    int k = 0;
    int end = c->pos;

    for (k = 0; k < c->n_cigar; ++k)
        {
            int op = 0;
            int len = 0;

            op  = bam_cigar_op(cigar[k]);
            len = bam_cigar_oplen(cigar[k]);
            if (op == BAM_CBACK)
                {
                    // move backward
                    int l = 0;
                    int u = 0;
                    int v = 0;

                    // skip trailing 'B'
                    if (k == c->n_cigar-1)
                        {
                            break;
                        }
                    for (l = k - 1, u = v = 0; l >= 0; --l)
                        {
                            int op1;
                            int len1;

                            op1  = bam_cigar_op(cigar[l]);
                            len1 = bam_cigar_oplen(cigar[l]);
                            if (bam_cigar_type(op1) & 1)
                                {
                                    // consume query
                                    if ((u + len1) >= len)
                                        {
                                            // stop
                                            if (bam_cigar_type(op1) & 2)
                                                {
                                                    v += len - u;
                                                }
                                            break;
                                        }
                                    else
                                        {
                                            u += len1;
                                        }
                                }
                            if (bam_cigar_type(op1) & 2)
                                {
                                    v += len1;
                                }
                        }
                    end = l < 0 ? c->pos : end-v;
                }
            else if (bam_cigar_type(op) & 2)
                {
                    end += bam_cigar_oplen(cigar[k]);
                }
        }
    return end;
}

bam_header_t *
bam_header_init(void)
{
    bam_is_be = bam_is_big_endian();
    return (bam_header_t*)calloc(1, sizeof(bam_header_t));
}

void
bam_header_destroy(bam_header_t *header)
{
    int i = 0;
    extern void bam_destroy_header_hash(bam_header_t *header);

    if (header == 0)
        {
            return;
        }
    if (header->target_name)
        {
            for (i = 0; i < header->n_targets; ++i)
                {
                    free(header->target_name[i]);
                }
            free(header->target_name);
            free(header->target_len);
        }
    free(header->text);
    if (header->dict)
        {
            sam_header_free(header->dict);
        }
    if (header->rg2lib)
        {
            sam_tbl_destroy(header->rg2lib);
        }
    bam_destroy_header_hash(header);
    free(header);
}

bam_header_t *
bam_header_read(bamFile fp)
{
    int magic_len = 0;
    int i = 1;
    int name_len = 0;
    char buf[4];
    bam_header_t *header;

    // check EOF
    i = bgzf_check_EOF(fp);

    if (i < 0)
        {
            // If the file is a pipe, checking the EOF marker will *always* fail
            // with ESPIPE.  Suppress the error message in this case.
            if (errno != ESPIPE)
                {
                    perror("[bam_header_read] bgzf_check_EOF");
                }
        }
    else if (i == 0)
        {
            fprintf(stderr, "[bam_header_read] EOF marker is absent. The input is probably truncated.\n");
        }

    // read "BAM1"
    magic_len = bam_read(fp, buf, 4);

    if ((magic_len != 4) || (strncmp(buf, "BAM\001", 4) != 0))
        {
            fprintf(stderr, "[bam_header_read] invalid BAM binary header (this is not a BAM file).\n");
            return 0;
        }
    header = bam_header_init();

    // read plain text and the number of reference sequences
    bam_read(fp, &header->l_text, 4);

    if (bam_is_be)
        {
            bam_swap_endian_4p(&header->l_text);
        }
    header->text = (char*)calloc(header->l_text + 1, 1);
    bam_read(fp, header->text, header->l_text);
    bam_read(fp, &header->n_targets, 4);
    if (bam_is_be)
        {
            bam_swap_endian_4p(&header->n_targets);
        }

    // read reference sequence names and lengths
    header->target_name = (char**)calloc(header->n_targets, sizeof(char*));
    header->target_len = (unsigned int*)calloc(header->n_targets, 4);

    for (i = 0; i != header->n_targets; ++i)
        {
            bam_read(fp, &name_len, 4);
            if (bam_is_be)
                {
                    bam_swap_endian_4p(&name_len);
                }
            header->target_name[i] = (char*)calloc(name_len, 1);
            bam_read(fp, header->target_name[i], name_len);
            bam_read(fp, &header->target_len[i], 4);
            if (bam_is_be)
                {
                    bam_swap_endian_4p(&header->target_len[i]);
                }
        }
    return header;
}

int
bam_header_write(bamFile fp, const bam_header_t *header)
{
    int i = 0;
    int name_len = 0;
    int x = 0;
    char buf[4];

    // write "BAM1"
    strncpy(buf, "BAM\001", 4);
    bam_write(fp, buf, 4);

    // write plain text and the number of reference sequences
    if (bam_is_be)
        {
            x = bam_swap_endian_4(header->l_text);
            bam_write(fp, &x, 4);
            if (header->l_text)
                {
                    bam_write(fp, header->text, header->l_text);
                }
            x = bam_swap_endian_4(header->n_targets);
            bam_write(fp, &x, 4);
        }
    else
        {
            bam_write(fp, &header->l_text, 4);
            if (header->l_text)
                {
                    bam_write(fp, header->text, header->l_text);
                }
            bam_write(fp, &header->n_targets, 4);
        }

    // write sequence names and lengths
    for (i = 0; i != header->n_targets; ++i)
        {
            char *p = header->target_name[i];

            name_len = strlen(p)+1;
            if (bam_is_be)
                {
                    x = bam_swap_endian_4(name_len);
                    bam_write(fp, &x, 4);
                }
            else
                {
                    bam_write(fp, &name_len, 4);
                }
            bam_write(fp, p, name_len);
            if (bam_is_be)
                {
                    x = bam_swap_endian_4(header->target_len[i]);
                    bam_write(fp, &x, 4);
                }
            else
                {
                    bam_write(fp, &header->target_len[i], 4);
                }
        }
    bgzf_flush(fp);
    return 0;
}

static void
swap_endian_data(const bam1_core_t *c, int data_len, unsigned char *data)
{
    unsigned int i = 0;
    unsigned int *cigar = (unsigned int*)(data + c->l_qname);
    unsigned char *s;

    s = data + c->n_cigar * 4 + c->l_qname + c->l_qseq + (c->l_qseq + 1) / 2;
    for (i = 0; i < c->n_cigar; ++i)
        {
            bam_swap_endian_4p(&cigar[i]);
        }
    while (s < data + data_len)
        {
            unsigned char type = 0;

            // skip key
            s += 2;
            type = toupper(*s);

            // skip type
            ++s;

            if (type == 'C' || type == 'A')
                {
                    ++s;
                }
            else if (type == 'S')
                {
                    bam_swap_endian_2p(s);
                    s += 2;
                }
            else if ((type == 'I') || (type == 'F'))
                {
                    bam_swap_endian_4p(s);
                    s += 4;
                }
            else if (type == 'D')
                {
                    bam_swap_endian_8p(s);
                    s += 8;
                }
            else if ((type == 'Z') || (type == 'H'))
                {
                    while (*s)
                        {
                            ++s;
                        }
                    ++s;
                }
            else if (type == 'B')
                {
                    int n = 0;
                    int Bsize = bam_aux_type2size(*s);

                    memcpy(&n, s + 1, 4);
                    if (1 == Bsize)
                        {
                        }
                    else if (2 == Bsize)
                        {
                            for (i = 0; i < n; i += 2)
                                {
                                    bam_swap_endian_2p(s + 5 + i);
                                }
                        }
                    else if (4 == Bsize)
                        {
                            for (i = 0; i < n; i += 4)
                                {
                                    bam_swap_endian_4p(s + 5 + i);
                                }
                        }
                    bam_swap_endian_4p(s + 1);
                }
        }
}

int
bam_read1(bamFile fp, bam1_t *b)
{
    int block_len = 0;
    int ret = 0;
    int i = 0;
    unsigned int x[8];
    bam1_core_t *c = &b->core;

    assert(BAM_CORE_SIZE == 32);
    if ((ret = bam_read(fp, &block_len, 4)) != 4)
        {
            if (ret == 0)
                {
                    return -1;
                }
            else
                {
                    return -2;
                }
        }
    if (bam_read(fp, x, BAM_CORE_SIZE) != BAM_CORE_SIZE)
        {
            return -3;
        }
    if (bam_is_be)
        {
            bam_swap_endian_4p(&block_len);
            for (i = 0; i < 8; ++i)
                {
                    bam_swap_endian_4p(x + i);
                }
        }
    c->tid = x[0];
    c->pos = x[1];
    c->bin = x[2] >> 16;
    c->qual = x[2] >> 8 & 0xff;
    c->l_qname = x[2] & 0xff;
    c->flag = x[3] >> 16;
    c->n_cigar = x[3] & 0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5];
    c->mpos = x[6];
    c->isize = x[7];
    b->data_len = block_len - BAM_CORE_SIZE;
    if (b->m_data < b->data_len)
        {
            b->m_data = b->data_len;
            kroundup32(b->m_data);
            b->data = (unsigned char*)realloc(b->data, b->m_data);
        }
    if (bam_read(fp, b->data, b->data_len) != b->data_len)
        {
            return -4;
        }
    b->l_aux = b->data_len - c->n_cigar * 4 - c->l_qname - c->l_qseq - (c->l_qseq + 1) / 2;
    if (bam_is_be)
        {
            swap_endian_data(c, b->data_len, b->data);
        }
    if (bam_no_B)
        {
            bam_remove_B(b);
        }
    return 4 + block_len;
}

int
bam_validate1(const bam_header_t *header, const bam1_t *b)
{
    char *s;

    if ((b->core.tid < -1) || (b->core.mtid < -1))
        {
            return 0;
        }
    if (header && ((b->core.tid >= header->n_targets) || (b->core.mtid >= header->n_targets)))
        {
            return 0;
        }
    if (b->data_len < b->core.l_qname)
        {
            return 0;
        }
    s = (char*)memchr(bam1_qname(b), '\0', b->core.l_qname);
    if (s != &bam1_qname(b)[b->core.l_qname-1])
        {
            return 0;
        }
    return 1;
}

int
bam_remove_B(bam1_t *b)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    int end_j = 0;
    int no_qual = 0;
    unsigned int *cigar;
    unsigned int *new_cigar;
    unsigned char *seq;
    unsigned char *qual;
    unsigned char *p;

    // test if removal is necessary
    if (b->core.flag & BAM_FUNMAP)
        {
            return 0;
        }
    cigar = bam1_cigar(b);
    for (k = 0; k < b->core.n_cigar; ++k)
        {
            if (bam_cigar_op(cigar[k]) == BAM_CBACK)
                {
                    break;
                }
        }

    // no 'B'
    if (k == b->core.n_cigar)
        {
            return 0;
        }

    if (bam_cigar_op(cigar[0]) == BAM_CBACK)
        {
            goto rmB_err;
        }

    // allocate memory for the new CIGAR
    if (b->data_len + (b->core.n_cigar + 1) * 4 > b->m_data)
        {
            // not enough memory
            b->m_data = b->data_len + b->core.n_cigar * 4;
            kroundup32(b->m_data);
            b->data = (unsigned char*)realloc(b->data, b->m_data);

            // after realloc, cigar may be changed
            cigar = bam1_cigar(b);
        }

    // from the end of b->data
    new_cigar = (unsigned int*)(b->data + (b->m_data - b->core.n_cigar * 4));

    // the core loop
    seq = bam1_seq(b);
    qual = bam1_qual(b);

    // test whether base quality is available
    no_qual = (qual[0] == 0xff);
    i = j = 0;
    end_j = -1;

    for (k = l = 0; k < b->core.n_cigar; ++k)
        {
            int op  = bam_cigar_op(cigar[k]);
            int len = bam_cigar_oplen(cigar[k]);

            if (op == BAM_CBACK)
                {
                    // the backward operation
                    int t = 0;
                    int u = 0;

                    // ignore 'B' at the end of CIGAR
                    if (k == b->core.n_cigar - 1)
                        {
                            break;
                        }

                    // an excessively long backward
                    if (len > j)
                        {
                            goto rmB_err;
                        }

                    for (t = l - 1, u = 0; t >= 0; --t)
                        {
                            // look back
                            int op1  = bam_cigar_op(new_cigar[t]);
                            int len1 = bam_cigar_oplen(new_cigar[t]);
                            if (bam_cigar_type(op1)&1)
                                {
                                    // consume the query
                                    if ((u + len1) >= len)
                                        {
                                            // stop
                                            new_cigar[t] -= (len - u) << BAM_CIGAR_SHIFT;
                                            break;
                                        }
                                    else
                                        {
                                            u += len1;
                                        }
                                }
                        }

                    // squeeze out the zero-length operation
                    if (bam_cigar_oplen(new_cigar[t]) == 0)
                        {
                            --t;
                        }
                    l = t + 1;
                    end_j = j;
                    j -= len;
                }
            else
                {
                    // other CIGAR operations
                    new_cigar[l++] = cigar[k];

                    if (bam_cigar_type(op) & 1)
                        {
                            // consume the query
                            if (i != j)
                                {
                                    // no need to copy if i == j
                                    int u = 0;
                                    int c = 0;
                                    int c0 = 0;

                                    for (u = 0; u < len; ++u)
                                        {
                                            // construct the consensus
                                            c = bam1_seqi(seq, i + u);
                                            if ((j + u) < end_j)
                                                {
                                                    // in an overlap
                                                    c0 = bam1_seqi(seq, j + u);
                                                    if (c != c0)
                                                        {
                                                            // a mismatch; choose the better base
                                                            if (qual[j + u] < qual[i + u])
                                                                {
                                                                    // the base in the 2nd segment is better
                                                                    bam1_seq_seti(seq, j + u, c);
                                                                    qual[j + u] = qual[i + u] - qual[j + u];
                                                                }
                                                            else
                                                                {
                                                                    qual[j + u] -= qual[i + u]; // the 1st is better; reduce base quality
                                                                }
                                                        }
                                                    else
                                                        {
                                                            qual[j + u] = qual[j + u] > qual[i + u] ? qual[j + u] : qual[i + u];
                                                        }
                                                }
                                            else
                                                {
                                                    // not in an overlap; copy over
                                                    bam1_seq_seti(seq, j + u, c);
                                                    qual[j + u] = qual[i + u];
                                                }
                                        }
                                }
                            i += len;
                            j += len;
                        }
                }
        }

    // in very rare cases, this may be modified
    if (no_qual)
        {
            qual[0] = 0xff;
        }

    // merge adjacent operations if possible
    for (k = 1; k < l; ++k)
        {
            if (bam_cigar_op(new_cigar[k]) == bam_cigar_op(new_cigar[k-1]))
                {
                    new_cigar[k] += new_cigar[k-1] >> BAM_CIGAR_SHIFT << BAM_CIGAR_SHIFT, new_cigar[k-1] &= 0xf;
                }
        }

    // kill zero length operations
    for (k = i = 0; k < l; ++k)
        {
            if (new_cigar[k] >> BAM_CIGAR_SHIFT)
                {
                    new_cigar[i++] = new_cigar[k];
                }
        }
    l = i;

    // update b

    // set CIGAR
    memcpy(cigar, new_cigar, l * 4);
    p = b->data + b->core.l_qname + l * 4;

    // set SEQ
    memmove(p, seq, (j + 1) >> 1);
    p += (j + 1) >> 1;

    // set QUAL
    memmove(p, qual, j);
    p += j;

    // set optional fields
    memmove(p, bam1_aux(b), b->l_aux);
    p += b->l_aux;

    // update CIGAR length and query length
    b->core.n_cigar = l, b->core.l_qseq = j;

    // update record length
    b->data_len = p - b->data;

    return 0;

rmB_err:
    b->core.flag |= BAM_FUNMAP;
    return -1;
}
