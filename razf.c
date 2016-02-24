/** \file razf.c
 *  \brief Functions for handling random access compressed files
 *  \author Daniel Garrigan
 *  \version 0.5
 *  Much of the code is adapted from samtools written by Heng Li
*/

#ifndef _NO_RAZF

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "razf.h"

#if ZLIB_VERNUM < 0x1221
struct _gz_header_s
{
    int     text;
    uLong   time;
    int     xflags;
    int     os;
    Bytef   *extra;
    uInt    extra_len;
    uInt    extra_max;
    Bytef   *name;
    uInt    name_max;
    Bytef   *comment;
    uInt    comm_max;
    int     hcrc;
    int     done;
};
#warning "zlib < 1.2.2.1; RAZF writing is disabled."
#endif

#define DEF_MEM_LEVEL 8

/* gzip flag byte */
#define ASCII_FLAG   0x01 /* bit 0 set: file probably ascii text */
#define HEAD_CRC     0x02 /* bit 1 set: header CRC present */
#define EXTRA_FIELD  0x04 /* bit 2 set: extra field present */
#define ORIG_NAME    0x08 /* bit 3 set: original file name present */
#define COMMENT      0x10 /* bit 4 set: file comment present */
#define RESERVED     0xE0 /* bits 5..7: reserved */

static inline unsigned int
byte_swap_4(unsigned int v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline unsigned long long
byte_swap_8(unsigned long long v)
{
    v = ((v & 0x00000000FFFFFFFFULL) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFULL) << 16) | ((v & 0xFFFF0000FFFF0000ULL) >> 16);
    return ((v & 0x00FF00FF00FF00FFULL) << 8) | ((v & 0xFF00FF00FF00FF00ULL) >> 8);
}

static inline int
is_big_endian(void)
{
    int x = 0x01;
    char *c = (char*)&x;
    return (c[0] != 0x01);
}

#ifndef _RZ_READONLY
static void
add_zindex(RAZF *rz, long long in, long long out)
{
    if (rz->index->size == rz->index->cap)
        {
            rz->index->cap = rz->index->cap * 1.5 + 2;
            rz->index->cell_offsets = (unsigned int*)realloc(rz->index->cell_offsets,
                                      sizeof(unsigned int) * rz->index->cap);
            rz->index->bin_offsets  = (long long*)realloc(rz->index->bin_offsets,
                                      sizeof(long long) * (rz->index->cap/RZ_BIN_SIZE + 1));
        }
    if ((rz->index->size % RZ_BIN_SIZE) == 0)
        {
            rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE] = out;
        }
    rz->index->cell_offsets[rz->index->size] = out -
            rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE];
    rz->index->size ++;
}

static void
save_zindex(RAZF *rz, int fd)
{
    int i = 0;
    int v32 = 0;
    int is_be = 0;

    is_be = is_big_endian();
    if (is_be)
        {
            write(fd, &rz->index->size, sizeof(int));
        }
    else
        {
            v32 = byte_swap_4((unsigned int)rz->index->size);
            write(fd, &v32, sizeof(unsigned int));
        }
    v32 = rz->index->size / RZ_BIN_SIZE + 1;
    if (!is_be)
        {
            for (i = 0; i < v32; i++)
                {
                    rz->index->bin_offsets[i]  = byte_swap_8((unsigned long long)
                                                 rz->index->bin_offsets[i]);
                }
            for (i = 0; i < rz->index->size; i++)
                {
                    rz->index->cell_offsets[i] = byte_swap_4((unsigned int)
                                                 rz->index->cell_offsets[i]);
                }
        }
    write(fd, rz->index->bin_offsets, sizeof(long long) * v32);
    write(fd, rz->index->cell_offsets, sizeof(int) * rz->index->size);
}
#endif

static void
load_zindex(RAZF *rz, int fd)
{
    int i = 0;
    int v32 = 0;
    int is_be = 0;

    if (!rz->load_index)
        {
            return;
        }
    if (rz->index == NULL)
        {
            rz->index = (ZBlockIndex*)malloc(sizeof(ZBlockIndex));
        }
    is_be = is_big_endian();
    read(fd, &rz->index->size, sizeof(int));
    if (!is_be)
        {
            rz->index->size = byte_swap_4((unsigned int)rz->index->size);
        }
    rz->index->cap = rz->index->size;
    v32 = rz->index->size / RZ_BIN_SIZE + 1;
    rz->index->bin_offsets  = (long long*)malloc(sizeof(long long) * v32);
    read(fd, rz->index->bin_offsets, sizeof(long long) * v32);
    rz->index->cell_offsets = (unsigned int*)malloc(sizeof(unsigned int) *
                              rz->index->size);
    read(fd, rz->index->cell_offsets, sizeof(int) * rz->index->size);
    if (!is_be)
        {
            for (i = 0; i < v32; i++)
                {
                    rz->index->bin_offsets[i] = byte_swap_8((unsigned long long)
                                                            rz->index->bin_offsets[i]);
                }
            for (i = 0; i < rz->index->size; i++)
                {
                    rz->index->cell_offsets[i] = byte_swap_4((unsigned int)
                                                 rz->index->cell_offsets[i]);
                }
        }
}

#ifdef _RZ_READONLY
static RAZF *
razf_open_w(int fd)
{
    fprintf(stderr,
            "[razf_open_w] Writing is not available with zlib ver < 1.2.2.1\n");
    return 0;
}
#else
static RAZF *
razf_open_w(int fd)
{
    RAZF *rz;

    rz = (RAZF*)calloc(1, sizeof(RAZF));
    rz->mode = 'w';
    rz->filedes = fd;
    rz->stream = (z_stream*)calloc(sizeof(z_stream), 1);
    rz->inbuf  = malloc(RZ_BUFFER_SIZE);
    rz->outbuf = malloc(RZ_BUFFER_SIZE);
    rz->index = (ZBlockIndex*)calloc(sizeof(ZBlockIndex), 1);
    deflateInit2(rz->stream, RZ_COMPRESS_LEVEL, Z_DEFLATED, WINDOW_BITS + 16,
                 DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    rz->stream->avail_out = RZ_BUFFER_SIZE;
    rz->stream->next_out = (Bytef*)rz->outbuf;
    rz->header = (gz_header*)calloc(sizeof(gz_header), 1);
    rz->header->os = 0x03;
    rz->header->text = 0;
    rz->header->time = 0;
    rz->header->extra = (Bytef*)malloc(7);
    strncpy((char*)rz->header->extra, "RAZF", 4);
    rz->header->extra[4] = 1;
    rz->header->extra[5] = RZ_BLOCK_SIZE >> 8;
    rz->header->extra[6] = RZ_BLOCK_SIZE & 0xFF;
    rz->header->extra_len = 7;
    rz->header->name = rz->header->comment  = 0;
    rz->header->hcrc = 0;
    deflateSetHeader(rz->stream, rz->header);
    rz->block_pos = rz->block_off = 0;
    return rz;
}

static void
_razf_write(RAZF* rz, const void *data, int size)
{
    int tout = 0;

    rz->stream->avail_in = size;
    rz->stream->next_in  = (Bytef*)data;
    while (1)
        {
            tout = rz->stream->avail_out;
            deflate(rz->stream, Z_NO_FLUSH);
            rz->out += tout - rz->stream->avail_out;
            if (rz->stream->avail_out)
                {
                    break;
                }
            write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
            rz->stream->avail_out = RZ_BUFFER_SIZE;
            rz->stream->next_out  = (Bytef*)rz->outbuf;
            if (rz->stream->avail_in == 0)
                {
                    break;
                }
        };
    rz->in += size - rz->stream->avail_in;
    rz->block_off += size - rz->stream->avail_in;
}

static void
razf_flush(RAZF *rz)
{
    unsigned int tout = 0;

    if (rz->buf_len)
        {
            _razf_write(rz, rz->inbuf, rz->buf_len);
            rz->buf_off = rz->buf_len = 0;
        }
    if (rz->stream->avail_out)
        {
            write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
            rz->stream->avail_out = RZ_BUFFER_SIZE;
            rz->stream->next_out  = (Bytef*)rz->outbuf;
        }
    while (1)
        {
            tout = rz->stream->avail_out;
            deflate(rz->stream, Z_FULL_FLUSH);
            rz->out += tout - rz->stream->avail_out;
            if (rz->stream->avail_out == 0)
                {
                    write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
                    rz->stream->avail_out = RZ_BUFFER_SIZE;
                    rz->stream->next_out  = (Bytef*)rz->outbuf;
                }
            else
                {
                    break;
                }
        }
    rz->block_pos = rz->out;
    rz->block_off = 0;
}

static void
razf_end_flush(RAZF *rz)
{
    unsigned int tout = 0;

    if (rz->buf_len)
        {
            _razf_write(rz, rz->inbuf, rz->buf_len);
            rz->buf_off = rz->buf_len = 0;
        }
    while (1)
        {
            tout = rz->stream->avail_out;
            deflate(rz->stream, Z_FINISH);
            rz->out += tout - rz->stream->avail_out;
            if(rz->stream->avail_out < RZ_BUFFER_SIZE)
                {
                    write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
                    rz->stream->avail_out = RZ_BUFFER_SIZE;
                    rz->stream->next_out  = (Bytef*)rz->outbuf;
                }
            else
                {
                    break;
                }
        }
}

static void
_razf_buffered_write(RAZF *rz, const void *data, int size)
{
    int i = 0;
    int n = 0;
    char *dd = (char*)data;

    while (1)
        {
            if (rz->buf_len == RZ_BUFFER_SIZE)
                {
                    _razf_write(rz, rz->inbuf, rz->buf_len);
                    rz->buf_len = 0;
                }
            if (size + rz->buf_len < RZ_BUFFER_SIZE)
                {
                    for (i = 0; i < size; i++)
                        {
                            ((char*)rz->inbuf + rz->buf_len)[i] = dd[i];
                        }
                    rz->buf_len += size;
                    return;
                }
            else
                {
                    n = RZ_BUFFER_SIZE - rz->buf_len;
                    for (i = 0; i < n; i++)
                        {
                            ((char*)rz->inbuf + rz->buf_len)[i] = dd[i];
                        }
                    size -= n;
                    dd += n;
                    rz->buf_len += n;
                }
        }
}

int
razf_write(RAZF* rz, const void *data, int size)
{
    int ori_size = 0;
    int n = 0;
    long long next_block = 0;
    char *dd = (char*)data;

    ori_size = size;
    next_block = ((rz->in / RZ_BLOCK_SIZE) + 1) * RZ_BLOCK_SIZE;
    while ((rz->in + rz->buf_len+size) >= next_block)
        {
            n = next_block - rz->in - rz->buf_len;
            _razf_buffered_write(rz, (void*)dd, n);
            dd += n;
            size -= n;
            razf_flush(rz);
            add_zindex(rz, rz->in, rz->out);
            next_block = ((rz->in/RZ_BLOCK_SIZE)+1)*RZ_BLOCK_SIZE;
        }
    _razf_buffered_write(rz, (void*)dd, size);
    return ori_size;
}
#endif

static int
_read_gz_header(unsigned char *data, int size, int *extra_off, int *extra_len)
{
    int method = 0;
    int flags = 0;
    int n = 0;
    int len = 0;

    if (size < 2)
        {
            return 0;
        }
    if ((data[0] != 0x1f) || (data[1] != 0x8b))
        {
            return 0;
        }
    if (size < 4)
        {
            return 0;
        }
    method = data[2];
    flags  = data[3];
    if ((method != Z_DEFLATED) || (flags & RESERVED))
        {
            return 0;
        }

    //skip 6 bytes
    n = 4 + 6;
    *extra_off = n + 2;
    *extra_len = 0;
    if (flags & EXTRA_FIELD)
        {
            if (size < n + 2)
                {
                    return 0;
                }
            len = ((int)data[n + 1] << 8) | data[n];
            n += 2;
            *extra_off = n;
            while (len)
                {
                    if (n >= size)
                        {
                            return 0;
                        }
                    n++;
                    len--;
                }
            *extra_len = n - (*extra_off);
        }
    if (flags & ORIG_NAME)
        while((n < size) && data[n++]);
    if (flags & COMMENT)
        while((n < size) && data[n++]);
    if (flags & HEAD_CRC)
        {
            if (n + 2 > size)
                {
                    return 0;
                }
            n += 2;
        }
    return n;
}

static RAZF *
razf_open_r(int fd, int _load_index)
{
    int ext_off = 0;
    int ext_len = 0;
    int n = 0;
    int is_be = 0;
    int ret = 0;
    long long end = 0;
    unsigned char c[] = "RAZF";
    RAZF *rz;

    rz = (RAZF*)calloc(1, sizeof(RAZF));
    rz->mode = 'r';
    rz->filedes = fd;
    rz->stream = (z_stream*)calloc(sizeof(z_stream), 1);
    rz->inbuf  = malloc(RZ_BUFFER_SIZE);
    rz->outbuf = malloc(RZ_BUFFER_SIZE);
    rz->end = rz->src_end = 0x7FFFFFFFFFFFFFFFLL;
    n = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
    ret = _read_gz_header((unsigned char*)rz->inbuf, n, &ext_off, &ext_len);
    if (ret == 0)
        {
PLAIN_FILE:
            rz->in = n;
            rz->file_type = FILE_TYPE_PLAIN;
            memcpy(rz->outbuf, rz->inbuf, n);
            rz->buf_len = n;
            free(rz->stream);
            rz->stream = NULL;
            return rz;
        }
    rz->header_size = ret;
    ret = inflateInit2(rz->stream, -WINDOW_BITS);
    if (ret != Z_OK)
        {
            inflateEnd(rz->stream);
            goto PLAIN_FILE;
        }
    rz->stream->avail_in = n - rz->header_size;
    rz->stream->next_in  = (Bytef*)rz->inbuf + rz->header_size;
    rz->stream->avail_out = RZ_BUFFER_SIZE;
    rz->stream->next_out  = (Bytef*)rz->outbuf;
    rz->file_type = FILE_TYPE_GZ;
    rz->in = rz->header_size;
    rz->block_pos = rz->header_size;
    rz->next_block_pos = rz->header_size;
    rz->block_off = 0;
    if ((ext_len < 7) || (memcmp((Bytef*)rz->inbuf + ext_off, c, 4) != 0))
        {
            return rz;
        }
    if (((((unsigned char*)rz->inbuf)[ext_off + 5] << 8) | ((
                unsigned char*)rz->inbuf)[ext_off + 6]) != RZ_BLOCK_SIZE)
        {
            fprintf(stderr,
                    " -- WARNING: RZ_BLOCK_SIZE is not %d, treat source as gz file.  in %s -- %s:%d --\n",
                    RZ_BLOCK_SIZE, __FUNCTION__, __FILE__, __LINE__);
            return rz;
        }
    rz->load_index = _load_index;
    rz->file_type = FILE_TYPE_RZ;
    if(lseek(fd, -16, SEEK_END) == -1)
        {
UNSEEKABLE:
            rz->seekable = 0;
            rz->index = NULL;
            rz->src_end = rz->end = 0x7FFFFFFFFFFFFFFFLL;
        }
    else
        {
            is_be = is_big_endian();
            rz->seekable = 1;
            read(fd, &end, sizeof(long long));
            if (!is_be)
                {
                    rz->src_end = (long long)byte_swap_8((unsigned long long)end);
                }
            else
                {
                    rz->src_end = end;
                }

            read(fd, &end, sizeof(long long));
            if (!is_be)
                {
                    rz->end = (long long)byte_swap_8((unsigned long long)end);
                }
            else
                {
                    rz->end = end;
                }
            if (n > rz->end)
                {
                    rz->stream->avail_in -= n - rz->end;
                    n = rz->end;
                }
            if (rz->end > rz->src_end)
                {
                    lseek(fd, rz->in, SEEK_SET);
                    goto UNSEEKABLE;
                }
            if (lseek(fd, rz->end, SEEK_SET) != rz->end)
                {
                    lseek(fd, rz->in, SEEK_SET);
                    goto UNSEEKABLE;
                }
            load_zindex(rz, fd);
            lseek(fd, n, SEEK_SET);
        }
    return rz;
}

RAZF *
razf_dopen(int fd, const char *mode)
{
    if (strstr(mode, "r"))
        {
            return razf_open_r(fd, 1);
        }
    else if (strstr(mode, "w"))
        {
            return razf_open_w(fd);
        }
    else
        {
            return NULL;
        }
}

RAZF *
razf_dopen2(int fd, const char *mode)
{
    if (strstr(mode, "r"))
        {
            return razf_open_r(fd, 0);
        }
    else if (strstr(mode, "w"))
        {
            return razf_open_w(fd);
        }
    else
        {
            return NULL;
        }
}

static inline RAZF *
_razf_open(const char *filename, const char *mode, int _load_index)
{
    int fd = 0;
    RAZF *rz;

    if (strstr(mode, "r"))
        {
            fd = open(filename, O_RDONLY);
            if (fd < 0)
                {
                    return NULL;
                }
            rz = razf_open_r(fd, _load_index);
        }
    else if (strstr(mode, "w"))
        {
            fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0666);
            if (fd < 0)
                {
                    return NULL;
                }
            rz = razf_open_w(fd);
        }
    else
        {
            return NULL;
        }
    return rz;
}

RAZF *
razf_open(const char *filename, const char *mode)
{
    return _razf_open(filename, mode, 1);
}

RAZF *
razf_open2(const char *filename, const char *mode)
{
    return _razf_open(filename, mode, 0);
}

int
razf_get_data_size(RAZF *rz, long long *u_size, long long *c_size)
{
    long long n = 0;

    if ((rz->mode != 'r') && (rz->mode != 'R'))
        {
            return 0;
        }
    switch (rz->file_type)
        {
        case FILE_TYPE_PLAIN:
            if (rz->end == 0x7fffffffffffffffLL)
                {
                    if ((n = lseek(rz->filedes, 0, SEEK_CUR)) == -1)
                        {
                            return 0;
                        }
                    rz->end = lseek(rz->filedes, 0, SEEK_END);
                    lseek(rz->filedes, n, SEEK_SET);
                }
            *u_size = *c_size = rz->end;
            return 1;
        case FILE_TYPE_GZ :
            return 0;
        case FILE_TYPE_RZ :
            if (rz->src_end == rz->end)
                {
                    return 0;
                }
            *u_size = rz->src_end;
            *c_size = rz->end;
            return 1;
        default:
            return 0;
        }
}

static int
_razf_read(RAZF* rz, void *data, int size)
{
    int ret = 0;
    int tin = 0;

    if (rz->z_eof || rz->z_err)
        {
            return 0;
        }
    if (rz->file_type == FILE_TYPE_PLAIN)
        {
            ret = read(rz->filedes, data, size);
            if (ret == 0)
                {
                    rz->z_eof = 1;
                }
            return ret;
        }
    rz->stream->avail_out = size;
    rz->stream->next_out  = (Bytef*)data;
    while (rz->stream->avail_out)
        {
            if (rz->stream->avail_in == 0)
                {
                    if (rz->in >= rz->end)
                        {
                            rz->z_eof = 1;
                            break;
                        }
                    if (rz->end - rz->in < RZ_BUFFER_SIZE)
                        {
                            rz->stream->avail_in = read(rz->filedes, rz->inbuf, rz->end -rz->in);
                        }
                    else
                        {
                            rz->stream->avail_in = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
                        }

                    if (rz->stream->avail_in == 0)
                        {
                            rz->z_eof = 1;
                            break;
                        }
                    rz->stream->next_in = (Bytef*)rz->inbuf;
                }
            tin = rz->stream->avail_in;
            ret = inflate(rz->stream, Z_BLOCK);
            rz->in += tin - rz->stream->avail_in;
            if ((ret == Z_NEED_DICT) || (ret == Z_MEM_ERROR) || (ret == Z_DATA_ERROR))
                {
                    fprintf(stderr, "[_razf_read] inflate error: %d %s (at %s:%d)\n", ret,
                            rz->stream->msg ? rz->stream->msg : "", __FILE__, __LINE__);
                    rz->z_err = 1;
                    break;
                }
            if (ret == Z_STREAM_END)
                {
                    rz->z_eof = 1;
                    break;
                }
            if ((rz->stream->data_type & 128) && !(rz->stream->data_type & 64))
                {
                    rz->buf_flush = 1;
                    rz->next_block_pos = rz->in;
                    break;
                }
        }
    return size - rz->stream->avail_out;
}

int
razf_read(RAZF *rz, void *data, int size)
{
    int ori_size = 0;
    int i = 0;
    char *dd = (char*)data;

    ori_size = size;

    while (size > 0)
        {
            if (rz->buf_len)
                {
                    if (size < rz->buf_len)
                        {
                            for (i = 0; i < size; i++)
                                {
                                    dd[i] = ((char*)rz->outbuf + rz->buf_off)[i];
                                }
                            rz->buf_off += size;
                            rz->buf_len -= size;
                            dd += size;
                            rz->block_off += size;
                            size = 0;
                            break;
                        }
                    else
                        {
                            for (i = 0; i < rz->buf_len; i++)
                                {
                                    dd[i] = ((char*)rz->outbuf + rz->buf_off)[i];
                                }
                            dd += rz->buf_len;
                            size -= rz->buf_len;
                            rz->block_off += rz->buf_len;
                            rz->buf_off = 0;
                            rz->buf_len = 0;
                            if (rz->buf_flush)
                                {
                                    rz->block_pos = rz->next_block_pos;
                                    rz->block_off = 0;
                                    rz->buf_flush = 0;
                                }
                        }
                }
            else if (rz->buf_flush)
                {
                    rz->block_pos = rz->next_block_pos;
                    rz->block_off = 0;
                    rz->buf_flush = 0;
                }
            if (rz->buf_flush)
                {
                    continue;
                }
            rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
            if (rz->z_eof && (rz->buf_len == 0))
                {
                    break;
                }
        }
    rz->out += ori_size - size;
    return ori_size - size;
}

int
razf_skip(RAZF* rz, int size)
{
    int ori_size = 0;

    ori_size = size;
    while (size > 0)
        {
            if (rz->buf_len)
                {
                    if (size < rz->buf_len)
                        {
                            rz->buf_off += size;
                            rz->buf_len -= size;
                            rz->block_off += size;
                            size = 0;
                            break;
                        }
                    else
                        {
                            size -= rz->buf_len;
                            rz->buf_off = 0;
                            rz->buf_len = 0;
                            rz->block_off += rz->buf_len;

                            if (rz->buf_flush)
                                {
                                    rz->block_pos = rz->next_block_pos;
                                    rz->block_off = 0;
                                    rz->buf_flush = 0;
                                }
                        }
                }
            else if (rz->buf_flush)
                {
                    rz->block_pos = rz->next_block_pos;
                    rz->block_off = 0;
                    rz->buf_flush = 0;
                }
            if (rz->buf_flush)
                {
                    continue;
                }
            rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
            if (rz->z_eof || rz->z_err)
                {
                    break;
                }
        }
    rz->out += ori_size - size;
    return ori_size - size;
}

static void
_razf_reset_read(RAZF *rz, long long in, long long out)
{
    lseek(rz->filedes, in, SEEK_SET);
    rz->in  = in;
    rz->out = out;
    rz->block_pos = in;
    rz->next_block_pos = in;
    rz->block_off = 0;
    rz->buf_flush = 0;
    rz->z_eof = rz->z_err = 0;
    inflateReset(rz->stream);
    rz->stream->avail_in = 0;
    rz->buf_off = rz->buf_len = 0;
}

long long
razf_jump(RAZF *rz, long long block_start, int block_offset)
{
    long long pos = 0;

    rz->z_eof = 0;
    if (rz->file_type == FILE_TYPE_PLAIN)
        {
            rz->buf_off = rz->buf_len = 0;
            pos = block_start + block_offset;
            pos = lseek(rz->filedes, pos, SEEK_SET);
            rz->out = rz->in = pos;
            return pos;
        }
    if ((block_start == rz->block_pos) && (block_offset >= rz->block_off))
        {
            block_offset -= rz->block_off;
            // needn't reset inflate
            goto SKIP;
        }

    // automaticly revist wrong block_start
    if (block_start == 0)
        {
            block_start = rz->header_size;
        }
    _razf_reset_read(rz, block_start, 0);
SKIP:
    if (block_offset)
        {
            razf_skip(rz, block_offset);
        }
    return rz->block_off;
}

long long
razf_seek(RAZF* rz, long long pos, int wher)
{
    long long idx = 0;
    long long seek_pos = 0;
    long long new_out = 0;

    rz->z_eof = 0;
    if (wher == SEEK_CUR)
        {
            pos += rz->out;
        }
    else if (wher == SEEK_END)
        {
            pos += rz->src_end;
        }
    if (rz->file_type == FILE_TYPE_PLAIN)
        {
            seek_pos = lseek(rz->filedes, pos, SEEK_SET);
            rz->buf_off = rz->buf_len = 0;
            rz->out = rz->in = seek_pos;
            return seek_pos;
        }
    else if (rz->file_type == FILE_TYPE_GZ)
        {
            if (pos >= rz->out)
                {
                    goto SKIP;
                }
            return rz->out;
        }
    if (pos == rz->out)
        {
            return pos;
        }
    if (pos > rz->src_end)
        {
            return rz->out;
        }
    if (!rz->seekable || !rz->load_index)
        if(pos >= rz->out)
            {
                goto SKIP;
            }
    idx = pos / RZ_BLOCK_SIZE - 1;
    seek_pos = (idx < 0) ? rz->header_size:(rz->index->cell_offsets[idx] +
                                            rz->index->bin_offsets[idx / RZ_BIN_SIZE]);
    new_out  = (idx + 1) * RZ_BLOCK_SIZE;
    if ((pos > rz->out) && (new_out <= rz->out))
        {
            goto SKIP;
        }
    _razf_reset_read(rz, seek_pos, new_out);
SKIP:
    razf_skip(rz, (int)(pos - rz->out));
    return rz->out;
}

unsigned long long
razf_tell2(RAZF *rz)
{
    return (unsigned long long)rz->block_pos << 16 | (rz->block_off & 0xffff);
}

long long
razf_seek2(RAZF *rz, unsigned long long voffset, int wher)
{
    if (wher != SEEK_SET)
        {
            return -1;
        }
    return razf_jump(rz, voffset >> 16, voffset & 0xffff);
}

void
razf_close(RAZF *rz)
{
    if (rz->mode == 'w')
        {
#ifndef _RZ_READONLY
            razf_end_flush(rz);
            deflateEnd(rz->stream);
            save_zindex(rz, rz->filedes);
            if (is_big_endian())
                {
                    write(rz->filedes, &rz->in, sizeof(long long));
                    write(rz->filedes, &rz->out, sizeof(long long));
                }
            else
                {
                    unsigned long long v64 = byte_swap_8((unsigned long long)rz->in);
                    write(rz->filedes, &v64, sizeof(long long));
                    v64 = byte_swap_8((unsigned long long)rz->out);
                    write(rz->filedes, &v64, sizeof(long long));
                }
#endif
        }
    else if (rz->mode == 'r')
        {
            if(rz->stream)
                {
                    inflateEnd(rz->stream);
                }
        }
    if (rz->inbuf)
        {
            free(rz->inbuf);
        }
    if (rz->outbuf)
        {
            free(rz->outbuf);
        }
    if (rz->header)
        {
            free(rz->header->extra);
            free(rz->header->name);
            free(rz->header->comment);
            free(rz->header);
        }
    if (rz->index)
        {
            free(rz->index->bin_offsets);
            free(rz->index->cell_offsets);
            free(rz->index);
        }
    free(rz->stream);
    close(rz->filedes);
    free(rz);
}

#endif
