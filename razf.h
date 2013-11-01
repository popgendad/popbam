 /*-
 * RAZF : Random Access compressed(Z) File
 * Version: 1.0
 * Release Date: 2008-10-27
 *
 * Copyright 2008, Jue Ruan <ruanjue@gmail.com>, Heng Li <lh3@sanger.ac.uk>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef __RAZF_RJ_H
#define __RAZF_RJ_H

#include <stdio.h>
#include <zlib.h>
#include <zconf.h>

#if ZLIB_VERNUM < 0x1221
#define _RZ_READONLY
struct _gz_header_s;
typedef struct _gz_header_s _gz_header;
#define gz_header _gz_header
#endif

#define WINDOW_BITS 15

#ifndef RZ_BLOCK_SIZE
#define RZ_BLOCK_SIZE (1 << WINDOW_BITS)
#endif

#ifndef RZ_BUFFER_SIZE
#define RZ_BUFFER_SIZE 4096
#endif

#ifndef RZ_COMPRESS_LEVEL
#define RZ_COMPRESS_LEVEL 6
#endif

#define RZ_BIN_SIZE ((1ULL << 32) / RZ_BLOCK_SIZE)

typedef struct
{
	unsigned int *cell_offsets;              // i
	long long *bin_offsets;                  // i / BIN_SIZE
	int size;
	int cap;
} ZBlockIndex;

// When storing index, output bytes in Big-Endian everywhere

#define FILE_TYPE_RZ	1
#define FILE_TYPE_PLAIN	2
#define FILE_TYPE_GZ	3

typedef struct RandomAccessZFile
{
	char mode;                    // 'w' : write mode; 'r' : read mode
	int file_type;                // plain file or rz file, razf_read support plain file as input too, in this case, razf_read work as buffered fread
	int filedes;                  // the file descriptor
	z_stream *stream;
	ZBlockIndex *index;
	long long in;                   // in: n bytes total in
	long long out;                  // out: n bytes total out
	long long end;                  // end: the end of all data blocks, while the start of index
	long long src_end;              // src_end: the true end position in uncompressed file
	int buf_flush;                // buffer should be flush, suspend inflate util buffer is empty
	long long block_pos;            // block_pos: the start postiion of current block in compressed file
	long long block_off;            // block_off: tell how many bytes have been read from current block
	long long next_block_pos;
	void *inbuf;
	void *outbuf;
	int header_size;
	gz_header *header;            // header is used to transfer inflate_state->mode from HEAD to TYPE after call inflateReset
	int buf_off;
	int buf_len;
	int z_err;
	int z_eof;
	int seekable;                 // Indice where the source is seekable
	int load_index;               // set has_index to 0 in mode 'w', then index will be discarded
} RAZF;

#ifdef __cplusplus
extern "C" {
#endif

RAZF* razf_dopen(int data_fd, const char *mode);
RAZF *razf_open(const char *fn, const char *mode);
int razf_write(RAZF* rz, const void *data, int size);
int razf_read(RAZF* rz, void *data, int size);
long long razf_seek(RAZF* rz, long long pos, int where);
void razf_close(RAZF* rz);

#define razf_tell(rz) ((rz)->out)

RAZF* razf_open2(const char *filename, const char *mode);
RAZF* razf_dopen2(int fd, const char *mode);
unsigned long long razf_tell2(RAZF *rz);
long long razf_seek2(RAZF *rz, unsigned long long voffset, int where);

#ifdef __cplusplus
}
#endif

#endif
