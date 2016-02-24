#ifndef BAM_ENDIAN_H
#define BAM_ENDIAN_H

static inline int
bam_is_big_endian()
{
    long int one = 0;
    one = 1;
    return !(*((char*)(&one)));
}

static inline unsigned short
bam_swap_endian_2(unsigned short v)
{
    return (unsigned short)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}

static inline void *
bam_swap_endian_2p(void *x)
{
    *(unsigned short*)x = bam_swap_endian_2(*(unsigned short*)x);
    return x;
}

static inline unsigned int
bam_swap_endian_4(unsigned int v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline void *
bam_swap_endian_4p(void *x)
{
    *(unsigned int*)x = bam_swap_endian_4(*(unsigned int*)x);
    return x;
}

static inline unsigned long long
bam_swap_endian_8(unsigned long long v)
{
    v = ((v & 0x00000000FFFFFFFFULL) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFULL) << 16) | ((v & 0xFFFF0000FFFF0000ULL) >> 16);
    return ((v & 0x00FF00FF00FF00FFULL) << 8) | ((v & 0xFF00FF00FF00FF00ULL) >> 8);
}

static inline void *
bam_swap_endian_8p(void *x)
{
    *(unsigned long long*)x = bam_swap_endian_8(*(unsigned long long*)x);
    return x;
}

#endif
