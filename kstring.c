/** \file kstring.c
 *  \brief Functions for string handling
 *  \author Daniel Garrigan
 *  \version 0.5
 *  Much of the code is adapted from samtools written by Heng Li
*/

#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "kstring.h"

int
ksprintf(kstring_t *s, const char *fmt, ...)
{
    va_list ap;
    int l;

    va_start(ap, fmt);

    // This line does not work with glibc 2.0. See `man snprintf'.
    l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap);
    va_end(ap);
    if (l + 1 > s->m - s->l)
        {
            s->m = s->l + l + 2;
            kroundup32(s->m);
            s->s = (char*)realloc(s->s, s->m);
            va_start(ap, fmt);
            l = vsnprintf(s->s + s->l, s->m - s->l, fmt, ap);
        }
    va_end(ap);
    s->l += l;
    return l;
}

char *
kstrtok(const char *str, const char *sep, ks_tokaux_t *aux)
{
    const char *p;
    const char *start;

    // set up the table
    if (sep)
        {
            // no need to set up if we have finished
            if ((str == 0) && (aux->tab[0] & 1))
                {
                    return 0;
                }
            aux->finished = 0;
            if (sep[1])
                {
                    aux->sep = -1;
                    aux->tab[0] = aux->tab[1] = aux->tab[2] = aux->tab[3] = 0;
                    for (p = sep; *p; ++p)
                        {
                            aux->tab[*p >> 6] |= 1ull << (*p & 0x3f);
                        }
                }
            else
                aux->sep = sep[0];
        }
    if (aux->finished)
        {
            return 0;
        }
    else if (str)
        {
            aux->p = str - 1;
            aux->finished = 0;
        }
    if (aux->sep < 0)
        {
            for (p = start = aux->p + 1; *p; ++p)
                {
                    if (aux->tab[*p >> 6] >> (*p & 0x3f) & 1)
                        {
                            break;
                        }
                }
        }
    else
        {
            for (p = start = aux->p + 1; *p; ++p)
                {
                    if (*p == aux->sep)
                        {
                            break;
                        }
                }
        }

    // end of token
    aux->p = p;

    // no more tokens
    if (*p == 0)
        {
            aux->finished = 1;
        }
    return (char*)start;
}

// s MUST BE a null terminated string; l = strlen(s)
int
ksplit_core(char *s, int delimiter, int *_max, int **_offsets)
{
    int i;
    int n;
    int max;
    int last_char;
    int last_start;
    int *offsets;
    int l;

    n = 0;
    max = *_max;
    offsets = *_offsets;
    l = strlen(s);

#define __ksplit_aux do {												\
		if (_offsets) {													\
			s[i] = 0;													\
			if (n == max) {												\
				max = max ? max << 1 : 2;									\
				offsets = (int*)realloc(offsets, sizeof(int) * max);	\
			}															\
			offsets[n++] = last_start;									\
		} else ++n;														\
	} while (0)

    for (i = 0, last_char = last_start = 0; i <= l; ++i)
        {
            if (delimiter == 0)
                {
                    if (isspace(s[i]) || (s[i] == 0))
                        {
                            // the end of a field
                            if (isgraph(last_char))
                                {
                                    __ksplit_aux;
                                }
                        }
                    else
                        {
                            if (isspace(last_char) || (last_char == 0))
                                {
                                    last_start = i;
                                }
                        }
                }
            else
                {
                    if ((s[i] == delimiter) || (s[i] == 0))
                        {
                            // the end of a field
                            if ((last_char != 0) && (last_char != delimiter))
                                {
                                    __ksplit_aux;
                                }
                        }
                    else
                        {
                            if ((last_char == delimiter) || (last_char == 0))
                                {
                                    last_start = i;
                                }
                        }
                }
            last_char = s[i];
        }
    *_max = max;
    *_offsets = offsets;
    return n;
}

/**********************
 * Boyer-Moore search *
 **********************/

// reference: http://www-igm.univ-mlv.fr/~lecroq/string/node14.html
static int *
ksBM_prep(const unsigned char *pat, int m)
{
    int i;
    int *suff;
    int *prep;
    int *bmGs;
    int *bmBc;

    prep = (int*)calloc((size_t)(m + 256), sizeof(int));
    bmGs = prep;
    bmBc = prep + m;
    {
        // preBmBc()
        for (i = 0; i < 256; ++i)
            {
                bmBc[i] = m;
            }
        for (i=0; i < m-1; ++i)
            {
                bmBc[pat[i]] = m - i - 1;
            }
    }
    suff = (int*)calloc((size_t)m, sizeof(int));
    {
        // suffixes()
        int f = 0;
        int g;
        suff[m - 1] = m;
        g = m - 1;

        for (i = m - 2; i >= 0; --i)
            {
                if ((i > g) && (suff[i + m - 1 - f] < i - g))
                    {
                        suff[i] = suff[i + m - 1 - f];
                    }
                else
                    {
                        if (i < g)
                            {
                                g = i;
                            }
                        f = i;
                        while ((g >= 0) && (pat[g] == pat[g + m - 1 - f]))
                            {
                                --g;
                            }
                        suff[i] = f - g;
                    }
            }
    }
    {
        // preBmGs()
        int j = 0;

        for (i = 0; i < m; ++i)
            {
                bmGs[i] = m;
            }
        for (i = m - 1; i >= 0; --i)
            {
                if (suff[i] == i + 1)
                    {
                        for (; j < m-1-i; ++j)
                            {
                                if (bmGs[j] == m)
                                    {
                                        bmGs[j] = m - 1 - i;
                                    }
                            }
                    }
            }
        for (i = 0; i <= m - 2; ++i)
            {
                bmGs[m - 1 - suff[i]] = m - 1 - i;
            }
    }
    free(suff);
    return prep;
}

void *
kmemmem(const void *_str, int n, const void *_pat, int m, int **_prep)
{
    int i = 0;
    int j = 0;
    int *prep = 0;
    int *bmGs;
    int *bmBc;
    const unsigned char *str;
    const unsigned char *pat;

    str = (const unsigned char*)_str;
    pat = (const unsigned char*)_pat;
    prep = ((_prep == 0) || (*_prep == 0)) ? ksBM_prep(pat, m) : *_prep;
    if (_prep && (*_prep == 0))
        {
            *_prep = prep;
        }
    bmGs = prep;
    bmBc = prep + m;
    j = 0;
    while (j <= (n - m))
        {
            for (i = m - 1; (i >= 0) && (pat[i] == str[i + j]); --i);

            if (i >= 0)
                {
                    int max = bmBc[str[i + j]] - m + 1 + i;
                    if (max < bmGs[i])
                        {
                            max = bmGs[i];
                        }
                    j += max;
                }
            else
                {
                    return (void*)(str + j);
                }
        }
    if (_prep == 0)
        {
            free(prep);
        }
    return 0;
}

char *
kstrstr(const char *str, const char *pat, int **_prep)
{
    return (char*)kmemmem(str, strlen(str), pat, strlen(pat), _prep);
}

char *kstrnstr(const char *str, const char *pat, int n, int **_prep)
{
    return (char*)kmemmem(str, n, pat, strlen(pat), _prep);
}
