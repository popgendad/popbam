/** \file pop_error.hpp
 *  \brief Header file for the pop_error.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
 */

#ifndef POP_ERROR_HPP
#define POP_ERROR_HPP

/*!
 * \struct errmod_coef_t
 * \brief A structure to hold the coefficients necessary in the error model
 */
typedef struct __errmod_coef_t
{
    double *fk;                       //!< Pointer to
    double *beta;                     //!< Pointer to
    double *lhet;                     //!< Pointer to
} errmod_coef_t;

/*!
 * \struct errmod_t
 * \brief A structure to hold data for the error model
 */
typedef struct __errmod_t
{
    double depcorr;                   //!< Dependency correlation
    errmod_coef_t *coef;              //!< Pre-computed coefficients
} errmod_t;

/*!
 * \struct call_aux_t
 * \brief A structure to hold auxiliary information for use in error model
 */
typedef struct __call_aux_t
{
    double fsum[16];                  //!< Array of
    double bsum[16];                  //!< Array of
    unsigned int c[16];               //!< Array of
} call_aux_t;

/*!
 * \fn static errmod_coef_t* cal_coef(double depcorr, double eta)
 * \brief Calculate coefficients in error model
 * \param depcorr The constant for the dependency correlation
 * \param eta
 */

extern errmod_coef_t* cal_coef (double depcorr, double eta);

/*!
 * \fn errmod_t *errmod_init(float depcorr)
 * \brief Initialize the error model data structure
 * \param depcorr The constant for the dependency correlation
 */

extern errmod_t* errmod_init (float depcorr);

/*!
 * \fn void errmod_destroy(errmod_t *em)
 * \brief Deallocate memory for error model data structure
 * \param em Pointer to error model data structure
 */

extern void errmod_destroy (errmod_t *em);

/*!
 * \fn int errmod_cal(const errmod_t *em, uint16_t n, int m, uint16_t *bases, float *q)
 * \brief Calculates probability for error model
 * \param em Pointer to the error model data structure
 * \param n The number of bases
 * \param m The maximum base
 * \param bases[i] qual:6, strand:1, base:4
 * \param q[i*m+j] Phred-scaled likelihood of (i,j)
 */

extern int errmod_cal (const errmod_t *em, uint16_t n, int m, uint16_t *bases,
                       float *q);

#endif
