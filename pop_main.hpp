/** \file pop_main.hpp
 *  \brief Header for the pop_main.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5b
 */

#ifndef POP_MAIN_H
#define POP_MAIN_H

#include <cstdint>
#include "pop_sample.hpp"
#include "pop_utils.hpp"
#include "pop_global.hpp"

///
/// Define classes
///

/*!
 * \class popbamData
 * \brief The abstract base class for passing parameters and data
 */

class popbamData
{
public:
    // default constructor
    popbamData();

    // destructor
    ~popbamData() {}

    // member functions
    int assignPops(const popbamOptions *p);

    // member variables
    int tid;                       //!< Reference chromosome/scaffold identifier
    int beg;                       //!< Reference coordinate of the beginning of the current region
    int end;                       //!< Reference coordinate of the end of current region
    int len;                       //!< Length of the reference sequence for current region
    int num_sites;                 //!< Total number of aligned sites
    int segsites;                  //!< Total number of segregating sites in entire sample
    char *ref_base;                //!< Reference sequence string for specified region
    uint8_t *pop_nsmpl;            //!< Sample size per population
    uint64_t *types;               //!< The site type for each aligned site
    uint64_t *pop_mask;            //!< Bit mask for which individuals are in which population
    popbam_func_t derived_type;    //!< Type of the derived class
};

const char version[] = "0.5b";

int popbam_usage(void);

#endif
