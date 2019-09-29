#ifndef MACRO_H_INCLUDED
#define MACRO_H_INCLUDED

#include <math.h>

#define are_equal(v1, v2) ((fabsl((v1) - (v2))) < eps)

#define is_zero(v) (fabsl((v)) < eps)

#endif  /* MACRO_H_INCLUDED */
