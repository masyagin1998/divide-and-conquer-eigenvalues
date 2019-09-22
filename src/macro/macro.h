#ifndef MACRO_H_INCLUDED
#define MACRO_H_INCLUDED

#define are_equal(v1, v2) ((fabs((v1) - (v2))) < eps)

#define is_zero(v) (fabs((v)) < eps)

#endif  /* MACRO_H_INCLUDED */
