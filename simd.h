/* ------- file: -------------------------- simd.h ------------------

	   RH_SIMD marks a loop for vectorization: #pragma omp simd when the
	   compiler takes -fopenmp-simd (set up by CMakeLists.txt, which then
	   defines RH_OPENMP_SIMD), nothing otherwise. Use it only on loops
	   whose iterations are independent and that contain no reduction,
	   and declare every scratch variable inside the loop body (the
	   pragma makes variables declared outside it shared between lanes).
	   Vectorizing such a loop reorders no sum; results then change only
	   where the compiler fuses multiply-adds differently (FMA).
	   --                                              -------------- */

#ifndef SIMD_H
#define SIMD_H

#ifdef RH_OPENMP_SIMD
#define RH_SIMD _Pragma("omp simd")
#else
#define RH_SIMD
#endif

#endif /* !SIMD_H */
