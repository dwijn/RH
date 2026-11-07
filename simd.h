/* ------- file: -------------------------- simd.h ------------------

	   RH_SIMD marks a loop for vectorization: #pragma omp simd when the
	   compiler takes -fopenmp-simd (set up by CMakeLists.txt, which then
	   defines RH_OPENMP_SIMD), nothing otherwise. Use it only on loops
	   whose iterations are independent and that contain no reduction,
	   so that vectorizing them leaves every result bit for bit the same.
	   --                                              -------------- */

#ifndef SIMD_H
#define SIMD_H

#ifdef RH_OPENMP_SIMD
#define RH_SIMD _Pragma("omp simd")
#else
#define RH_SIMD
#endif

#endif /* !SIMD_H */
