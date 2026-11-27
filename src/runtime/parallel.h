/* ------- file: -------------------------- parallel.h --------------
	   --------------------------                      ----------RH-- */

#ifndef __RH_PARALLEL_H__
#define __RH_PARALLEL_H__

/* --- Work on the half-open range [begin, end) of items. nt is the
	   index of the worker, in 0..nthreads-1; rh is the state of the
	   caller of parallelFor. A work function that fails records the
	   error with Fail and returns: parallelFor then hands out no more
	   work and returns the recorded status (Q-51) --  -------------- */

typedef void (*ParallelRange)(RHState *rh, int begin, int end, int nt,
							  void *arg);

RH_NODISCARD rh_status parallelFor(RHState *rh, int N, int nthreads, ParallelRange work,
								   void *arg);
RH_NODISCARD rh_status parallelForDynamic(RHState *rh, int N, int chunk, int nthreads,
										  ParallelRange work, void *arg);

#endif /* __RH_PARALLEL_H__ */
