/* ------- file: -------------------------- test_parallel.c ---------

	   parallelFor and parallelForDynamic (parallel.c) with a worker
	   that fails (Q-51):

	   -- without a failure every item is done once and the status is
		  RH_OK, with one thread and with several;
	   -- a worker that records an error (Fail) and returns: the loop
		  returns that status, and parallelForDynamic hands out no more
		  work after it (of 100000 items in chunks of 1, far fewer
		  than all are done);
	   -- an error recorded before the call: no work at all, and the
		  recorded status.
	   --                                              -------------- */

#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"
#include "state.h"
#include "check.h"

#define NITEM 100000

typedef struct {
	atomic_int done[NITEM];
	atomic_int count;
	int        failAt; /* the item at which the worker fails, or -1 */
} Job;

static void work(RHState *rh, int begin, int end, int nt, void *arg)
{
	Job *job = (Job *) arg;

	(void) nt;
	for (int n = begin; n < end; n++) {
		if (n == job->failAt) {
			Fail(rh, RH_E_NUMERIC, "work", "item %d failed", n);
			return;
		}
		atomic_fetch_add(&job->done[n], 1);
		atomic_fetch_add(&job->count, 1);
	}
}

static Job *newJob(int failAt)
{
	Job *job = calloc(1, sizeof(Job));

	job->failAt = failAt;
	return job;
}

static int allOnce(Job *job, int N)
{
	for (int n = 0; n < N; n++)
		if (atomic_load(&job->done[n]) != 1) return 0;
	return 1;
}

int main(void)
{
	Job      *job;
	rh_status s;

	pthread_attr_init(&rh->input.thread_attr);

	/* --- No failure --                                 -------------- */

	for (int nthreads = 1; nthreads <= 4; nthreads += 3) {
		job = newJob(-1);
		s   = parallelFor(rh, NITEM, nthreads, work, job);
		CHECK(s == RH_OK && allOnce(job, NITEM), "parallelFor, %d threads: status %d, %d done",
			  nthreads, s, atomic_load(&job->count));
		free(job);

		job = newJob(-1);
		s   = parallelForDynamic(rh, NITEM, 7, nthreads, work, job);
		CHECK(s == RH_OK && allOnce(job, NITEM),
			  "parallelForDynamic, %d threads: status %d, %d done", nthreads, s,
			  atomic_load(&job->count));
		free(job);
	}

	/* --- A worker fails --                             -------------- */

	for (int nthreads = 1; nthreads <= 4; nthreads += 3) {
		job = newJob(10);
		s   = parallelForDynamic(rh, NITEM, 1, nthreads, work, job);
		CHECK(s == RH_E_NUMERIC && strcmp(rh->error.message, "item 10 failed") == 0,
			  "parallelForDynamic, %d threads: status %d", nthreads, s);
		CHECK(atomic_load(&job->count) < NITEM / 10,
			  "parallelForDynamic, %d threads: %d items done after the failure", nthreads,
			  atomic_load(&job->count));
		ClearError(rh);
		free(job);

		job = newJob(NITEM / 2);
		s   = parallelFor(rh, NITEM, nthreads, work, job);
		CHECK(s == RH_E_NUMERIC && Failed(rh), "parallelFor, %d threads: status %d", nthreads, s);
		ClearError(rh);
		free(job);
	}

	/* --- An error recorded before the call --          -------------- */

	Fail(rh, RH_E_DATA, "earlier", "an earlier error");
	job = newJob(-1);
	s   = parallelFor(rh, NITEM, 4, work, job);
	CHECK(s == RH_E_DATA && atomic_load(&job->count) == 0, "parallelFor after an error: %d, %d done",
		  s, atomic_load(&job->count));
	s = parallelForDynamic(rh, NITEM, 7, 4, work, job);
	CHECK(s == RH_E_DATA && atomic_load(&job->count) == 0,
		  "parallelForDynamic after an error: %d, %d done", s, atomic_load(&job->count));
	ClearError(rh);
	free(job);

	return check_summary("test_parallel");
}
