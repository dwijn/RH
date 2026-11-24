/* ------- file: -------------------------- parallel.c --------------
	   --------------------------                      ----------RH-- */

/* --- Splits a loop over N independent items into contiguous chunks
	   and runs them concurrently, one chunk per thread. Meant for
	   loops over grid points, lines, etc. whose iterations do not
	   share state; the chunks are static, so this suits loops whose
	   items cost about the same. The calling thread does chunk 0
	   itself. With nthreads <= 1 (or N <= 1) it simply calls work()
	   once for the whole range, with no threads involved.

	   parallelForDynamic does the same with dynamic scheduling, for
	   loops whose items differ a lot in cost: the workers take
	   contiguous ranges of chunk items, one after the other, from a
	   shared atomic counter until the counter passes N. One worker
	   therefore calls work() many times; per-worker scratch space
	   should be indexed by nt, not by range. With nthreads <= 1 it
	   calls work(0, N, 0, arg) once, i.e. in item order.
	   --                                              -------------- */

#include <stdlib.h>
#include <stdatomic.h>
#include <pthread.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "parallel.h"

extern InputData input;
extern char      messageStr[];

typedef struct {
	ParallelRange work;
	void         *arg;
	int           begin, end, nt;
} Chunk;

static void *runChunk(void *argument)
{
	Chunk *c = (Chunk *) argument;

	c->work(c->begin, c->end, c->nt, c->arg);
	return NULL;
}

void parallelFor(int N, int nthreads, ParallelRange work, void *arg)
{
	int        nt, size, Nchunk;
	Chunk     *chunk;
	pthread_t *id;
	bool_t    *started;

	if (nthreads > N) nthreads = N;
	if (nthreads <= 1) {
		work(0, N, 0, arg);
		return;
	}

	size    = (N + nthreads - 1) / nthreads;
	Nchunk  = (N + size - 1) / size;
	chunk   = (Chunk *) malloc(Nchunk * sizeof(Chunk));
	id      = (pthread_t *) malloc(Nchunk * sizeof(pthread_t));
	started = (bool_t *) calloc(Nchunk, sizeof(bool_t));

	for (nt = 0; nt < Nchunk; nt++) {
		chunk[nt].work  = work;
		chunk[nt].arg   = arg;
		chunk[nt].nt    = nt;
		chunk[nt].begin = nt * size;
		chunk[nt].end   = (nt + 1) * size < N ? (nt + 1) * size : N;
	}
	for (nt = 1; nt < Nchunk; nt++) {
		if (pthread_create(&id[nt], &input.thread_attr, runChunk, &chunk[nt])) {

			/* --- Could not start the thread: do its share here instead --- */

			Error(WARNING, "parallelFor",
				  "Could not create a thread, running its share in the caller");
			runChunk(&chunk[nt]);
		} else
			started[nt] = TRUE;
	}
	runChunk(&chunk[0]);
	for (nt = 1; nt < Nchunk; nt++)
		if (started[nt]) pthread_join(id[nt], NULL);

	free(started);
	free(id);
	free(chunk);
}

/* ------- begin -------------------------- parallelForDynamic.c ---- */

typedef struct {
	ParallelRange work;
	void         *arg;
	int           N, chunk;
	atomic_llong  next;
} DynamicJob;

typedef struct {
	DynamicJob *job;
	int         nt;
} DynamicWorker;

static void *runDynamic(void *argument)
{
	DynamicWorker *w   = (DynamicWorker *) argument;
	DynamicJob    *job = w->job;
	long long      begin;

	/* --- Take the next chunk until none are left. The counter is
		   64-bit so that it cannot overflow when N is close to INT_MAX */

	while ((begin = atomic_fetch_add(&job->next, job->chunk)) < job->N) {
		long long end = begin + job->chunk;

		if (end > job->N) end = job->N;
		job->work((int) begin, (int) end, w->nt, job->arg);
	}
	return NULL;
}

void parallelForDynamic(int N, int chunk, int nthreads,
						ParallelRange work, void *arg)
{
	int            nt, Nchunk;
	DynamicJob     job;
	DynamicWorker *worker;
	pthread_t     *id;
	bool_t        *started;

	if (chunk < 1) chunk = 1;
	Nchunk = (int) (((long long) N + chunk - 1) / chunk);
	if (nthreads > Nchunk) nthreads = Nchunk;
	if (nthreads <= 1) {
		work(0, N, 0, arg);
		return;
	}

	job.work  = work;
	job.arg   = arg;
	job.N     = N;
	job.chunk = chunk;
	atomic_init(&job.next, 0);

	worker  = (DynamicWorker *) malloc(nthreads * sizeof(DynamicWorker));
	id      = (pthread_t *) malloc(nthreads * sizeof(pthread_t));
	started = (bool_t *) calloc(nthreads, sizeof(bool_t));

	for (nt = 0; nt < nthreads; nt++) {
		worker[nt].job = &job;
		worker[nt].nt  = nt;
	}
	for (nt = 1; nt < nthreads; nt++) {
		if (pthread_create(&id[nt], &input.thread_attr, runDynamic, &worker[nt])) {

			/* --- Could not start the thread: the others will simply
				   take its chunks. Worker nt stays unused --  ---------- */

			Error(WARNING, "parallelForDynamic",
				  "Could not create a thread, its chunks go to the others");
		} else
			started[nt] = TRUE;
	}
	runDynamic(&worker[0]);
	for (nt = 1; nt < nthreads; nt++)
		if (started[nt]) pthread_join(id[nt], NULL);

	free(started);
	free(id);
	free(worker);
}
/* ------- end ---------------------------- parallelForDynamic.c ---- */
