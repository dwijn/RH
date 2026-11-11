/* ------- file: -------------------------- test_ng.c ---------------

	   Ng acceleration (Accelerate in accelerate.c) and its positivity
	   check (N-33).

	   On a linear contraction x_{k+1} = lambda x_k + (1 - lambda) x*
	   with one lambda for every component, Ng of order 1 is exact: its
	   first extrapolation lands on x*. The test checks:

	   -- x* positive: the step is accepted (Accelerate returns TRUE),
		  the solution is x* and the history holds it;
	   -- x* with a negative component while G's value there is still
		  positive: the step is rejected (returns FALSE, one value
		  counted), the solution and the history keep G's value, and
		  NgStatus says so;
	   -- a component where G's value is itself negative: not counted,
		  so the step is accepted;
	   -- the step before the delay and the messages of NgStatus.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "accelerate.h"
#include "check.h"

#define N 2

/* --- Runs the contraction from x0 to x* with lambda = 0.5 through
	   Ng of order 1 (delay 3, the minimum): NgInit with x0, then
	   Accelerate on x1 (no step) and x2 (the step). Returns the
	   result of the second call, x2 in g, the result in x -- --- */

static bool_t run(const double x0[N], const double xstar[N], double x[N],
				  double g[N], struct Ng **Ngs)
{
	const double lambda = 0.5;
	bool_t       accel;
	int          k;

	for (k = 0; k < N; k++) x[k] = x0[k];
	*Ngs = NgInit(N, 0, 1, 1, x);

	for (k = 0; k < N; k++) x[k] = lambda * x[k] + (1.0 - lambda) * xstar[k];
	accel = Accelerate(*Ngs, x);
	CHECK(!accel && (*Ngs)->rejected == 0, "no step before the delay");

	for (k = 0; k < N; k++) {
		x[k] = lambda * x[k] + (1.0 - lambda) * xstar[k];
		g[k] = x[k];
	}
	return Accelerate(*Ngs, x);
}

/* --- The last solution in the history (what MaxChange compares) -- */

static double *last(struct Ng *Ngs)
{
	return Ngs->previous[(Ngs->count - 1) % (Ngs->Norder + 2)];
}

int main(void)
{
	struct Ng *Ngs;
	double     x[N], g[N];
	bool_t     accel;
	int        k;

	/* --- x* positive: accepted, exact --             -------------- */
	{
		const double x0[N] = {3.0, 5.0}, xstar[N] = {1.0, 2.0};

		accel = run(x0, xstar, x, g, &Ngs);
		CHECK(accel, "positive limit: step accepted");
		CHECK(Ngs->rejected == 0, "rejected = %ld", Ngs->rejected);
		for (k = 0; k < N; k++) {
			CHECK_CLOSE(x[k], xstar[k], 1e-14, 0.0, "x[%d] = x*", k);
			CHECK(last(Ngs)[k] == x[k], "history[%d] holds the step", k);
		}
		CHECK(strcmp(NgStatus(Ngs, accel), " (accelerated)\n") == 0,
			  "status \"%s\"", NgStatus(Ngs, accel));
		NgFree(Ngs);
	}
	/* --- x*[1] negative, G's value positive: rejected -- ---------- */
	{
		const double x0[N] = {3.0, 1.0}, xstar[N] = {1.0, -0.2};

		accel = run(x0, xstar, x, g, &Ngs);
		CHECK(!accel, "negative limit: step rejected");
		CHECK(Ngs->rejected == 1, "rejected = %ld, want 1", Ngs->rejected);
		for (k = 0; k < N; k++) {
			CHECK(x[k] == g[k], "x[%d] = %g, G's value %g kept", k, x[k], g[k]);
			CHECK(last(Ngs)[k] == g[k], "history[%d] holds G's value", k);
		}
		CHECK(strstr(NgStatus(Ngs, accel), "rejected: 1 value not positive") != NULL,
			  "status \"%s\"", NgStatus(Ngs, accel));

		NgFree(Ngs);
	}
	/* --- Two values going negative: counted, plural -- ------------ */
	{
		const double x0[N] = {1.0, 1.0}, xstar[N] = {-0.2, -0.2};

		accel = run(x0, xstar, x, g, &Ngs);
		CHECK(!accel && Ngs->rejected == 2, "rejected = %ld, want 2", Ngs->rejected);
		CHECK(strstr(NgStatus(Ngs, accel), "rejected: 2 values not positive") != NULL,
			  "status \"%s\"", NgStatus(Ngs, accel));
		NgFree(Ngs);
	}
	/* --- G's value itself negative: not counted, accepted -- ------- */
	{
		const double x0[N] = {3.0, -1.0}, xstar[N] = {1.0, -2.0};

		accel = run(x0, xstar, x, g, &Ngs);
		CHECK(accel, "negative G value: step accepted");
		CHECK(g[1] < 0.0, "G's value %g is negative", g[1]);
		for (k = 0; k < N; k++)
			CHECK_CLOSE(x[k], xstar[k], 1e-14, 0.0, "x[%d] = x*", k);
		NgFree(Ngs);
	}
	/* --- No step: empty status --                     -------------- */
	{
		const double x0[N] = {3.0, 5.0};

		for (k = 0; k < N; k++) x[k] = x0[k];
		Ngs   = NgInit(N, 10, 1, 1, x);
		accel = Accelerate(Ngs, x);
		CHECK(!accel && strcmp(NgStatus(Ngs, accel), "\n") == 0,
			  "status before the delay \"%s\"", NgStatus(Ngs, accel));
		NgFree(Ngs);
	}
	return check_summary("test_ng");
}

/* ------- end ---------------------------- test_ng.c --------------- */
