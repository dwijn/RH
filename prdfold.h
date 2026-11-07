/* ------- file: -------------------------- prdfold.h --------------- */

/* --- Fine absorption grid of the PRD scattering integrals and the
	   folding of linear interpolation into the weights (prdfold.c) -- */

#ifndef PRDFOLD_H
#define PRDFOLD_H

/* --- Needs rh.h (bool_t) included first, like RH's other headers -- */

int    prdGrid(double q_emit, double waveratio, double *q0, double *qN);
double prdQuadWeight(int lap, int Np);
void   prdFineGrid(double q0, int Np, bool_t symmetric, double *qp);
void   prdFoldIndex(int N, double *qtab, double x, int *j, double *fx);
int    prdFoldRange(int N, double *qtab, double q0, int Np, bool_t symmetric, int *m0);
void   prdFold(int N, double *qtab, int Np, double *qp, double *gw, int m0, int len,
			   double *W);

#endif /* !PRDFOLD_H */
