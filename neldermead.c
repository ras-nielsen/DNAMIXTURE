/* Nelder-Mead simplex minimization.
 * Implemented from the published algorithm (Nelder & Mead 1965; standard
 * parameter choices as in Lagarias et al. 1998). Contains no code derived
 * from Numerical Recipes or any other implementation. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "neldermead.h"

#define NM_RHO   1.0   /* reflection coefficient */
#define NM_CHI   2.0   /* expansion coefficient */
#define NM_GAMMA 0.5   /* contraction coefficient */
#define NM_SIGMA 0.5   /* shrink coefficient */

int nelder_mead(double *x, int n, double step, double ftol, int maxeval,
                double (*func)(const double *), double *fmin_out)
{
        int npts = n + 1;
        int i, j, ilo, ihi, inhi, neval = 0, converged = 0;
        double *v, *f, *cent, *xr, *xc;
        double fr, fc, rtol;

        v    = (double *)malloc((size_t)npts * n * sizeof(double)); /* vertices, row i = v + i*n */
        f    = (double *)malloc((size_t)npts * sizeof(double));
        cent = (double *)malloc((size_t)n * sizeof(double));
        xr   = (double *)malloc((size_t)n * sizeof(double));
        xc   = (double *)malloc((size_t)n * sizeof(double));
        if (!v || !f || !cent || !xr || !xc) {
                fprintf(stderr, "Error: nelder_mead: allocation failure\n");
                exit(1);
        }

        /* Initial simplex: the starting point, plus one vertex displaced by
         * step along each coordinate axis. */
        for (i = 0; i < npts; i++) {
                for (j = 0; j < n; j++)
                        v[i*n + j] = x[j];
                if (i > 0)
                        v[i*n + (i-1)] += step;
                f[i] = func(&v[i*n]);
                neval++;
        }

        for (;;) {
                /* Identify best (ilo), worst (ihi), second-worst (inhi). */
                ilo = 0;
                for (i = 1; i < npts; i++)
                        if (f[i] < f[ilo]) ilo = i;
                ihi = (ilo == 0) ? 1 : 0;
                for (i = 0; i < npts; i++)
                        if (f[i] > f[ihi]) ihi = i;
                inhi = (ihi == 0) ? 1 : 0;
                for (i = 0; i < npts; i++)
                        if (i != ihi && f[i] > f[inhi]) inhi = i;

                rtol = 2.0 * fabs(f[ihi] - f[ilo]) /
                       (fabs(f[ihi]) + fabs(f[ilo]) + 1e-10);
                if (rtol < ftol) {
                        converged = 1;
                        break;
                }
                if (neval >= maxeval) {
                        fprintf(stderr, "Warning: Nelder-Mead: maximum number of function "
                                "evaluations (%d) reached before convergence tolerance; "
                                "returning best point found. This is expected when the "
                                "optimum lies on a boundary of the parameter space.\n",
                                maxeval);
                        break;
                }

                /* Centroid of all vertices except the worst. */
                for (j = 0; j < n; j++) {
                        cent[j] = 0.0;
                        for (i = 0; i < npts; i++)
                                if (i != ihi)
                                        cent[j] += v[i*n + j];
                        cent[j] /= n;
                }

                /* Reflection. */
                for (j = 0; j < n; j++)
                        xr[j] = cent[j] + NM_RHO * (cent[j] - v[ihi*n + j]);
                fr = func(xr);
                neval++;

                if (fr < f[ilo]) {
                        /* Expansion. */
                        for (j = 0; j < n; j++)
                                xc[j] = cent[j] + NM_CHI * (xr[j] - cent[j]);
                        fc = func(xc);
                        neval++;
                        if (fc < fr) {
                                for (j = 0; j < n; j++) v[ihi*n + j] = xc[j];
                                f[ihi] = fc;
                        } else {
                                for (j = 0; j < n; j++) v[ihi*n + j] = xr[j];
                                f[ihi] = fr;
                        }
                } else if (fr < f[inhi]) {
                        /* Accept the reflected point. */
                        for (j = 0; j < n; j++) v[ihi*n + j] = xr[j];
                        f[ihi] = fr;
                } else {
                        /* Contraction: outside if the reflected point improves
                         * on the worst vertex, inside otherwise. */
                        if (fr < f[ihi]) {
                                for (j = 0; j < n; j++)
                                        xc[j] = cent[j] + NM_GAMMA * (xr[j] - cent[j]);
                        } else {
                                for (j = 0; j < n; j++)
                                        xc[j] = cent[j] - NM_GAMMA * (cent[j] - v[ihi*n + j]);
                        }
                        fc = func(xc);
                        neval++;
                        if (fc < fr || (fr >= f[ihi] && fc < f[ihi])) {
                                for (j = 0; j < n; j++) v[ihi*n + j] = xc[j];
                                f[ihi] = fc;
                        } else {
                                /* Shrink all vertices toward the best. */
                                for (i = 0; i < npts; i++) {
                                        if (i == ilo) continue;
                                        for (j = 0; j < n; j++)
                                                v[i*n + j] = v[ilo*n + j] +
                                                        NM_SIGMA * (v[i*n + j] - v[ilo*n + j]);
                                        f[i] = func(&v[i*n]);
                                        neval++;
                                }
                        }
                }
        }

        for (j = 0; j < n; j++)
                x[j] = v[ilo*n + j];
        *fmin_out = f[ilo];

        free(v); free(f); free(cent); free(xr); free(xc);
        return converged ? 0 : 1;
}
