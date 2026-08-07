/* Nelder-Mead simplex minimization.
 * Implemented from the published algorithm (Nelder & Mead 1965; standard
 * parameter choices as in Lagarias et al. 1998). Contains no code derived
 * from Numerical Recipes or any other implementation. */

#ifndef NELDERMEAD_H
#define NELDERMEAD_H

/* Minimize func over n parameters.
 * x[0..n-1]: starting point on input, best point found on output.
 * step:      initial simplex edge length (added to each coordinate in turn).
 * ftol:      relative tolerance on the spread of function values across the
 *            simplex; convergence when 2|f_worst - f_best| /
 *            (|f_worst| + |f_best| + 1e-10) < ftol.
 * maxeval:   maximum number of function evaluations. If reached before the
 *            tolerance is met, a warning is printed to stderr and the best
 *            point found is returned (this is expected when the optimum lies
 *            on a boundary of the parameter space).
 * fmin_out:  best function value found.
 * Returns 0 on convergence, 1 if maxeval was reached. */
int nelder_mead(double *x, int n, double step, double ftol, int maxeval,
                double (*func)(const double *), double *fmin_out);

#endif /* NELDERMEAD_H */
