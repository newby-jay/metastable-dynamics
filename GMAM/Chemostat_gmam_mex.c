// Written by Jay Newby 2019
// Send questions and comments to jnewby@ualberta.ca

#include "mex.h"
#include "math.h"
#include "Chemostat.h"
#include "gmam.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  //// input: (x, y, p, q, Npnts, Niter, dt, abs_err, pars)
  //// output: sol, error
  double *x = mxGetPr(prhs[0]);
  double *y = mxGetPr(prhs[1]);
  double *p = mxGetPr(prhs[2]);
  double *q = mxGetPr(prhs[3]);
  double *Npnts0 = mxGetPr(prhs[4]);
  double *Niter0 = mxGetPr(prhs[5]);
  double *dt0 = mxGetPr(prhs[6]);
  double *abs_err0 = mxGetPr(prhs[7]);
  double *pars = mxGetPr(prhs[8]);



  double dt = *dt0;
  double abs_err = *abs_err0;
  int Npnts = *Npnts0;
  int Niter = *Niter0;
  mwSize nnz = Npnts + 1;

  // double x[Npnts + 1], y[Npnts + 1], p[Npnts + 1], q[Npnts + 1];
  // for (int n=0; n<Npnts+1; n++) {
  //   x[n] = x0[n];
  //   y[n] = y0[n];
  //   p[n] = p0[n];
  //   q[n] = q0[n];
  // }


  // 0, 1, 2, 3, 4, 5, 6
  // t, s, x, y, p, q, phi
  plhs[0] = mxCreateDoubleMatrix(nnz, 7, mxREAL);
  // convergence rate
  plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *output = mxGetPr(plhs[0]);
  double *c_rate = mxGetPr(plhs[1]);

  double t[Npnts+1], s[Npnts+1], phi[Npnts+1];
  struct model_params params = getParams(pars);
  gmam(t, s, x, y, p, q, phi, Npnts, Niter, dt, abs_err, c_rate, &params);

  for (int n=0; n<Npnts+1; n++) {
    output[0*nnz + n] = t[n];
    output[1*nnz + n] = s[n];
    output[2*nnz + n] = x[n];
    output[3*nnz + n] = y[n];
    output[4*nnz + n] = p[n];
    output[5*nnz + n] = q[n];
    output[6*nnz + n] = phi[n];
  }

  return;
}
