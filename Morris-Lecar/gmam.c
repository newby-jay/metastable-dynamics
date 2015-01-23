// Written by Jay Newby 2014
// Send questions and comments to newby@math.utah.edu.

#include "math.h"
#include "iostream.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_odeiv2.h"
using namespace std;


struct SlKAP_params {
  // struct to contain the model parameters to pass to various functions
  double epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2;
};
struct HESS_params {
  struct SlKAP_params *model_params;
  double delt, t0;
  double x0, y0, p0, q0;
  double xf, yf, pf, qf;
  double xt0, yt0, pt0, qt0;
  double xtf, ytf, ptf, qtf;
};
////////////////////////////////////////////////////////////////
// Various partial derivatives of the Hamiltonian, H(x, y, p, q).
////////////////////////////////////////////////////////////////
double H(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double tfev = $tfev, minf = $minf;
  double fna = $fna;
  double g = $g;
  double h = $h;
  /* double sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N)); */
  return $H;
}
double Hp(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double N = (params->N), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double tfev = $tfev, minf = $minf;
  double fna = $fna;
  double g = $g;
  return $Hp;
}
double Hq(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), M = (params->M), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2); // all the paremeters
  double hgam = 1./(epsilon*M);
  return $hq;
}
double Hpp(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double N = (params->N), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double tfev = $tfev, minf = $minf;
  double fna = $fna;
  double g = $g;
  return $Hpp;
}
double Hqq(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), M = (params->M), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2); // all the paremeters
  double hgam = 1./(epsilon*M);
  return $hqq;
}
double Hx(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double tfev = $tfev, minf = $minf, dminf = $dminf;
  double fna = $fna, dfna = -gna;
  double g = $g, gv = $gv;
  double hv = $hv;
  return $Hx;
}
double Hy(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double tfev = $tfev, minf = $minf;
  double fna = $fna;
  double g = $g, gw = $gw;
  double hw = $hw;
  return $Hy;
}
double Hpx(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double N = (params->N), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double tfev = $tfev, minf = $minf, dminf = $dminf;
  double fna = $fna, dfna = -gna;
  double g = $g, gv = $gv;
  /* double sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N)); */
  return $Hpx;
}
double Hpy(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double N = (params->N), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double tfev = $tfev, minf = $minf;
  double fna = $fna;
  double g = $g, gw = $gw;
  return $Hpy;
}
double Hqx(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), M = (params->M), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2); // all the paremeters
  double hgam = 1./(epsilon*M);
  return $hqv;
}
double Hqy(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), M = (params->M), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2); // all the paremeters
  double hgam = 1./(epsilon*M);
  return $hqw;
}
double Hxx(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double tfev = $tfev, minf = $minf, dminf = $dminf, ddminf = $ddminf;
  double fna = $fna, dfna = $dfna;
  double g = $g, gv = $gv;
  double hvv = $hvv;
  return $Hxx;
}
double Hxy(double v, double w, double p, double q, void *pars){
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2);
  double hgam = 1./(epsilon*M);
  double tfev = $tfev, minf = $minf, dminf = $dminf;
  double fna = $fna, dfna = $dfna;
  double g = $g, gv = $gv, gw = $gw, gvw = $gvw;
  double hvw = $hvw;
  return $Hxy;
}
////////////////////////////////////////////////////////////////
//////////////// Aux functions /////////////////////////////////
////////////////////////////////////////////////////////////////
void action_solve(double x, double y, double dx, double dy, double *ppnt, double *qpnt, void *pars){
  // computes dPhidt using variational formula
  // find max_{(p,q) \in R^2 | H(p, q)=0}(lambda*(dx*p + dy*q) - H(p, q))
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double v = x, w = y;
  // to be filled in with template
  double minf = (1. + tanh(2.*(thetana1*v + thetana2)))/2;
  double fna = $fna;
  double g = $g;
  double dp, dq, dpmag, thresh_abs = 1e-15, thresh_rel = 1e-10;
  int j = 0, max_iter = 100;
  do { // compute p and q satisfying Hp(x_bar, y_bar, p, q) = lambda*dx, Hq(x_bar, y_bar, p, q) = lambda*dy, H(x_bar, y_bar, p, q) = 0
    double p = *ppnt, q = *qpnt; // newton raphson method
    double h = $h; // more template stuff
    double Hpp =  $Hpp;
    double Hqq =  $hqq;
    double H = $H;
    double Hp = $Hp;
    double Hq = $hq;
    double a = pow(Hp, 2)/Hpp + pow(Hq, 2)/Hqq;
    double b = pow(dx, 2)/Hpp + pow(dy, 2)/Hqq;
    double c = (a - 2*H)/b;
    if (c < 0){c = 0;}
    double lambda = sqrt(c); // Lagrange mult
    dp = (lambda*dx - Hp)/Hpp;
    dq = (lambda*dy - Hq)/Hqq;
    dpmag = sqrt(pow(dp, 2) + pow(dq, 2));
    *ppnt += dp;
    *qpnt += dq;
    j += 1;
  } while((j < max_iter) && (dpmag > thresh_abs) && (dpmag/sqrt(pow(*ppnt, 2) + pow(*qpnt, 2)) > thresh_rel));
}
////////////////////////////////////////////////////////////////
//////////////// geometric minimum action method ///////////////
////////////////////////////////////////////////////////////////
void gmam(double t[], double s[], double x[], double y[], double p[], double q[], double phi[], int Npnts, int Niter, double dt, double abs_err, double pars[]) {
  // input initial guess x, y of length Npnts+1
  // input space for solution t, p, q, phi, k of length Npnts+1
  // input initial value of hessian at the stable fixed point, Z0
  // input dt, timestep for iteration
  // input dt2, timestep for prefactor integration
  // input Niter the number of iterations
  // output solution t, x, y, p, q, phi, k
  double epsilon = pars[0], N = pars[1], M = pars[2], gna = pars[3], vna = pars[4], gk = pars[5], vk = pars[6], gl = pars[7], vl = pars[8], Iapp = pars[9], betak = pars[10], thetak1 = pars[11], thetak2 = pars[12], thetana1 = pars[13], thetana2 = pars[14]; // unpack parameters
  // double x0 = x[0], y0 = y[0], xn = x[Npnts], yn = y[Npnts];
  int i, j=0, n;
  double ds = 1./double(Npnts); // arclength step size
  struct SlKAP_params params = {epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2}; // parameter struct
  double dx[Npnts+1], dy[Npnts+1]; // to store dx/ds, dy/ds
  double l[Npnts+1], dl[Npnts+1]; // to store ds/dt, and dl/ds?
  double pi, qi;
  p[0] = 0, q[0] = 0;
  double xt[Npnts+1], yt[Npnts+1]; // to store iterated x, y before interpolation
  xt[0] = x[0], yt[0] = y[0], xt[Npnts] = x[Npnts], yt[Npnts] = y[Npnts];
  double ds2, convcheck = 0;
  //// first interpolate to equally spaced points
  s[0] = 0;
  for (i=1;i<Npnts+1;i++) { // compute arclength vector s[i]
    double dsi = sqrt(pow(x[i] - x[i-1], 2) + pow(y[i] - y[i-1], 2)); // distance between endpoints of the bin (i-1, i)
    s[i] = s[i-1] + dsi;
  }
  n = 1;
  ds2 = s[Npnts]*ds; // scale ds by the current total length of the curve
  for (i=1;i<Npnts+1;i++) { // loop over bins (i-1, i)
    double dsi = s[i] - s[i-1];
    while ((n*ds2 < s[i]) && (n < Npnts)) { // linear interpolation for points in the bin (i-1, i)
      x[n] = (x[i] - x[i-1])*(n*ds2 - s[i-1])/dsi + x[i-1];
      y[n] = (y[i] - y[i-1])*(n*ds2 - s[i-1])/dsi + y[i-1];
      n ++;
    }
  }
  do {
    convcheck = 0;
    //////////////////////////////// step 1: initialization ///////////////////////////////
    // compute derivatives for step 2
    double l0 = 1;
    for (i=1;i<Npnts;i++){
      dx[i] = (x[i+1] - x[i-1])/(2*ds); // centered difference
      dy[i] = (y[i+1] - y[i-1])/(2*ds);
      pi = p[i-1], qi = q[i-1];
      action_solve(x[i], y[i], dx[i], dy[i], &pi, &qi, &params); // Legendre xform with H(p)=0 solve
      p[i] = pi, q[i] = qi; // store gradient of quasipotential
      double Hpi = Hp(x[i], y[i], p[i], q[i], &params), Hqi = Hq(x[i], y[i], p[i], q[i], &params); // elements of dx/dt
      l[i] = (Hpi*dx[i] + Hqi*dy[i])/(dx[i]*dx[i] + dy[i]*dy[i]); // ds/dt
      l0 = l[i];
    }
    l[0] = 3*l[1] - 3*l[2] + l[3];
    l[Npnts] = 3*l[Npnts-1] - 3*l[Npnts-2] + l[Npnts-3]; // lambda end points
    pi = p[Npnts-1], qi = q[Npnts-1];
    action_solve(x[i], y[i], dx[i], dy[i], &pi, &qi, &params);
    p[Npnts] = pi, q[Npnts] = qi;
    for (i=1;i<Npnts;i++) { // compute lambda'
      dl[i] = (l[i+1] - l[i-1])/(2*ds);
    }
    //////////////////////////////// step 2: implicit time step ///////////////////////////
    double a[Npnts+1], b[Npnts+1], c[Npnts+1]; // a is the off diagonal(s) and b is the main diagonal
    for (i=1;i<Npnts;i++) {
      a[i] = -l[i]*l[i]*dt/(ds*ds);
      c[i] = a[i];
      b[i] = 1 - 2*a[i];
      double Hpxi = Hpx(x[i], y[i], p[i], q[i], &params), Hpyi = Hpy(x[i], y[i], p[i], q[i], &params), Hqxi = Hqx(x[i], y[i], p[i], q[i], &params);
      double Hqyi = Hqy(x[i], y[i], p[i], q[i], &params), Hppi = Hpp(x[i], y[i], p[i], q[i], &params);
      double Hqqi = Hqq(x[i], y[i], p[i], q[i], &params), Hxi = Hx(x[i], y[i], p[i], q[i], &params), Hyi = Hy(x[i], y[i], p[i], q[i], &params);
      xt[i] = x[i] - dt*(l[i]*(Hpxi*dx[i] + Hpyi*dy[i]) - Hppi*Hxi - l[i]*dl[i]*dx[i]); // store the rhs in the solution vector initially
      yt[i] = y[i] - dt*(l[i]*(Hqxi*dx[i] + Hqyi*dy[i]) - Hqqi*Hyi - l[i]*dl[i]*dy[i]);
    }
    xt[1] -= a[1]*x[0];
    yt[1] -= a[1]*y[0]; // boundary condition
    xt[Npnts-1] -= a[Npnts-1]*x[Npnts];
    yt[Npnts-1] -= a[Npnts-1]*y[Npnts]; // boundary condition
    // begin Thomas algorithm
    c[1] = c[1]/b[1];
    xt[1] = xt[1]/b[1];
    yt[1] = yt[1]/b[1];
    for (i=2;i<Npnts-1;i++) { // forward sweep
      double m = 1./(b[i] - a[i]*c[i-1]);
      c[i] *= m;
      xt[i] = (xt[i] - a[i]*xt[i-1])*m;
      yt[i] = (yt[i] - a[i]*yt[i-1])*m;
    }
    xt[Npnts-1] = (xt[Npnts-1] - a[Npnts-1]*xt[Npnts-2])/(b[Npnts-1] - a[Npnts-1]*c[Npnts-2]);
    yt[Npnts-1] = (yt[Npnts-1] - a[Npnts-1]*yt[Npnts-2])/(b[Npnts-1] - a[Npnts-1]*c[Npnts-2]);
    for (i = Npnts-1; i-- > 1;) { // back substitution sweep
      xt[i] -= c[i]*xt[i+1];
      yt[i] -= c[i]*yt[i+1];
    }
    //////////////////////////////// step 3: interpolation ////////////////////////////////
    // interpolate solution to equally spaced (in arclength) points
    s[0] = 0;
    for (i=1;i<Npnts+1;i++) { // compute arclength vector s[i]
      double dsi = sqrt(pow(xt[i] - xt[i-1], 2) + pow(yt[i] - yt[i-1], 2)); // distance between endpoints of the bin (i-1, i)
      s[i] = s[i-1] + dsi;
    }
    n = 1;
    ds2 = s[Npnts]*ds; // scale ds by the current total length of the curve
    for (i=1;i<Npnts+1;i++) { // loop over bins (i-1, i)
      double dsi = s[i] - s[i-1];
      while ((n*ds2 < s[i]) && (n < Npnts)) { // linear interpolation for points in the bin (i-1, i)
        double s3 = (n*ds2 - s[i-1]);
        double xn = (xt[i] - xt[i-1])*s3/dsi + xt[i-1];
        double yn = (yt[i] - yt[i-1])*s3/dsi + yt[i-1];
        double convcheckn = pow(x[n] - xn, 2) + pow(y[n] - yn, 2);
        if (convcheck < convcheckn) {convcheck = convcheckn;}
        x[n] = xn, y[n] = yn;
        n ++;
      }
    }
    j++;
  } while ((j < Niter) && (sqrt(convcheck) > abs_err));
  cout << "abs_err = " << sqrt(convcheck) << "\n";
  //////////////////////////////// recover time ////////////////////////////////
  // time discretization required for prefactor
  t[0] = 0, t[1] = t[0] + ds/(2*l[1]); // should have l[0] = 0, so start at l[1]
  for (i=2;i<Npnts+1;i++) {
    t[i] = t[i-1] + ds*(1/(2*l[i-1]) + 1/(2*l[i]));
  }
  //////////////////// compute quasipotential ///////////////////////////////////
  phi[0] = 0;
  for (i=1;i<Npnts;i++) {
    phi[i] = phi[i-1] + p[i]*(x[i] - x[i-1]) + q[i]*(y[i] - y[i-1]); // quasipotential integration does not depend on a timestep
  }
  phi[Npnts] = phi[Npnts-1];
}
