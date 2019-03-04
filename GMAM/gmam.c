// Written by Jay Newby 2019
// Send questions and comments to jnewby@ualberta.ca

#include "math.h"
#include "gmam.h"
//using namespace std;


////////////////////////////////////////////////////////////////
//////////////// Aux functions /////////////////////////////////
////////////////////////////////////////////////////////////////
void action_solve(double x, double y, double dx, double dy, double *ppnt,
    double *qpnt, void *pars) {
  // computes dPhidt using variational formula
  // find max_{(p,q) \in R^2 | H(p, q)=0}(lambda*(dx*p + dy*q) - H(p, q))
  double dp, dq, dpmag, thresh_abs = 1e-15, thresh_rel = 1e-10;
  int j = 0, max_iter = 100;
  // compute p and q satisfying Hp(x_bar, y_bar, p, q) = lambda*dx,
  // Hq(x_bar, y_bar, p, q) = lambda*dy, H(x_bar, y_bar, p, q) = 0
  do {
    double p = *ppnt, q = *qpnt; // newton raphson method
    double Hi = H(x, y, p, q, pars);
    double Hpi = Hp(x, y, p, q, pars);
    double Hqi = Hq(x, y, p, q, pars);
    double Hppi = Hpp(x, y, p, q, pars);
    double Hqqi = Hqq(x, y, p, q, pars);
    double Hpqi = Hpq(x, y, p, q, pars);
    double detHpp = Hppi*Hqqi - pow(Hpqi, 2);
    double a = (Hpi*(Hpi*Hqqi - Hqi*Hpqi) + Hqi*(-Hpi*Hpqi + Hqi*Hppi));
    double b = (dx*(dx*Hqqi - dy*Hpqi) + dy*(-dx*Hpqi + dy*Hppi));
    double c = (a - 2*Hi*detHpp)/b;
    if (c < 0){c = 0;}
    double lambda = sqrt(c); // Lagrange mult
    dp = ((lambda*dx - Hpi)*Hqqi - (lambda*dy - Hqi)*Hpqi)/detHpp;
    dq = (-(lambda*dx - Hpi)*Hpqi + (lambda*dy - Hqi)*Hppi)/detHpp;
    dpmag = sqrt(pow(dp, 2) + pow(dq, 2));
    *ppnt += dp;
    *qpnt += dq;
    j += 1;
  } while (
      (j < max_iter)
      && (dpmag > thresh_abs)
      && (dpmag/sqrt(pow(*ppnt, 2) + pow(*qpnt, 2)) > thresh_rel));
}
////////////////////////////////////////////////////////////////
//////////////// geometric minimum action method ///////////////
////////////////////////////////////////////////////////////////
void gmam(double t[], double s[], double x[], double y[], double p[],
    double q[], double phi[], int Npnts, int Niter, double dt, double abs_err,
    double error[], void *pars) {
  //// input initial guess x, y of length Npnts+1
  //// input space for solution t, p, q, phi, k of length Npnts+1
  //// input initial value of hessian at the stable fixed point, Z0
  //// input dt, timestep for iteration
  //// input dt2, timestep for prefactor integration
  //// input Niter the number of iterations
  //// output solution t, x, y, p, q, phi, k
  ////////////////////////////////////////////////////////////////
  // double x0 = x[0], y0 = y[0], xn = x[Npnts], yn = y[Npnts];
  int i, j=0, n;
  double ds = 1./Npnts; // arclength step size
  double dx[Npnts+1], dy[Npnts+1]; // to store dx/ds, dy/ds
  double l[Npnts+1], dl[Npnts+1]; // to store ds/dt, and dl/ds?
  double pi, qi;
  p[0] = 0, q[0] = 0;
  // to store iterated x, y before interpolation
  double xt[Npnts+1], yt[Npnts+1];
  xt[0] = x[0], yt[0] = y[0], xt[Npnts] = x[Npnts], yt[Npnts] = y[Npnts];
  double ds2, convcheck = 0;
  //// first interpolate to equally spaced points
  s[0] = 0;
  for (i=1;i<Npnts+1;i++) { // compute arclength vector s[i]
    // distance between endpoints of the bin (i-1, i)
    double dsi = sqrt(pow(x[i] - x[i-1], 2) + pow(y[i] - y[i-1], 2));
    s[i] = s[i-1] + dsi;
  }
  n = 1;
  ds2 = s[Npnts]*ds; // scale ds by the current total length of the curve
  for (i=1;i<Npnts+1;i++) { // loop over bins (i-1, i)
    double dsi = s[i] - s[i-1];
    // linear interpolation for points in the bin (i-1, i)
    while ((n*ds2 < s[i]) && (n < Npnts)) {
      x[n] = (x[i] - x[i-1])*(n*ds2 - s[i-1])/dsi + x[i-1];
      y[n] = (y[i] - y[i-1])*(n*ds2 - s[i-1])/dsi + y[i-1];
      n ++;
    }
  }
  do {
    convcheck = 0;
    ////////////////////////// step 1: initialization //////////////////////////
    // compute derivatives for step 2
    double l0 = 1;
    for (i=1;i<Npnts;i++){
      dx[i] = (x[i+1] - x[i-1])/(2*ds); // centered difference
      dy[i] = (y[i+1] - y[i-1])/(2*ds);
      pi = p[i-1], qi = q[i-1];
      // Legendre xform with H(p)=0
      action_solve(x[i], y[i], dx[i], dy[i], &pi, &qi, pars);
      p[i] = pi, q[i] = qi; // store gradient of quasipotential
      double Hpi = Hp(x[i], y[i], p[i], q[i], pars);
      double Hqi = Hq(x[i], y[i], p[i], q[i], pars); // elements of dx/dt
      l[i] = (Hpi*dx[i] + Hqi*dy[i])/(dx[i]*dx[i] + dy[i]*dy[i]); // ds/dt
      l0 = l[i];
    }
    l[0] = 3*l[1] - 3*l[2] + l[3];
    l[Npnts] = 3*l[Npnts-1] - 3*l[Npnts-2] + l[Npnts-3]; // lambda end points
    pi = p[Npnts-1], qi = q[Npnts-1];
    action_solve(x[i], y[i], dx[i], dy[i], &pi, &qi, pars);
    p[Npnts] = pi, q[Npnts] = qi;
    for (i=1;i<Npnts;i++) { // compute lambda'
      dl[i] = (l[i+1] - l[i-1])/(2*ds);
    }
    //////////////////////// step 2: implicit time step ////////////////////////
    // a is the off diagonal(s) and b is the main diagonal
    double a[Npnts+1], b[Npnts+1], c[Npnts+1];
    for (i=1;i<Npnts;i++) {
      a[i] = -l[i]*l[i]*dt/(ds*ds);
      c[i] = a[i];
      b[i] = 1 - 2*a[i];
      double Hpxi = Hpx(x[i], y[i], p[i], q[i], pars);
      double Hpyi = Hpy(x[i], y[i], p[i], q[i], pars);
      double Hqxi = Hqx(x[i], y[i], p[i], q[i], pars);
      double Hqyi = Hqy(x[i], y[i], p[i], q[i], pars);
      double Hppi = Hpp(x[i], y[i], p[i], q[i], pars);
      double Hqqi = Hqq(x[i], y[i], p[i], q[i], pars);
      double Hpqi = Hqq(x[i], y[i], p[i], q[i], pars);
      double Hxi = Hx(x[i], y[i], p[i], q[i], pars);
      double Hyi = Hy(x[i], y[i], p[i], q[i], pars);
      // store the rhs in the solution vector initially
      xt[i] = x[i]
        - dt*(l[i]*(Hpxi*dx[i] + Hpyi*dy[i])
          - (Hppi*Hxi + Hpqi*Hyi) - l[i]*dl[i]*dx[i]);
      yt[i] = y[i]
        - dt*(l[i]*(Hqxi*dx[i] + Hqyi*dy[i])
          - (Hpqi*Hxi + Hqqi*Hyi) - l[i]*dl[i]*dy[i]);
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
    xt[Npnts-1] = (xt[Npnts-1] - a[Npnts-1]*xt[Npnts-2])
      /(b[Npnts-1] - a[Npnts-1]*c[Npnts-2]);
    yt[Npnts-1] = (yt[Npnts-1] - a[Npnts-1]*yt[Npnts-2])
      /(b[Npnts-1] - a[Npnts-1]*c[Npnts-2]);
    for (i = Npnts-1; i-- > 1;) { // back substitution sweep
      xt[i] -= c[i]*xt[i+1];
      yt[i] -= c[i]*yt[i+1];
    }
    ///////////////////////// step 3: interpolation ////////////////////////////
    // interpolate solution to equally spaced (in arclength) points
    s[0] = 0;
    for (i=1;i<Npnts+1;i++) { // compute arclength vector s[i]
      // distance between endpoints of the bin (i-1, i)
      double dsi = sqrt(pow(xt[i] - xt[i-1], 2) + pow(yt[i] - yt[i-1], 2));
      s[i] = s[i-1] + dsi;
    }
    n = 1;
    ds2 = s[Npnts]*ds; // scale ds by the current total length of the curve
    for (i=1;i<Npnts+1;i++) { // loop over bins (i-1, i)
      double dsi = s[i] - s[i-1];
      // linear interpolation for points in the bin (i-1, i)
      while ((n*ds2 < s[i]) && (n < Npnts)) {
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
  //cout << "abs_err = " << sqrt(convcheck) << "\n";
  error[0] = sqrt(convcheck);
  //////////////////////////////// recover time ////////////////////////////////
  // time discretization required for prefactor
  t[0] = 0, t[1] = t[0] + ds/(2*l[1]); // should have l[0] = 0, so start at l[1]
  for (i=2;i<Npnts+1;i++) {
    t[i] = t[i-1] + ds*(1/(2*l[i-1]) + 1/(2*l[i]));
  }
  ////////////////////////// compute quasipotential ////////////////////////////
  phi[0] = 0;
  for (i=1;i<Npnts;i++) {
    // quasipotential integration does not depend on a timestep
    phi[i] = phi[i-1] + p[i]*(x[i] - x[i-1]) + q[i]*(y[i] - y[i-1]);
  }
  phi[Npnts] = phi[Npnts-1];
}
