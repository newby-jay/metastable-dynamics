#include "math.h"
#include "gsl/gsl_rng.h"
#include "time.h"
/* #include "gsl/gsl_histogram.h" */
#include "gsl/gsl_math.h"
#include "gsl/gsl_sf_expint.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_odeiv2.h"

#include "iostream.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double absf(double x){
  if (x>=0){return x;}
  else {return -x;}
}
double Eifun(double x){
  if (x > 1000) {
    double v = 1;
    double rv = 1;
    int i = 1;
    // Ei \sim exp(x)/x*sum(n!/x^n, n=0..N-1) error is O(N!x^(-N)) with 6 terms x = 1000 error is O(1e-16)
    for (i=1; i<=6; i++){
      v *= i/x;
      rv += v;
    }
    return exp(x)/x*rv;
  }
  else {
    gsl_set_error_handler_off();
    //cout << "\n loop 1 \n";
    return gsl_sf_expint_Ei(x);
  }
}
void rf_fn(double *tp, double u, double c, double z0, double Eiz0, double q1, double thresh) {
  // Root finding function for computing random times in the MC simulation algorithm
  // -c*(gsl_sf_expint_Ei(z0*exp(-q1*t3)) - Eiz0)/q1 + u3 = 0 goes from u3 to inf
  int nt = 0.;
  double dt, t = *tp, ta = 0; // left side of bracket is ta, right side is t
  double Fa = u - 1, Fb = u - c*(Eifun(z0*exp(-q1*t)) - Eiz0)/q1; // bracket root
  while ((Fb < 0) && (nt < 1e3)) { // if Fup < 0 shift the bracket up
    ta = t;
    Fa = Fb;
    t *= 2;
    Fb = -c*(Eifun(z0*exp(-q1*t)) - Eiz0)/q1 + u;
    nt ++;
  }
  if (nt == 1e3) {cout << "\n loop 1 \n";}
  nt = 0;
  do { // use secant method to get close to the root
    double tr = ta - Fa*(t - ta)/(Fb - Fa);
    double Fr = -c*(Eifun(z0*exp(-q1*tr)) - Eiz0)/q1 + u;
    if (Fr < 0) {
      ta = tr;
      Fa = Fr;
    }
    else {
      t = tr;
      Fb = Fr;
    }
    nt += 1;
  } while (((t - ta)/t > 1) && (nt < 1e3));
  nt = 0;
  do { // Newton method to polish
    double z = z0*exp(-q1*t);
    dt = (-c*(Eifun(z) - Eiz0)/q1 + u)/(c*exp(z0*exp(-q1*t)));
    t -= dt;
    nt += 1;
  } while ((absf(dt) > thresh) && (gsl_finite(t) == 1) && (nt < 1e4));
  /* if (nt == 1e4) {cout << "\n newton loop. t=" << t << " dt=" << absf(dt) << "\n";} */
  *tp = t;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////Sodium and Potassium ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void monte_carlo_2D_traj(double ti, double ni, double mi, double xi, int Npoints,
                         double epsilon, double Nna, double Nk, double gna, double vna, double gk, double vk, double gl, 
                         double vl, double Iapp, double betana, double betak, double thetak1, double thetak2, double thetana1, double thetana2, 
                         double solt[], double solx[], double soln[], double solm[]){
  double thresh = 1e-8;
  //////////////// RNG initialization ////////////////
  struct timeval tv;  // seeding with the current time
  const gsl_rng_type * Trng;
  gettimeofday(&tv, NULL);
  long int sd = tv.tv_sec+tv.tv_usec;
  gsl_rng * r;  // initializtion of the GSL RNG generator
  gsl_rng_env_setup();
  Trng = gsl_rng_default;
  r = gsl_rng_alloc(Trng);
  gsl_rng_set(r, sd);
  //////////////// Notes about the rate functions ////////////////
  // K1 = betana*n
  // K2 = betana*(Nna-n)*exp(4*(thetana1*x + thetana2))
  // K3 = betak*m*exp(thetak1*x + thetak2)
  // K4 = betak*(M - m)*exp(-thetak1*x - thetak2)
  // xdot = -(n/Nna*gna + m/Nk*gk + gl)*x + (n/Nna*gna*vna + m/Nk*gk*vk + gl*vl + Iapp)
  // xdot = -q1(n, m)*x + q2(n, m)
  // x(t) = (x0 - q2/q1)*exp(-q1*t) + q2/q1
  // t = - log((x-q2/q1)/(x0-q2/q1))/q1
  ////////////////////////////////////////////////////////////////
  // initialization of the Gillespie algorithm
  double t=ti, x=xi, n=ni, m=mi;
  int j = 0;
  solt[j] = t;
  solx[j] = x;
  soln[j] = n;
  solm[j] = m;
  for (j=0;j<Npoints;j++) {   // Evolution loop
    double t0=t, x0=x, n0=n, m0=m;
    // generate 4 RVs, 3 of which need newton solve, then pick the smallest time
    double u1 = gsl_rng_uniform_pos(r);
    double u2 = gsl_rng_uniform_pos(r);
    double u3 = gsl_rng_uniform_pos(r);
    double u4 = gsl_rng_uniform_pos(r);
    u1 = log(u1);
    u2 = log(u2);
    u3 = log(u3);
    u4 = log(u4);
    double t1 = -u1/(betana/epsilon*n0);
    // following three are initial guesses for Newton interation, based on time independent rates fixed at the initial time
    double t2 = fmin(-u2/(betana/epsilon*(Nna - n0)*exp(4*(thetana1*x0 + thetana2))), 500);
    double t3 = -u3/(betak*m0*exp(thetak1*x0 + thetak2));
    double t4 = -u4/(betak*(Nk - m0)*exp(-(thetak1*x0 + thetak2)));
    //////////////// Newton iterations ////////////////
    double q1 = n0/Nna*gna + m0/Nk*gk + gl;
    double q2 = n0/Nna*gna*vna + m0/Nk*gk*vk + gl*vl + Iapp;
    // Newton iteration  for t2...
    double c = betana/epsilon*(Nna - n0)*exp(4*(thetana1*q2/q1 + thetana2));
    double z00 = x0 - q2/q1;
    double z0 = 4*thetana1*z00;
    double Eiz0, dt = 0;
    if ((Nna - n0) == 0){t2 = 1./0.;}
    else {
      Eiz0 = Eifun(z0); // note that the maple Ei function and the gsl version differ slightly Ei(n=1, -x) = -Ei(x)
      rf_fn(&t2, u2, c, z0, Eiz0, q1, thresh);
    }
    // Root finding iteration  for t3...
    if (gsl_finite(t3) == 1){ // when m0 = 0, t3 = inf, i.e., reflecting boundary at m=0.  The t3 transition of for m -> m-1.
      c = betak*m0*exp(thetak1*q2/q1 + thetak2);
      z0 = thetak1*z00;
      Eiz0 = Eifun(z0);
      rf_fn(&t3, u3, c, z0, Eiz0, q1, thresh);
    }
    // Root finding iteration for t4...
    if (gsl_finite(t4) == 1) {
      c = betak*(Nk - m0)*exp(-thetak1*q2/q1 - thetak2);
      z0 = -thetak1*z00;
      Eiz0 = Eifun(z0);
      rf_fn(&t4, u4, c, z0, Eiz0, q1, thresh);
    }
    // pick smallest time as next event
    double t12 = fmin(t1, t2);
    double t34 = fmin(t3, t4);
    double tnext = fmin(t12, t34);
    ///////////////////////// update time and state ///////////////////////////////////
    t = t0 + tnext; // update time
    x = (x0 - q2/q1)*exp(-q1*tnext) + q2/q1;
    if (tnext == t1) {n-=1;}
    else if (tnext == t2) {n+=1;}
    else if (tnext == t3) {m-=1;}
    else {m+=1;}
    solt[j] = t;
    solx[j] = x;
    soln[j] = n;
    solm[j] = m;
  }
  gsl_rng_free(r);  // free memory
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void monte_carlo_2D_traj_sep(double ti, double ni, double mi, double xi, int Npoints, double sep[],
                             double epsilon, double Nna, double Nk, double gna, double vna, double gk, double vk, double gl, 
                             double vl, double Iapp, double betana, double betak, double thetak1, double thetak2, double thetana1, double thetana2, 
                             double solt[], double solx[], double soln[], double solm[]){
  double thresh = 1e-8;
  //////////////// RNG initialization ////////////////
  struct timeval tv;  // seeding with the current time
  const gsl_rng_type * Trng;
  gettimeofday(&tv, NULL);
  long int sd = tv.tv_sec+tv.tv_usec;
  gsl_rng * r;  // initializtion of the GSL RNG generator
  gsl_rng_env_setup();
  Trng = gsl_rng_default;
  r = gsl_rng_alloc(Trng);
  gsl_rng_set(r, sd);
  //////////////// Notes about the rate functions ////////////////
  // K = betana*(n + (Nna-n)*exp(thetana1*x + thetana2)) + betak*(m*exp(thetak1*x) + (Nk-m)*exp(thetak2*x))
  // K1 = betana*n
  // K2 = betana*(Nna-n)*exp(thetana1*x + thetana2)
  // K3 = betak*m*exp(thetak1*x)
  // K4 = betak*m*exp(thetak1*x)
  // xdot = -(n/Nna*gna + m/Nk*gk + gl)*x + (n/Nna*gna*vna + m/Nk*gk*vk + gl*vl + Iapp)
  // xdot = -q1(n, m)*x + q2(n, m)
  // x(t) = (x0 - q2/q1)*exp(-q1*t) + q2/q1
  // t = - log((x-q2/q1)/(x0-q2/q1))/q1
  ////////////////////////////////////////////////////////////////
  // initialization of the Gillespie algorithm
  double vside = (xi - sep[int(mi)])/absf(xi - sep[int(mi)]);
  double t=ti, x=xi, n=ni, m=mi;
  int j = 0, exit = 0;
  solt[j] = t;
  solx[j] = x;
  soln[j] = n;
  solm[j] = m;
  do {   // Evolution loop
    double t0=t, x0=x, n0=n, m0=m;
    // generate 4 RVs, 3 of which need newton solve, then pick the smallest time
    double u1 = gsl_rng_uniform_pos(r);
    double u2 = gsl_rng_uniform_pos(r);
    double u3 = gsl_rng_uniform_pos(r);
    double u4 = gsl_rng_uniform_pos(r);
    u1 = log(u1);
    u2 = log(u2);
    u3 = log(u3);
    u4 = log(u4);
    double t1 = -u1/(betana/epsilon*n0);
    // following three are initial guesses for fixed point interation, based on time independent rates fixed at the initial time
    double t2 = fmin(-u2/(betana/epsilon*(Nna - n0)*exp(4*(thetana1*x0 + thetana2))), 500);
    double t3 = fmin(-u3/(betak*m0*exp(thetak1*x0 + thetak2)), 560);
    double t4 = -u4/(betak*(Nk - m0)*exp(-(thetak1*x0 + thetak2)));
    //////////////// Root finding iterations ////////////////
    double q1 = n0/Nna*gna + m0/Nk*gk + gl;
    double q2 = n0/Nna*gna*vna + m0/Nk*gk*vk + gl*vl + Iapp;
    // Root finding iteration  for t2...
    double c = betana/epsilon*(Nna - n0)*exp(4*(thetana1*q2/q1 + thetana2));
    double z00 = x0 - q2/q1;
    double z0 = 4*thetana1*z00;
    double Eiz0, dt = 0;
    if ((Nna - n0) == 0){t2 = 1./0.;}
    else {
      Eiz0 = Eifun(z0); // note that the maple Ei function and the gsl version differ slightly Ei(n=1, -x) = -Ei(x)
      rf_fn(&t2, u2, c, z0, Eiz0, q1, thresh);
    }
    // Root finding iteration  for t3...
    if (gsl_finite(t3) == 1){ // when m0 = 0, t3 = inf, i.e., reflecting boundary at m=0.  The t3 transition of for m -> m-1.
      c = betak*m0*exp(thetak1*q2/q1 + thetak2);
      z0 = thetak1*z00;
      Eiz0 = Eifun(z0);
      rf_fn(&t3, u3, c, z0, Eiz0, q1, thresh);
    }
    // Root finding iteration for t4...
    if (gsl_finite(t4) == 1) {
      c = betak*(Nk - m0)*exp(-thetak1*q2/q1 - thetak2);
      z0 = -thetak1*z00;
      Eiz0 = Eifun(z0);
      rf_fn(&t4, u4, c, z0, Eiz0, q1, thresh);
    }
    // pick smallest time as next event
    double t12 = fmin(t1, t2);
    double t34 = fmin(t3, t4);
    double tnext = fmin(t12, t34);
    ///////////////////////// update time and state ///////////////////////////////////
    t = t0 + tnext; // update time
    x = (x0 - q2/q1)*exp(-q1*tnext) + q2/q1;
    if (tnext == t1) {n-=1;}
    else if (tnext == t2) {n+=1;}
    else if (tnext == t3) {m-=1;}
    else {m+=1;}
    j += 1;
    if (j >= Npoints){j = 0;}
    solt[j] = t;
    solx[j] = x;
    soln[j] = n;
    solm[j] = m;
    double vsep = sep[int(m)];
    if ((x - vsep)*vside < 0){ // check to see if separatrix is crossed
      exit = 1;
      solt[j] = t0 - log((vsep - q2/q1)/(x0-q2/q1))/q1; // the time at which the separatrix was crossed
      solx[j] = vsep;
      soln[j] = n0;
      solm[j] = m0;
    }
  } while (exit == 0);
  gsl_rng_free(r);  // free memory
  return;
}
