// Written by Jay Newby 2013
// Send questions and comments to newby@math.utah.edu.

#include "math.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_min.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_sort_double.h"
#include "gsl/gsl_sort_vector.h"
#include "gsl/gsl_errno.h"
#include "iostream.h"
using namespace std;

ostream& clear_line(ostream& out, int n) {
  return out<<'\r'<< string(n, ' ') <<'\r';
}
struct SlKAP_params {
  double epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2;
};
struct minfun_params {
  struct SlKAP_params *model_params;
  double x, y, x1, x2, y1, y2, p1, p2, q1, q2, phi1, phi2;
  double *p, *q;
};
void action_solve(double x, double y, double dx, double dy, double *ppnt, double *qpnt, void *pars) {
  // computes p, q using variational formula
  // find max_{(p,q) \in R^2 | H(p, q)=0}(lambda*(dx*p + dy*q) - H(p, q))
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2); // all the paremeters
  double hgam = 1./(epsilon*M);
  double v = x, w = y;
  // to be filled in with template
  double minf = (1. + tanh(2.*(thetana1*v + thetana2)))/2;
  double fna = $fna;
  double g = $g;
  double dp, dq, dpmag, thresh_rel = 1e-12, thresh_abs = 1e-10;
  int j = 0, max_iter = 50;
  do { // compute p and q satisfying Hp(x_bar, y_bar, p, q) = lambda*dx, Hq(x_bar, y_bar, p, q) = lambda*dy, H(x_bar, y_bar, p, q) = 0
    double p = *ppnt, q = *qpnt; // newton raphson method
    double h = $h; // more template stuff
    double sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N));
    double Hpp =  $Hpp;
    double Hqq =  $hqq;
    double H = $H;
    double Hp = $Hp;
    double Hq = $hq;
    double a = pow(Hp, 2) + pow(Hq, 2)*Hpp/Hqq;
    double b = pow(dx, 2) + pow(dy, 2)*Hpp/Hqq;
    double c = (a - 2*H*Hpp)/b;
    if (c < 0){c = 0;}
    double lambda = sqrt(c); // Lagrange mult
    dp = (lambda*dx - Hp)/Hpp;
    dq = (lambda*dy - Hq)/Hqq;
    dpmag = sqrt(pow(dp, 2) + pow(dq, 2));
    double dcheck = 10.;
    if (dpmag > dcheck) { // relaxation for stability
      double omega0 = 10.;
      double omega = omega0/(omega0 + dpmag - dcheck);
      dp = omega*dp; // if dpmag is too large reduce the magnitude of the update
      dq = omega*dq;
      dpmag = sqrt(pow(dp, 2) + pow(dq, 2));
    }
    *ppnt += dp;
    *qpnt += dq;
    j += 1;
  } while((gsl_finite(*ppnt + *qpnt) == 1) && (j < max_iter) && // absolute error
          (dpmag/sqrt(pow(*ppnt, 2) + pow(*qpnt, 2)) > thresh_rel)); // relative error
}
double minfun(double theta, void *pars){
  // function to minimize, i.e., compute min_{theta}(phi_bar + p*dx + q*dy), where p and q are determined by the 'action_solve' function
  struct minfun_params *params = (struct minfun_params *)pars; // struct has all the values at update point and the two AF points
  struct SlKAP_params *model_params = (params->model_params); // ...and the model prameters
  double x = (params->x), y = (params->y), x1 = (params->x1), x2 = (params->x2), y1 = (params->y1), y2 = (params->y2), p1 = (params->p1), p2 = (params->p2), q1 = (params->q1), q2 = (params->q2), phi1 = (params->phi1), phi2 = (params->phi2);
  double x_bar = (theta*x1 + (1. - theta)*x2), y_bar = (theta*y1 + (1. - theta)*y2); // computing average point based on theta
  double phi_bar = theta*phi1 + (1 - theta)*phi2;
  double *p = (params->p), *q = (params->q); // pointers to p, q
  //xtheta and ytheta undefined
  double dx = x - x_bar, dy = y - y_bar; // update differential
  *p = theta*p1 + (1 - theta)*p2, *q = theta*q1 + (1 - theta)*q2;  //initial guess for p, q
  // stability depends on the initial guess
  double pmag = sqrt(pow(*p, 2) + pow(*q, 2)); // make the initial guess a unit vector
  if (pmag == 0){
    *p = 1, *q = 1;
  }
  else {
    *p /= pmag, *q /= pmag;
  }
  action_solve(x, y, dx, dy, p, q, model_params); // p and q must satisfy stuff
  return phi_bar + (*p)*dx + (*q)*dy; // return the proposed value of the action
}
void bracket_interval(int *init_status, double *a, double *b, double *theta, void *pars){
  // bracket the minimum before polishing with a faster GSL minimization routine
  // interval is (a, b)
  // min arg is theta
  int j;
  struct minfun_params *mfparams = (struct minfun_params *)pars;
  *init_status = 0;
  struct SlKAP_params *params = (mfparams->model_params);
  double epsilon = (params->epsilon), N = (params->N), M = (params->M), gna = (params->gna), vna = (params->vna), gk = (params->gk), vk = (params->vk), gl = (params->gl), vl = (params->vl), Iapp = (params->Iapp), betak = (params->betak), thetak1 = (params->thetak1), thetak2 = (params->thetak2), thetana1 = (params->thetana1), thetana2 = (params->thetana2);
  double x = (mfparams->x), y = (mfparams->y), x1 = (mfparams->x1), x2 = (mfparams->x2), y1 = (mfparams->y1), y2 = (mfparams->y2);
  double hgam = 1./(epsilon*M);
  double v = x, w = y;
  // to be filled in with template
  double tfev = $tfev;
  double minf = $minf;
  double fna = $fna;
  double g = $g;
  double fna2g = fna + 2*g; // these terms show up often in H, Hp, and Hpp
  double binf = 1 - minf;
  double Iinf = minf*fna + g;
  double gfnag = g*(fna + g);
  double c1 = fna2g/2; // these terms are for phat and show up in Hp
  double c3 = N/4*2.*(pow(fna2g/N, 2) - 4/N/N*gfnag); // N/4*2*(fna2g/N*fna2g/N - 4/N/2*2*gfnag/N)
  double c6 = -4/N*gfnag/N + pow(fna2g/N, 2);
  double alphaa = min(c1 - c3/sqrt(c6), c1 + c3/sqrt(c6)), alphab = max(c1 - c3/sqrt(c6), c1 + c3/sqrt(c6));
  double theta0 = (x - x2)/(x1 - x2); // should not be possible for theta0 = 0/0
  double dx0 = x - x2, dy0 = y - y2;
  double dx1 = x - x1, dy1 = y - y1;
  if (((dx0 == 0) && (dx1 == 0)) || ((dy0 == 0) && (dy1 == 0)) || (dy0/dx0 == dy1/dx1)) { // the three points should not be on a single line
    if ((pow(dx1, 2) + pow(dy1, 2)) < (pow(dx0, 2) + pow(dy0, 2))){
      *theta = 1;
    }
    else {*theta = 0;}
    return;
  }
  //// can refine the interval before we even start by considering the sign of dx and the interval (alphaa, alphab)
  if (alphaa > 0){
    if ((0 < theta0) && (theta0 < 1)) {
      if (dx0 < 0) {*a = theta0;} // dx must be positive
      else {*b = theta0;}
    }
    else if (dx0 < 0) { // no soln possible, move on
      return;
    }
  }
  if (alphab < 0) {
    if ((0 < theta0) && (theta0 < 1)) {
      if (dx0 < 0) {*b = theta0;} // dx must be negative
      else {*a = theta0;}
    }
    else if (dx0 > 0) { // no soln possible, move on
      return;
    }
  }
  // Bracket the interval
  double iter = 0, max_iter = 50; // maximum # of iterations of the slower bracketing algorithm
                                  // if this is too low, holes can show up unevaluated
  *theta = (*a) + ((*b) - (*a))*0.3189660; // initial guess, golden ratio is fastest possible.  0<theta<1 is used as weighted average of the two AF points for the fd.  It is chosen based on minimizing the action over theta
  double fa, fb, fm;
  double dxa = x - ((*a)*x1 + (1 - (*a))*x2);
  double dxb = x - ((*b)*x1 + (1 - (*b))*x2);
  fa = minfun(*a, mfparams);
  fb = minfun(*b, mfparams);
  do { // bracket the interval containing the minimum, if it isn't located after <max_iter> iterations, give up
    iter += 1;
    fm = minfun(*theta, mfparams);
    if (fa < fb) { // check for the smaller end point
      if (fm < fa) { // check to see if fm < smaller endpoint
        *init_status = 1;
        return;
      }
      else { // if not, set larger end point to be fm, and let the new guess by (a+b)/2
        fb = fm;
        *b = (*theta);
        *theta = (*a) + ((*b) - (*a))*0.3189660;
      }
    }
    else {
      if (fm < fb) {
        *init_status = 1;
        return;
      }
      else {
        fa = fm;
        *a = (*theta);
        *theta = (*a) + ((*b) - (*a))*0.3189660;
      }
    }
  } while (iter < max_iter);
}
void upwindfd(int ij, int ij1, int ij2, int nx, int ny, 
              double *xpoints, double *ypoints, double *phic, double *phi, 
              double *p, double *q, double *dx, double *dy, void *pars) {
  // the upwind finite difference approximation for phi
  struct SlKAP_params *params = (struct SlKAP_params *)pars;
  int i = (ij/nx); // compute square index from linear index
  int j = ij - i*nx;
  int i1 = (ij1/nx);
  int j1 = ij1 - i1*nx;
  int i2 = (ij2/nx);
  int j2 = ij2 - i2*nx;
  double x = xpoints[j], y = ypoints[i]; // three points are needed: the point to update (x, y) and two nearby AF points
  double x1 = xpoints[j1], y1 = ypoints[i1];
  double x2 = xpoints[j2], y2 = ypoints[i2];
  double phi1 = phi[ij1], phi2 = phi[ij2]; // need all the relevent values at the two AF points
  double p1 = p[ij1], p2 = p[ij2];
  double q1 = q[ij1], q2 = q[ij2];
  double pn, qn;
  struct minfun_params mfparams = {params, x, y, x1, x2, y1, y2, p1, p2, q1, q2, phi1, phi2, &pn, &qn}; // pass to the minimization function
  double theta, a = 0., b = 1.; // initial minimization interval
  int init_status = 0;  // this is a flag used later to see if the minimum is bracketed
  bracket_interval(&init_status, &a, &b, &theta, &mfparams);
  // if the interval was found, compute theta using a minimization routine
  // if not then just use theta
  if (gsl_finite(minfun(theta, &mfparams)) == 0) { // abandon the update if minfun returns nan or inf
    return;
  }
  if (init_status == 1){ // if the minimum was bracketed, switch to the gsl minimizer
    const gsl_min_fminimizer_type *T; // for minimizing the action
    gsl_min_fminimizer *smin; 
    gsl_function F;
    F.function = &minfun; // the function to minimize
    F.params = &mfparams;
    // T = gsl_min_fminimizer_brent; // the fastest method?
    T = gsl_min_fminimizer_quad_golden;
    smin = gsl_min_fminimizer_alloc(T); // gsl minimizer memory allocation
    double epsabs = 1e-12, epsrel = 1e-15; // tolerances for the minimizer
    int status, iter = 0, max_iter = 100; // maximum number of iterations
    gsl_min_fminimizer_set(smin, &F, theta, a, b); // initialize (shit will blow up if the min was not brackted)
    gsl_set_error_handler_off();
    do { // gsl minimizer loop
      iter += 1;
      a = gsl_min_fminimizer_x_lower(smin); // continue to refine the interval
      b = gsl_min_fminimizer_x_upper(smin);
      theta = gsl_min_fminimizer_x_minimum(smin); // set theta to the computed min
      status = gsl_min_fminimizer_iterate(smin);
      status = gsl_min_test_interval(a, b, epsabs, epsrel); // checks for convergence and to make sure 0<theta<1
    } while ((status == GSL_CONTINUE) && (iter < max_iter));
    gsl_min_fminimizer_free(smin); // dealocate minimizer memory
  }
  double phi_bar = theta*phi1 + (1. - theta)*phi2; // compute all of the relevent quantities averaged with theta over both points
  double x_bar = theta*x1 + (1. - theta)*x2;
  double y_bar = theta*y1 + (1. - theta)*y2;
  double dx_bar = x - x_bar, dy_bar = y - y_bar;
  double dPhi = pn*dx_bar + qn*dy_bar; // update for phi
  double phicnew = phi_bar + dPhi; // new value of phi
  double gradphisize = sqrt(pow(pn, 2) + pow(qn, 2));
  if ((phicnew < phic[ij]) && (dPhi > 0) && (gradphisize > 0.) && (gradphisize < 5000)){ // if the tentative value of phi is less than the current considered value, update solution data
    // also check if dPhi > 0 as the action can only increase
    phic[ij] = phicnew; // if everything checks out, set the tentative value of phi to the computed value
    p[ij] = pn; // tentative values of p and q
    q[ij] = qn;
    dx[ij] = dx_bar; // update length
    dy[ij] = dy_bar; // also used as a check, should be the approximate vector field for characteristics,
  }
}
void compute_considered(int nx, int ny, 
                        double epsilon, double N, double M, double gna, double vna, double gk, double vk, double gl,
                        double vl, double Iapp, double betak, double thetak1, double thetak2, double thetana1, double thetana2,
                        double *xpoints, double *ypoints,
                        double *iadj, double *jadj, double *igamma_range, double *jgamma_range, double *ijconsidered,
                        int nadj, int ngamma_range, int nconsidered,
                        double *phic, double *phi, double *p, double *q, double *dx, double *dy,
                        bool *accepted_front){
  int i, j, ij, k, i1, j1, ij1, k1, i2, j2, ij2, k2;
  struct SlKAP_params params = {epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2};
  // update considered points within gamma_range of the new accepted front
  for (k=0;k<nconsidered;k++){ // for each considered point...
    ij = ijconsidered[k];
    i = ij/nx;
    j = ij - i*nx;
    for (k1=0;k1<ngamma_range;k1++) { // ...find adjacent AF points
      i1 = i + igamma_range[k1];
      j1 = j + jgamma_range[k1];
      ij1 = i1*nx + j1;
      if ((i1>=0) && (i1<ny) && (j1>=0) && (j1<nx) && (accepted_front[ij1] == 1)){ // make sure its within the domain,
                                                                                   // and make sure its an AF point
        for (k2=0;k2<nadj;k2++){ // find another AF point adjacent to the first for the fd formula
                                 // Loop over all suitable second AF points, and compute a tentative value of phi
                                 // Some points are computed twice, but this does not matter since this is only for initialization
          i2 = i1 + iadj[k2];
          j2 = j1 + jadj[k2];
          ij2 = i2*nx + j2;
          if ((i2>=0) && (i2<ny) && (j2>=0) && (j2<nx) && (accepted_front[ij2] == 1)){ // make sure the second point is an AF point
            upwindfd(ij, ij1, ij2, nx, ny, xpoints,  ypoints,  phic,  phi,  p,  q,  dx,  dy, &params); // fd approx
          }
        }
      }
    }
  }
  return;
}
void advance_soln(int nups, int nx, int ny, 
                  double epsilon, double N, double M, double gna, double vna, double gk, double vk, double gl, 
                  double vl, double Iapp, double betak, double thetak1, double thetak2, double thetana1, double thetana2, 
                  double *xpoints, double *ypoints,
                  double *iadj, double *jadj, double *inn, double *jnn, double *igamma_range, double *jgamma_range,
                  int nadj, int nnn, int ngamma_range,
                  double *phic, double *phi, double *p, double *q, double *dx, double *dy,
                  bool *unconsidered, bool *considered, bool *accepted_front, bool *accepted){
  //// xpoints and ypoints are linear arrays containing all the grid points
  //// iadj, jadj are the i and j indices of the adjacency stencil (i.e., (-1, 1), (0, 1), ... )
  //// inn, jnn are the i and j indices of the nearest neigbor stencil
  //// igamma_range, jgamma_range are the i and j indices for stencil for points within distance gamma
  //// nadj, nnn, ngamma_range are the number of elements in each of the arrays
  //// phic is a linear array of tentative values of phi at considered points
  //// phi is a linear array for the final solution
  //// p and q are the gradient elements
  //// dx and dy store the update direction for each computed grid point (excluding the initial data)
  //// unconsidered, considered, accepted_front, and accepted are linear boolian arrays (1 if the grid point is in the set and 0 otherwise)
  int m, i, j, ij, k, i1, j1, ij1, k1, s;
  int ijmin, imin, jmin;
  double phimin;
  size_t meshsize = nx*ny;
  gsl_vector_view phic_vecview = gsl_vector_view_array(phic, meshsize); // create a vector wrapper for the phi considered
                                                                        // to find the minimum value
  gsl_permutation *perm = gsl_permutation_alloc(meshsize); // to store the sorted indices
  gsl_vector *vphic = &phic_vecview.vector; // to store the tentative values of phi
  struct SlKAP_params params = {epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2};
  cout << "\n\n\n\n";
  cout << "%" << 0.0 << flush;
  for (m=0;m<nups;m++){
    if (m % 100 == 0){ // terminal output to see progress
      clear_line(cout, 20);
      cout << "%" << double(m*100)/double(nups) << flush;
    }
    gsl_sort_vector_index(perm, vphic);// choose the smallest considered phi value, (does this use heap sort the right way?)
    ijmin = gsl_permutation_get(perm, 0);// its linear index
    phimin = gsl_vector_get(vphic, ijmin);// min value
    imin = (ijmin/nx);
    jmin = ijmin - imin*nx;
    if (gsl_isinf(phimin)){ // if all phi considered values are GSL_POSINF, no further updates are possible
      return;
    }
    gsl_vector_set(vphic, ijmin, GSL_POSINF); // remove it from the set of considered phi values
    phi[ijmin] = phimin; // set the final solution point
    considered[ijmin] = 0; // remove it from considered...
    accepted_front[ijmin] = 1; //... to accepted front
    // move unconsidered points adjacent to the new accepted front point to considered
    for (k=0;k<nadj;k++){
      i = imin + iadj[k];
      j = jmin + jadj[k];
      ij = i*nx + j;
      if ((i>=0) && (i<ny) && (j>=0) && (j<nx) && (unconsidered[ij] == 1)){
        unconsidered[ij] = 0;
        considered[ij] = 1;
      }
    }
    // quick check to see if the new accepted front point should be immediately moved to accepted
    s = 1;
    for (k=0;k<nnn;k++){ // for all nearest neighbors of the new accepted front point
      i = imin + inn[k];
      j = jmin + jnn[k];
      ij = i*nx + j;
      if ((i>=0) && (i<ny) && (j>=0) && (j<nx) && (considered[ij] == 1)){ // check if the new accepted_front
                                                                          // point has a nearest neigbor considered point
        s = 0;
      }
    }
    if (s == 1){ // if it does not then move it to accepted
      accepted_front[ijmin] = 0;
      accepted[ijmin] = 1;
    }
    // update the accepted front points, moving accepted front points with
    // no nearest neighbor considered points to accepted
    for (k=0;k<nadj;k++){ // for all points adjacent to the new accepted front point
      i = imin + iadj[k];
      j = jmin + jadj[k];
      ij = i*nx + j;
      if ((i>=0) && (i<ny) && (j>=0) && (j<nx) && (accepted_front[ij] == 1)){ // check to see if it is also an
                                                                              // accepted front point
        s = 1;
        for (k1=0;k1<nnn;k1++){ // for all nearest neighbors of the adjacent accepted front point
          i1 = i + inn[k1];
          j1 = j + jnn[k1];
          ij1 = i1*nx + j1;
          if ((i1>=0) && (i1<ny) && (j1>=0) && (j1<nx) && (considered[ij1] == 1)){ // check if the adjacent
                                                                                   // accepted_front point has
                                                                                   // a nearest neigbor considered point
            s = 0;  // if it does then s = 0 and the point wont be moved
          }
        }
        if (s == 1){ // if it does not then s = 1 and the point is moved to accepted
          accepted_front[ij] = 0;
          accepted[ij] = 1;
        }
      }
    }
    // update considered points within gamma_range of the new accepted front point
    for (k=0;k<ngamma_range;k++){ // loop over each point in gamma_range of the new AF point...
      i = imin + igamma_range[k];
      j = jmin + jgamma_range[k];
      ij = i*nx + j;
      if ((i>=0) && (i<ny) && (j>=0) && (j<nx) && (considered[ij] == 1)){ // ...check to see if it is a considered point. If it is, update
        for (k1=0;k1<nadj;k1++){ // Form pairs from the new AF point and one of its adjacent AF points
          i1 = imin + iadj[k1];
          j1 = jmin + jadj[k1];
          ij1 = i1*nx + j1;
          if ((i1>=0) && (i1<ny) && (j1>=0) && (j1<nx) && (accepted_front[ij1] == 1)){ // make sure it is also an AF point
            upwindfd(ij, ij1, ijmin, nx, ny, xpoints,  ypoints,  phic,  phi,  p,  q,  dx,  dy, &params); // compute the fd
          }
        }
      }
    }
  }
  gsl_permutation_free(perm); // free memory
  return;
}
