// Written by Jay Newby 2019
// Send questions and comments to jnewby@ualberta.ca

#include "math.h"
#include "Chemostat.h"
//using namespace std;


struct model_params getParams(double pars[]) {
  struct model_params params = {
    pars[0],
    pars[1],
    pars[2]
  }; // parameter struct
  return params;
}

////////////////////////////////////////////////////////////////
// Various partial derivatives of the Hamiltonian, H(x, y, p, q).
////////////////////////////////////////////////////////////////
double H(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double h =
    (exp(p) - 1)*pars.alpha1*(y/(1+y))*x + (exp(-p) - 1)*x
    + q*(-(y/(1+y))*x - y + pars.alpha2)
    + q*q/2*pars.ey*((y/(1+y))*x + y + pars.alpha2);
  return h;
}
double Hp(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hp = exp(p)*pars.alpha1*(y/(1+y))*x - exp(-p)*x;
  return hp;
}
double Hq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hq = -(y/(1+y))*x - y + pars.alpha2
              + q*pars.ey*((y/(1+y))*x + y + pars.alpha2);
  return hq;
}
double Hx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hx = (exp(p) - 1)*pars.alpha1*(y/(1+y)) + (exp(-p) - 1)
              + q*(-(y/(1+y)))
              + q*q/2*pars.ey*((y/(1+y)));
  return hx;
}
double Hy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hy = (exp(p) - 1)*pars.alpha1*(1/(1+y)/(1+y))*x
              + q*(-(1/(1+y)/(1+y))*x - 1.0)
              + q*q/2*pars.ey*((1/(1+y)/(1+y))*x + 1);
  return hy;
}
double Hpp(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpp = exp(p)*pars.alpha1*(y/(1+y))*x + exp(-p)*x;
  return hpp;
}
double Hqq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  return pars.ey*((y/(1+y))*x + y + pars.alpha2);
}
double Hpq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  return 0.0;
}
double Hpx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpx = exp(p)*pars.alpha1*(y/(1+y)) - exp(-p);
  return hpx;
}
double Hpy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpy = exp(p)*pars.alpha1*(1/(1+y)/(1+y))*x;
  return hpy;
}
double Hqx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hqx = -(y/(1+y))
              + q*pars.ey*((y/(1+y)));
  return hqx;
}
double Hqy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hqy = -(1/(1+y)/(1+y))*x - 1
              + q*pars.ey*((1/(1+y)/(1+y))*x + 1);
  return hqy;
}
