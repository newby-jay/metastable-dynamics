// Written by Jay Newby 2019
// Send questions and comments to jnewby@ualberta.ca

#include "math.h"
#include "Lotka_Volterra.h"
//using namespace std;


struct model_params getParams(double pars[]) {
  struct model_params params = {
    pars[0],
    pars[1],
    pars[2],
    pars[3],
    pars[4]
  }; // parameter struct
  return params;
}

////////////////////////////////////////////////////////////////
// Various partial derivatives of the Hamiltonian, H(x, y, p, q).
////////////////////////////////////////////////////////////////
double H(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double h =
    (exp(p) - 1)*pars.alpha*x + (exp(-p) - 1)*(pars.beta*x*y + pars.rho*x*x)
    + (exp(q) - 1)*pars.delta*x*y + (exp(-q) - 1)*pars.gamma*y;
  return h;
}
double Hp(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hp = exp(p)*pars.alpha*x - exp(-p)*(pars.beta*x*y + pars.rho*x*x);
  return hp;
}
double Hq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hq = exp(q)*pars.delta*x*y - exp(-q)*pars.gamma*y;
  return hq;
}
double Hx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hx =
    (exp(p) - 1)*pars.alpha + (exp(-p) - 1)*(pars.beta*y + 2*pars.rho*x)
    + (exp(q) - 1)*pars.delta*y;
  return hx;
}
double Hy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hy =
    (exp(-p) - 1)*pars.beta*x
    + (exp(q) - 1)*pars.delta*x + (exp(-q) - 1)*pars.gamma;
  return hy;
}
double Hpp(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpp = exp(p)*pars.alpha*x + exp(-p)*(pars.beta*x*y + pars.rho*x*x);
  return hpp;
}
double Hqq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hqq = exp(q)*pars.delta*x*y + exp(-q)*pars.gamma*y;
  return hqq;
}
double Hpq(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  return 0.0;
}
double Hpx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpx = exp(p)*pars.alpha - exp(-p)*(pars.beta*y + 2*pars.rho*x);
  return hpx;
}
double Hpy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hpy = -exp(-p)*pars.beta*x;
  return hpy;
}
double Hqx(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hqx = exp(q)*pars.delta*y;
  return hqx;
}
double Hqy(double x, double y, double p, double q, void *params) {
  struct model_params pars = *(struct model_params *)params;
  double hqy = exp(q)*pars.delta*x - exp(-q)*pars.gamma;
  return hqy;
}
