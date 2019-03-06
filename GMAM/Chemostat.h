// Written by Jay Newby 2019
// Send questions and comments to jnewby@ualberta.ca
#ifndef CSTAT_H
#define CSTAT_H

#ifdef __cplusplus
extern "C" {
#endif

struct model_params {
  // Model parameters
  double alpha1;
  double alpha2;
  double ey;
};

struct model_params getParams(double *pars);
double H(double x, double y, double p, double q, void *params);
double Hp(double x, double y, double p, double q, void *params);
double Hq(double x, double y, double p, double q, void *params);
double Hx(double x, double y, double p, double q, void *params);
double Hy(double x, double y, double p, double q, void *params);
double Hpp(double x, double y, double p, double q, void *params);
double Hqq(double x, double y, double p, double q, void *params);
double Hpq(double x, double y, double p, double q, void *params);
double Hpx(double x, double y, double p, double q, void *params);
double Hpy(double x, double y, double p, double q, void *params);
double Hqx(double x, double y, double p, double q, void *params);
double Hqy(double x, double y, double p, double q, void *params);

#ifdef __cplusplus
}
#endif

#endif
