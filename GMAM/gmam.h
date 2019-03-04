#ifndef GMAM_H
#define GMAM_H

#ifdef __cplusplus
extern "C" {
#endif
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

void action_solve(double x, double y, double dx, double dy, double *ppnt,
    double *qpnt, void *pars);
void gmam(double *t, double *s, double *x, double *y, double *p,
    double *q, double *phi, int Npnts, int Niter, double dt, double abs_err,
    double *error, void *pars);

#ifdef __cplusplus
}
#endif

#endif
