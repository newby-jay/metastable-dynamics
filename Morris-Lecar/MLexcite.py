# Written by Jay Newby 2013
# Send questions and comments to newby@math.utah.edu.

from pylab import *
import time
from scipy.optimize import fsolve
import weave
from string import Template

class MLexcite:
    def __init__(self, args):
        """N, sodium channels.
        M, potassium channels.
        gna, sodium conductance.
        vna, sodium reversal potential.
        gk, postassium conductance.
        vk, potassium reversal potential.
        gl, leak conductance.
        vl, leak reversal potential.
        Iapp, applied current.
        betak, potassium channel rate constant.
        thetak, potassium channel activation constant.
        thetana1, sodium channel activation constant.
        thetana2, sodium channel activation constant.
        """
        for k in args:
            setattr(self, k, args[k])
        self.hgam = 1./(self.epsilon*self.M)
        self.N = args['N']/args['taum']
        self.parameters = array([self.__dict__[k] for k in ['epsilon', 'N', 'M', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2']])
        self.x0, self.y0 = self._find_fixedpoint(self.v00)
        self._create_functions()
    def _find_fixedpoint(self, v0):
        "finds fixed points given initial guess v0"
        def fpfun(v):
            minf = (1. + tanh(2.*(self.thetana1*v + self.thetana2)))/2.
            alpha = exp(-(self.thetak1*v + self.thetak2))
            beta = 1./alpha
            winf = alpha/(alpha + beta)
            return winf + (self.gna*minf*(self.vna - v) + self.gl*(self.vl - v) + self.Iapp)/(self.gk*(self.vk - v))
        x0 = fsolve(fpfun, v0, xtol=1e-16, maxfev = 1000)
        if not isscalar(x0):
            x0 = x0[0]
        y0 = exp(-(self.thetak1*x0 + self.thetak2))/(exp(-(self.thetak1*x0 + self.thetak2)) + exp(self.thetak1*x0 + self.thetak2))
        return x0, y0
    def find_lc(self, x0, y0, tinit, tperiod, Npts):
        """ find stable limit cycle
            x0, y0: the initial condition to use for numerically computing the limit cycle.
            tinit: time to allow initial convergence to limit cycle from (x0, y0) to xc.
            tperiod: a guess for the period of the lc.
            Npts: the number of discrete points on the lc."""
        template_file = open('gslode/gslode_template_newham.c')
        code_template = Template(template_file.read())
        scode = code_template.substitute(self.func_string)
        def deterministic_traj(x0, y0, tmax, Ncols):
            "use GSL ode integrator to compute deterministic trajectories"
            code = '//'+"""
            //deterministic V0.00
            double *solsub = sol; // For some reason, I can't pass sol to gsl_odesolver()
            gsl_odesolver_det(int(method), Ncols, x0, y0, epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2, tmax, dt, rel_err, solsub);
            """
            dt, rel_err, method = tmax/Ncols, 1e-12, 10
            epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2 = self.parameters
            sol = zeros((4, Ncols))*nan
            weave.inline(code,
               ['method', 'Ncols', 'sol', 'x0', 'y0', 'epsilon', 'N', 'M', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'tmax', 'dt', 'rel_err'],
               support_code = scode,
               libraries=['gsl'],
               include_dirs=['/opt/local/include'],
               library_dirs=['/opt/local/lib'])
            traj = dict((k, sol[j, :]) for k, j in zip(['t', 'x', 'y', 's'], arange(4)))
            return traj
        xcinit = deterministic_traj(x0, y0, tinit, Npts)
        xc0, yc0 = xcinit['x'][-1], xcinit['y'][-1]
        xc = deterministic_traj(xc0, yc0, tperiod, Npts)
        NT, T, L = nan, nan, nan
        for n in arange(4, Npts):
            if (((xc0 - xc['x'][n])*(xc0 - xc['x'][n-1]) < 0) and ((yc0 - xc['y'][n])*(yc0 - xc['y'][n-1]) < 0)):
                NT = n
                T = xc['t'][n]
                L = xc['s'][n]
                break
        xc = dict((k, xc[k][:NT+1]) for k in ['t', 'x', 'y', 's']) #truncated trajectory
        self.xc = xc
    def _create_functions(self):
        "Store the many partial derivatives of H as strings to copy into cfiles.  I do it this way so that if I need to change something, I only have to change it once."
        abcs = {'b' : '(-1./(1. - minf) + (2.*g + fna)/N*p)',
                'dbdx': '(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)',
                'dbdy': '2.*gw*p/N',
                'dbdp': '(2.*g + fna)/N',
                'c': '(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N)',
                'dcdx': '(-(dminf*fna + minf*dfna + gv)/(1. - minf)*p - (minf*fna + g)*dminf/pow((1. - minf), 2)*p + (dfna + gv)*g*pow(p, 2)/N + (fna + g)*gv*pow(p, 2)/N)',
                'dcdy': '(-gw/(1. - minf)*p + (gw*g + (fna + g)*gw)/N*pow(p, 2))',
                'dcdp': '(-(minf*fna + g)/(1. - minf) + 2.*(fna + g)*g/N*p)',
                'd2bdp': '0.',
                'd2cdp': '2.*(fna + g)*g/N',
                'ddbdxdx': '(-ddminf/pow((1. - minf), 2) - 2.*pow(dminf, 2)/pow((1. - minf), 3))',
                'ddbdxdy': '2.*gvw*p/N',
                'ddbdydy': '0.',
                'ddcdxdx': '(-((ddminf*fna + dminf*dfna + dminf*dfna)/(1. - minf) + (dminf*fna + minf*dfna + gv)*dminf/pow((1. - minf), 2))*p - (((dminf*fna + minf*dfna + gv)*dminf + (minf*fna + g)*ddminf)/pow((1. - minf), 2) + 2.*pow(dminf, 2)*(minf*fna + g)/pow((1. - minf), 3))*p + 2.*(dfna + gv)*gv*pow(p, 2)/N)',
                'ddcdxdy': '(-(gvw/(1. - minf) + gw*dminf/pow((1. - minf), 2))*p + (gvw*g + gw*gv + (dfna + gv)*gw + (fna + g)*gvw)/N*pow(p, 2))',
                'ddcdydy': '(2*pow(gw, 2)/N*pow(p, 2))',
                'ddbdpdx': '(2.*gv + dfna)/N',
                'ddcdpdx': '(-(dminf*fna + minf*dfna + gv)/(1. - minf) - dminf*(minf*fna + g)/pow((1. - minf), 2) + 2./N*p*((dfna + gv)*g + (fna + g)*gv))',
                'ddbdpdy': '2.*gw/N',
                'ddcdpdy': '(-gw/(1. - minf) + 2./N*p*(gw*g + (fna + g)*gw))',
                'sqrttrm': 'sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4./N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N))'}
        H = Template('(h + N/2*($b + $sqrttrm))').substitute(abcs)
        Hx = Template('(hv + N/2*$dbdx + N/4*(2.*$dbdx*$b - 4/N*$dcdx)/$sqrttrm)').substitute(abcs)
        Hy = Template('(hw + N/2*$dbdy + N/4*(2.*$dbdy*$b - 4/N*$dcdy)/$sqrttrm)').substitute(abcs)
        Hp =      Template('(N/2*$dbdp + N/4*(2.*$dbdp*$b - 4/N*$dcdp)/$sqrttrm)').substitute(abcs)
        Hpp =                 Template('(N/4.*(2.*pow($dbdp, 2) - 4/N*$d2cdp)/$sqrttrm - N/8.*pow((2.*$dbdp*$b - 4/N*$dcdp), 2)/pow($sqrttrm, 3))').substitute(abcs)
        Hq = '(hq)'
        Hpq = '0.'
        Hqq = '(hqq)'
        Hpx = Template('(N/2*$ddbdpdx + N/4*((2.*$ddbdpdx*$b + 2.*$dbdp*$dbdx - 4/N*$ddcdpdx)/$sqrttrm - 0.5*(2.*$dbdp*$b - 4/N*$dcdp)*(2.*$dbdx*$b - 4/N*$dcdx)/pow($sqrttrm, 3)))').substitute(abcs)
        Hpy = Template('(N/2*$ddbdpdy + N/4*((2.*$ddbdpdy*$b + 2.*$dbdp*$dbdy - 4/N*$ddcdpdy)/$sqrttrm - 0.5*(2.*$dbdp*$b - 4/N*$dcdp)*(2.*$dbdy*$b - 4/N*$dcdy)/pow($sqrttrm, 3)))').substitute(abcs)
        Hqx = 'hqv'
        Hqy = 'hqw'
        Hxx = Template('(hvv + N/2*$ddbdxdx + N/4*((2.*$ddbdxdx*$b + 2.*$dbdx*$dbdx - 4/N*$ddcdxdx) - 0.5*pow((2.*$dbdx*$b - 4/N*$dcdx), 2)/pow($sqrttrm, 2))/$sqrttrm)').substitute(abcs)
        Hxy = Template('(hvw + N/2*$ddbdxdy + N/4*((2.*$ddbdxdy*$b + 2.*$dbdx*$dbdy - 4/N*$ddcdxdy)/$sqrttrm - 0.5*(2.*$dbdx*$b - 4/N*$dcdx)*(2.*$dbdy*$b - 4/N*$dcdy)/pow($sqrttrm, 3)))').substitute(abcs)
        Hyy = '0.'
        self.func_string = {'H': H, # Hamiltonian and its derivatives
                            'Hp': Hp, 'Hq': Hq, 'Hx': Hx, 'Hy': Hy,
                            'Hpp': Hpp, 'Hpq': Hpq, 'Hqq': Hqq,
                            'Hpx': Hpx, 'Hpy': Hpy, 'Hqx': Hqx, 'Hqy': Hqy,
                            'Hxx': Hxx, 'Hxy': Hxy, 'Hyy': Hyy,
                            'J00': Hpx, 'J01': Hpy, 'J02': Hpp, 'J03': Hpq, # Jacobian of hamiltonian system
                            'J10': Hqx, 'J11': Hqy, 'J12': Hpq, 'J13': Hqq,
                            'J20': str('-(%s)' % Hxx), 'J21': str('-(%s)' % Hxy), 'J22': str('-(%s)' % Hpx), 'J23': str('-(%s)' % Hqx),
                            'J30': str('-(%s)' % Hxy), 'J31': str('-(%s)' % Hyy), 'J32': str('-(%s)' % Hpy), 'J33': str('-(%s)' % Hqy),
                            'minf': '(1. + tfev)/2.', # aux functions
                            'tfev': 'tanh(2.*(thetana1*v + thetana2))',
                            'dminf': '(1. - pow(tfev, 2))*thetana1',
                            'ddminf': '-2.*tfev*(1. - pow(tfev, 2))*2.*pow(thetana1, 2)',
                            'fna': 'gna*(vna - v)',
                            'dfna': '-gna',
                            'g': 'w*gk*(vk - v) + gl*(vl - v) + Iapp',
                            'gv': '-(w*gk + gl)',
                            'gw': 'gk*(vk - v)',
                            'gvw': '-gk',
                            'h': 'betak*(w*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1)/hgam + (1-w)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1)/hgam)',
                            'hv': 'betak*(w*thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1)/hgam - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1)/hgam)',
                            'hw': 'betak*(exp((thetak1*v + thetak2))*(exp(-hgam*q)-1)/hgam - exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1)/hgam)',
                            'hq': 'betak*(-w*exp((thetak1*v + thetak2))*exp(-hgam*q) + (1-w)*exp(-(thetak1*v + thetak2))*exp(hgam*q))',
                            'hqq':'betak*hgam*(w*exp((thetak1*v + thetak2))*exp(-hgam*q) + (1-w)*exp(-(thetak1*v + thetak2))*exp(hgam*q))',
                            'hvv':'betak*(w*pow(thetak1, 2)*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1)/hgam + (1-w)*pow(thetak1, 2)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1)/hgam)',
                            'hvw': 'betak*(thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1)/hgam + thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1)/hgam)',
                            'hww': '0.',
                            'hqv': 'betak*(-w*thetak1*exp((thetak1*v + thetak2))*exp(-hgam*q) - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*exp(hgam*q))',
                            'hqw': 'betak*(-exp((thetak1*v + thetak2))*exp(-hgam*q) - exp(-(thetak1*v + thetak2))*exp(hgam*q))'}
    def g(self, v, w):
        return w*self.gk*(self.vk - v) + self.gl*(self.vl - v) + self.Iapp
    def fna(self, v):
        return self.gna*(self.vna - v)
    def minf(self, v):
        return (1. + tanh(2.*(self.thetana1*v + self.thetana2)))/2.
    def winf(self, v):
        return exp(-(self.thetak1*v + self.thetak2))/(exp(-(self.thetak1*v + self.thetak2)) + exp(self.thetak1*v + self.thetak2))
    def vinf(self, v):
        return -(self.gna*self.minf(v)*(self.vna - v) + self.gl*(self.vl - v) + self.Iapp)/(self.gk*(self.vk - v))
    def detjacxx(self, v, w):
        "xx element of the deterministic Jacobian"
        minf = (1. + tanh(2.*(self.thetana1*v + self.thetana2)))/2.
        dminf = (1. - tanh(2.*(self.thetana1*v + self.thetana2))**2)*self.thetana1
        gv = w*(-self.gk) - self.gl
        fna = self.gna*(self.vna - v)
        dfna = -self.gna
        return dminf*fna + minf*dfna + gv
    def detjacxy(self, v, w):
        "xy element of the deterministic Jacobian"
        return self.gk*(self.vk - v)
    def detjacyx(self, v, w):
        "yx element of the deterministic Jacobian"
        alpha = exp(-(self.thetak1*v + self.thetak2))
        return self.betak*(-w*self.thetak1/alpha - (1. - w)*self.thetak1*alpha)
    def detjacyy(self, v, w):
        "yy element of the deterministic Jacobian"
        alpha = exp(-(self.thetak1*v + self.thetak2))
        return -self.betak*(1./alpha + alpha)
    def h(self, v, w, p, q):
        return self.betak/self.hgam*(w*exp((self.thetak1*v + self.thetak2))*(exp(-self.hgam*q) - 1.) + (1. - w)*exp(-(self.thetak1*v + self.thetak2))*(exp(self.hgam*q) - 1.))
    def Z(self, v, w):
        "Hessian matrix at a fixed point.  Solution of the Ricatti equation."
        minf = (1. + tanh(2.*(self.thetana1*v + self.thetana2)))/2.
        dminf = (1. - tanh(2.*(self.thetana1*v + self.thetana2))**2)*self.thetana1
        g = w*self.gk*(self.vk - v) + self.gl*(self.vl - v) + self.Iapp
        gv = w*(-self.gk) - self.gl
        gw = self.gk*(self.vk - v)
        fna = self.gna*(self.vna - v)
        dfna = -self.gna
        h = 0.
        alpha = exp(-(self.thetak1*v + self.thetak2))
        hq = self.betak*(-w/alpha + (1. - w)*alpha)
        hqq = self.betak*self.hgam*(w/alpha + (1. - w)*alpha)
        hx = 0.
        hy = 0.
        hqx = self.betak*(-w*self.thetak1/alpha - (1. - w)*self.thetak1*alpha)
        hqy = self.betak*(-1./alpha - alpha)
        Hpp = -2.*(1. - minf)*g*(fna + g)/self.N
        Hpq = -(1. - minf)*(2.*g + fna)*hq/self.N
        Hqq = -(-hqq*self.N + 2.*hq**2 - 2.*hq**2*minf)/self.N
        Hpx = (dminf*fna + minf*dfna + gv)
        Hpy = gw
        Hqx = hqx
        Hqy = hqy
        Det = 2.*(Hpx + Hqy)*(Hpx*Hqy - Hqx*Hpy)
        Q11 = (Hqx*Hpy - Hqy*(Hpx + Hqy))*Hpp + 2.*Hpy*Hqy*Hpq - Hqq*Hpy**2
        Q12 = (-2.*Hpx*Hpq + Hqx*Hpp)*Hqy + Hqq*Hpx*Hpy
        Q22 = -Hqx**2*Hpp + (Hqq*Hpy + 2.*Hpx*Hpq)*Hqx - Hpx*Hqq*(Hpx + Hqy)
        Z = matrix([[Q22, -Q12], [-Q12, Q11]])*Det/(Q11*Q22 - Q12**2)
        return Z
    def Zr1(self, v, w):
        "Hessian matrix at a fixed point.  Rank 1 solutions of the Ricatti equation, assuming that Hpq = 0."
        minf = (1. + tanh(2.*(self.thetana1*v + self.thetana2)))/2.
        dminf = (1. - tanh(2.*(self.thetana1*v + self.thetana2))**2)*self.thetana1
        g = w*self.gk*(self.vk - v) + self.gl*(self.vl - v) + self.Iapp
        gv = w*(-self.gk) - self.gl
        gw = self.gk*(self.vk - v)
        fna = self.gna*(self.vna - v)
        dfna = -self.gna
        h = 0.
        alpha = exp(-(self.thetak1*v + self.thetak2))
        hq = self.betak*(-w/alpha + (1. - w)*alpha)
        hqq = self.betak*self.hgam*(w/alpha + (1. - w)*alpha)
        hx = 0.
        hy = 0.
        hqx = self.betak*(-w*self.thetak1/alpha - (1. - w)*self.thetak1*alpha)
        hqy = self.betak*(-1./alpha - alpha)
        Hpp = (self.N/2*(2.*g + fna)/self.N + self.N/4*(2.*(2.*g + fna)/self.N*(-1./(1. - minf)) - 4/self.N*(-(minf*fna + g)/(1. - minf)))/sqrt(pow(-1./(1. - minf), 2)))
        Hqq = hqq
        Hpx = self.detjacxx(v, w)
        Hpy = self.detjacxy(v, w)
        Hqx = hqx
        Hqy = hqy
        a = Hqx**2*Hpp**2+Hqq*(Hpx**2-2.*Hpx*Hqy+2.*Hpy*Hqx+Hqy**2)*Hpp+Hpy**2*Hqq**2
        b = (2.*Hqy**3-4.*Hpx*Hqy**2+(2.*Hpx**2+6.*Hpy*Hqx)*Hqy-2.*Hpx*Hpy*Hqx)*Hpp+2*Hpy**2*Hqq*(Hpx+Hqy)
        c = 4.*Hpx*Hpy**2*Hqy-4.*Hpy**3*Hqx
        R1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
        R2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)
        Z1xx = -((Hpp*Hqx**2+Hpx**2*Hqq-Hpx*Hqq*Hqy+Hpy*Hqq*Hqx)*R1 + 2*Hpx**2*Hqy-2*Hpx*Hpy*Hqx-2*Hpx*Hqy**2+2*Hpy*Hqx*Hqy)/(Hpp*Hpx*Hqy-Hpp*Hpy*Hqx-Hpp*Hqy**2-Hpy**2*Hqq)
        Z1xy = -((Hpp*Hqx*Hqy+Hpx*Hpy*Hqq)*R1 + 2*Hpx*Hpy*Hqy-2*Hpy**2*Hqx)/(Hpp*Hpx*Hqy-Hpp*Hpy*Hqx-Hpp*Hqy**2-Hpy**2*Hqq)
        Z1yy = R1
        Z2xx = -((Hpp*Hqx**2+Hpx**2*Hqq-Hpx*Hqq*Hqy+Hpy*Hqq*Hqx)*R2 + 2*Hpx**2*Hqy-2*Hpx*Hpy*Hqx-2*Hpx*Hqy**2+2*Hpy*Hqx*Hqy)/(Hpp*Hpx*Hqy-Hpp*Hpy*Hqx-Hpp*Hqy**2-Hpy**2*Hqq)
        Z2xy = -((Hpp*Hqx*Hqy+Hpx*Hpy*Hqq)*R2 + 2*Hpx*Hpy*Hqy-2*Hpy**2*Hqx)/(Hpp*Hpx*Hqy-Hpp*Hpy*Hqx-Hpp*Hqy**2-Hpy**2*Hqq)
        Z2yy = R2
        return matrix([[Z1xx, Z1xy], [Z1xy, Z1yy]]), matrix([[Z2xx, Z2xy], [Z2xy, Z2yy]])
    def H(self, v, w, p, q):
        "Hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        fna = self.fna(v)
        g = self.g(v, w)
        h = betak/hgam*(w*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) + (1-w)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        return (h - N/2.*(-(-1./(1. - minf) + (2.*g + fna)/N*p) - sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4./N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N))))
    def Hx(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        dminf = (1. - tanh(2.*(self.thetana1*v + self.thetana2))**2)*self.thetana1
        fna = self.fna(v)
        dfna = -self.gna
        g = self.g(v, w)
        gv  = -(w*self.gk + self.gl)
        hv = betak/hgam*(w*thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        return ((hv + N/2*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p) + N/4*(2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4/N*(-(dminf*fna + minf*dfna + gv)/(1. - minf)*p - (minf*fna + g)*dminf/pow((1. - minf), 2)*p + (dfna + gv)*g*pow(p, 2)/N + (fna + g)*gv*pow(p, 2)/N))/sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N))))
    def Hy(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        fna = self.fna(v)
        g = self.g(v, w)
        gw  = self.gk*(self.vk - v)
        hw = betak/hgam*(exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) - exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        return ((hw + N/2*2.*gw*p/N + N/4*(2.*2.*gw*p/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4/N*(-gw/(1. - minf)*p + (gw*g + (fna + g)*gw)/N*pow(p, 2)))/sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N))))
    def Hp(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        fna = self.fna(v)
        g = self.g(v, w)
        return (N/2*(2.*g + fna)/N + N/4*(2.*(2.*g + fna)/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4/N*(-(minf*fna + g)/(1. - minf) + 2.*(fna + g)*g/N*p))/sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N)))
    def Hq(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        return betak*(-w*exp((thetak1*v + thetak2))*exp(-hgam*q) + (1-w)*exp(-(thetak1*v + thetak2))*exp(hgam*q))
    def Hpp(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        fna = self.fna(v)
        g = self.g(v, w)
        return (N/4.*(2.*pow((2.*g + fna)/N, 2.) - 4./N*2.*(fna + g)*g/N)/sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2.) - 4./N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2.)/N))
                - N/8.*pow((2.*(2.*g + fna)/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(minf*fna + g)/(1. - minf) + 2.*(fna + g)*g/N*p)), 2.)/pow(sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2.) - 4./N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N)), 3.))
    def Hqq(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        return betak*hgam*(w*exp((thetak1*v + thetak2))*exp(-hgam*q) + (1-w)*exp(-(thetak1*v + thetak2))*exp(hgam*q))
    def Hpq(self, v, w, p, q):
        "partial derivative of the hamiltonian"
        return 0.
    def Hpx(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam; N = self.N
        minf = self.minf(v)
        dminf = (1. - tanh(2.*(self.thetana1*v + self.thetana2))**2)*self.thetana1
        fna = self.fna(v)
        dfna = -self.gna
        g = self.g(v, w)
        gv = w*(-self.gk) - self.gl
        sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N));
        return (N/2.*(2.*gv + dfna)/N + N/4.*((2.*(2.*gv + dfna)/N*(-1./(1. - minf) + (2.*g + fna)/N*p) + 2.*(2.*g + fna)/N*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p) - 4./N*(-(dminf*fna + minf*dfna + gv)/(1. - minf) - dminf*(minf*fna + g)/pow((1. - minf), 2) + 2./N*p*((dfna + gv)*g + (fna + g)*gv)))/sqrttrm - 0.5*(2.*(2.*g + fna)/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(minf*fna + g)/(1. - minf) + 2.*(fna + g)*g/N*p))*(2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(dminf*fna + minf*dfna + gv)/(1. - minf)*p - (minf*fna + g)*dminf/pow((1. - minf), 2)*p + (dfna + gv)*g*pow(p, 2)/N + (fna + g)*gv*pow(p, 2)/N))/pow(sqrttrm, 3)))
    def Hpy(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam; N = self.N
        minf = self.minf(v)
        fna = self.fna(v)
        g = self.g(v, w)
        gw = self.gk*(self.vk - v)
        sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N));
        return (N/2.*2.*gw/N + N/4.*((2.*2.*gw/N*(-1./(1. - minf) + (2.*g + fna)/N*p) + 2.*(2.*g + fna)/N*2.*gw*p/N - 4./N*(-gw/(1. - minf) + 2./N*p*(gw*g + (fna + g)*gw)))/sqrttrm - 0.5*(2.*(2.*g + fna)/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(minf*fna + g)/(1. - minf) + 2.*(fna + g)*g/N*p))*(2.*2.*gw*p/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-gw/(1. - minf)*p + (gw*g + (fna + g)*gw)/N*pow(p, 2)))/pow(sqrttrm, 3)))
    def Hqx(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam;
        thetak1 = self.thetak1; thetak2 = self.thetak2;
        return betak*(-w*thetak1*exp((thetak1*v + thetak2))*exp(-hgam*q) - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*exp(hgam*q))
    def Hqy(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam;
        thetak1 = self.thetak1; thetak2 = self.thetak2;
        return betak*(-exp((thetak1*v + thetak2))*exp(-hgam*q) - exp(-(thetak1*v + thetak2))*exp(hgam*q))
    def Hxx(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        dminf = (1. - tanh(2.*(thetana1*v + thetana2))**2)*thetana1
        ddminf = -2.*tanh(2.*(thetana1*v + thetana2))*(1. - pow(tanh(2.*(thetana1*v + thetana2)), 2))*2.*pow(thetana1, 2)
        fna = self.fna(v)
        dfna = -self.gna
        g = self.g(v, w)
        gv = w*(-self.gk) - self.gl
        sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N));
        h = betak/hgam*(w*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) + (1-w)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        hv = betak/hgam*(w*thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        hvv = betak/hgam*(w*pow(thetak1, 2)*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) + (1-w)*pow(thetak1, 2)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        return hvv + N/2.*(-ddminf/pow((1. - minf), 2) - 2.*pow(dminf, 2)/pow((1. - minf), 3)) + N/4.*((2.*(-ddminf/pow((1. - minf), 2) - 2.*pow(dminf, 2)/pow((1. - minf), 3))*(-1./(1. - minf) + (2.*g + fna)/N*p) + 2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p) - 4./N*(-((ddminf*fna + dminf*dfna + dminf*dfna)/(1. - minf) + (dminf*fna + minf*dfna + gv)*dminf/pow((1. - minf), 2))*p - (((dminf*fna + minf*dfna + gv)*dminf + (minf*fna + g)*ddminf)/pow((1. - minf), 2) + 2.*pow(dminf, 2)*(minf*fna + g)/pow((1. - minf), 3))*p + 2.*(dfna + gv)*gv*pow(p, 2)/N)) - 0.5*pow((2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(dminf*fna + minf*dfna + gv)/(1. - minf)*p - (minf*fna + g)*dminf/pow((1. - minf), 2)*p + (dfna + gv)*g*pow(p, 2)/N + (fna + g)*gv*pow(p, 2)/N)), 2)/pow(sqrttrm, 2))/sqrttrm
    def Hxy(self, v, w, p, q):
        betak = self.betak; hgam = self.hgam; N = self.N
        thetana1 = self.thetana1; thetana2 = self.thetana2; thetak1 = self.thetak1; thetak2 = self.thetak2;
        minf = self.minf(v)
        dminf = (1. - tanh(2.*(thetana1*v + thetana2))**2)*thetana1
        fna = self.fna(v)
        dfna = -self.gna
        g = self.g(v, w)
        gv = w*(-self.gk) - self.gl
        gw = self.gk*(self.vk - v)
        gvw = -self.gk
        sqrttrm = sqrt(pow((-1./(1. - minf) + (2.*g + fna)/N*p), 2) - 4/N*(-(minf*fna + g)/(1. - minf)*p + (fna + g)*g*pow(p, 2)/N));
        h = betak/hgam*(w*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) + (1-w)*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        hv = betak/hgam*(w*thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) - (1-w)*thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        hw = betak/hgam*(exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) - exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        hvw = betak/hgam*(thetak1*exp((thetak1*v + thetak2))*(exp(-hgam*q)-1) + thetak1*exp(-(thetak1*v + thetak2))*(exp(hgam*q)-1))
        return (hvw + N/2.*2.*gvw*p/N + N/4.*((2.*2.*gvw*p/N*(-1./(1. - minf) + (2.*g + fna)/N*p) + 2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*2.*gw*p/N - 4./N*(-(gvw/(1. - minf) + gw*dminf/pow((1. - minf), 2))*p + (gvw*g + gw*gv + (dfna + gv)*gw + (fna + g)*gvw)/N*pow(p, 2)))/sqrttrm - 0.5*(2.*(-dminf/pow((1. - minf), 2) + (2.*gv + dfna)/N*p)*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-(dminf*fna + minf*dfna + gv)/(1. - minf)*p - (minf*fna + g)*dminf/pow((1. - minf), 2)*p + (dfna + gv)*g*pow(p, 2)/N + (fna + g)*gv*pow(p, 2)/N))*(2.*2.*gw*p/N*(-1./(1. - minf) + (2.*g + fna)/N*p) - 4./N*(-gw/(1. - minf)*p + (gw*g + (fna + g)*gw)/N*pow(p, 2)))/pow(sqrttrm, 3)))
    def Hyy(self, v, w, p, q):
        return 0.
