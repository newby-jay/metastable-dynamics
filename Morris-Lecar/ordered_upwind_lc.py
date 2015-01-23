# Written by Jay Newby 2013
# Send questions and comments to newby@math.utah.edu.

from pylab import *
from string import Template
# import scipy.weave
from scipy.weave import inline
from itertools import product
from scipy.interpolate import interp1d

class OUM:
    ## needs Z(), x0, y0, Hx(), Hp(), Hy(), Hq() from NeuronAP
    def __init__(self, model, x, y, args, gamma0 = 1, sx = 1, sy = 1, debug = False):
        """
        Input
        -----
        model : class containing methods for Z() (the Hessian of Phi), and the four derivatives of H.
        x, y : arays of x and y grid points
        args: filename and delta, the size of the initial data region
        gamma0: the maximum update range for each point
        optional prameters : sx = 1, sy = 1, scaling parameters for the relative distance between nearest neighbor grid points
        """
        self.model = model
        for k, v in model.__dict__.iteritems():
            setattr(self, k, v)
        for k in args:
            setattr(self, k, args[k])
        # self.Z = model.Z
        self.Hx, self.Hy, self.Hp, self.Hq = model.Hx, model.Hy, model.Hp, model.Hq
        if debug:
            self._id = '//' + str(rand(4))
        else:
            self._id = '//new2V0.00'
        cfile = open('oum.c')
        code_template = Template(cfile.read())
        code = code_template.substitute(self.func_string)
        self.scode = code
        cfile.close()
        self.x, self.y  = x, y
        self.nx, self.ny = x.size, y.size
        self.meshdim = (self.ny, self.nx)
        dx, dy = x[1] - x[0], y[1] - y[0]
        self.X, self.Y = meshgrid(x, y)
        self.J, self.I = meshgrid(arange(self.nx), arange(self.ny))
        self.adj = array([array(xi) for xi in product((-1, 0, 1), (-1, 0, 1)) if (xi[0] != 0 or xi[1] != 0)])
        self.nn = array([[0, 1], [1, 0], [0, -1], [-1, 0]])
        self.nn = self.adj
        self.set_gamma(gamma0, sx, sy)
        self.unconsidered = ones(self.meshdim).astype(bool)
        self.dx, self.dy = nan*ones(self.meshdim), nan*ones(self.meshdim)
        # self.initialize()
        # self.set_structure()
    def set_gamma(self, gamma, sx = 1., sy = 1.):
        self.gamma = gamma
        if gamma == 1:
            G = array([xi for xi in product(arange(-gamma, gamma+1), arange(-gamma, gamma+1))])
        else:
            G = array([xi for xi in product(arange(-gamma, gamma+1), arange(-gamma, gamma+1)) if sqrt((xi[0]/sx)**2 + (xi[1]/sy)**2) <= gamma])
        gamma_range = array([[i, j]  for i, j in G if i != 0 or j != 0])
        d = array([sqrt((xi[0]/sx)**2 + (xi[1]/sy)**2) for xi in gamma_range])
        self.gamma_range = gamma_range[d.argsort()]
    def initialize(self):
        Npts = 2000
        xc = self.xc
        L = xc['s'][-1]
        taux, tauy = self.Hp(xc['x'], xc['y'], 0, 0), self.Hq(xc['x'], xc['y'], 0, 0) # tangent vector
        B = sqrt(taux**2 + tauy**2) # speed on lc
        thatx, thaty = taux/B, tauy/B # unit tangent
        nhatx, nhaty = thaty, -thatx # unit normal (outward)
        Jxx = self.model.detjacxx(xc['x'], xc['y'])  # elements of the deterministic jacobian
        Jxy = self.model.detjacxy(xc['x'], xc['y'])
        Jyx = self.model.detjacyx(xc['x'], xc['y'])
        Jyy = self.model.detjacyy(xc['x'], xc['y'])
        Cc = Jxx*nhatx**2 + (Jxy + Jyx)*nhatx*nhaty + Jyy*nhaty**2 # inward speed normal to lc
        Dxx = self.model.Hpp(xc['x'], xc['y'], 0, 0) # diffusion tensor
        Dxy = self.model.Hpq(xc['x'], xc['y'], 0, 0)
        Dyy = self.model.Hqq(xc['x'], xc['y'], 0, 0)
        Dc = Dxx*nhatx**2 + 2.*Dxy*nhatx*nhaty + Dyy*nhaty**2 # normal diffusivity
        # Implicit solver to get periodic solution to a first order linear problem, using centered finite difference
        # 1/2*B*(d/ds)alpha - Cc*alpha = Dc, alpha(0) = alpha(L)
        s = linspace(0, L, Npts) # discrete points on s
        ds = s[1] - s[0]
        Bf = interp1d(xc['s'], B) # everything with f or I at the end is interpolated to s points
        BI = Bf(s)
        Ccf = interp1d(xc['s'], Cc)
        CcI = Ccf(s)
        Dcf = interp1d(xc['s'], Dc)
        DcI = Dcf(s)
        tcf = interp1d(xc['s'], xc['t'])
        tcI = tcf(s)
        xcf, ycf = interp1d(xc['s'], xc['x']), interp1d(xc['s'], xc['y'])
        xcI, ycI = xcf(s), ycf(s)
        nhatxf, nhatyf = interp1d(xc['s'], nhatx), interp1d(xc['s'], nhaty)
        nhatxI, nhatyI = nhatxf(s), nhatyf(s)
        thatxf, thatyf = interp1d(xc['s'], thatx), interp1d(xc['s'], thaty)
        D1 = diag(ones(Npts-1), 1) - diag(ones(Npts-1), -1) # centered difference matrix
        D1[0, -1] = -1. # periodic BC
        D1[-1, 0] = 1.
        D1 = D1/(2.*ds)
        M = 0.5*dot(diag(BI), D1) - diag(CcI) # operator for alpha = 1/phi
        phirr = 1./solve(M, DcI) # phi = 1/alpha
        phirrf = interp1d(s, phirr)
        dphirrds = dot(D1, phirr)
        dphirrdsf = interp1d(s, dphirrds)
        # 'x': xcI - r0*nhatxI,
        # 'y': ycI - r0*nhatyI,
        def srofxy(x, y):
            d = sqrt((xcI - x)**2 + (ycI - y)**2)
            ns = d.argmin()
            sn = s[ns]
            rn = -d[ns]*sign(nhatxI[ns]*(xcI[ns] - x) + nhatyI[ns]*(ycI[ns] - y))
            return sn, rn
        S, R = inf*ones(self.meshdim), inf*ones(self.meshdim)
        for i in arange(self.ny):
            for j in arange(self.nx):
                S[i, j], R[i, j] = srofxy(self.X[i, j], self.Y[i, j])
        self.phi, self.p, self.q = nan*ones(self.meshdim), nan*ones(self.meshdim), nan*ones(self.meshdim)
        Phi1 = 0.5*R**2*phirrf(S)
        inds = (Phi1 <= self.delta)# + (R < 0)
        self.phi[inds] = Phi1[inds]
        phirr2 = phirrf(S[inds])
        dphirrds2 = dphirrdsf(S[inds])
        self.p[inds] = R[inds]*phirr2*nhatxf(S[inds]) + 0.5*R[inds]**2*dphirrds2*thatxf(S[inds])
        self.q[inds] = R[inds]*phirr2*nhatyf(S[inds]) + 0.5*R[inds]**2*dphirrds2*thatyf(S[inds])
        self.accepted_front = zeros(self.meshdim).astype(bool)
        self.accepted = zeros(self.meshdim).astype(bool)
        self.accepted[inds] = True
        self.set_structure()
        
        # f4 = figure(4)
        # clf()
        # f4.add_subplot(131)
        # plot(s, phirr)
        # f4.add_subplot(132)
        # pcolor(self.X, self.Y, R)
        # colorbar()
        # f4.add_subplot(133)
        # Phi1[Phi1 > 0.5*self.delta] = 0
        # pcolor(self.X, self.Y, Phi1)
        # colorbar()
    def set_structure(self):
        self.unconsidered[self.accepted] = False
        self.phic = inf*ones(self.meshdim)
        self.considered = zeros(self.meshdim).astype(bool)
        ## find considered points and tentative accepted_front points
        for ij, dij in product(zip(self.I[self.accepted].flatten(), self.J[self.accepted].flatten()), self.adj):
            k, l = ij + dij
            ind = (k, l)
            if 0 < k < self.ny and 0 < l < self.nx and self.unconsidered[ind]:
                self.considered[ind] = 1
                self.accepted_front[ij] = 1
        self.unconsidered[self.considered] = 0
        self.accepted[self.accepted_front] = 0  # remove accepted_front points from accepted
        ## refine accepted_front: find accepted_front points that do not have a nearest neigbor considered point and move them back to accepted
        for ij in zip(self.I[self.accepted_front].flatten(), self.J[self.accepted_front].flatten()):
            i, j = ij
            s = 1
            for k, l in ij + self.nn:
                if 0 < k < self.ny and 0 < l < self.nx and self.considered[(k, l)]:
                    s = 0
                    break
            if s and 0 < i < self.ny and 0 < j < self.nx:
                self.accepted_front[ij] = 0
                self.accepted[ij] = 1
        self.compute_considered(self.considered)
    def compute_considered(self, inds):
        ## weave wrapper for compute_considered
        code = self._id + """
        compute_considered(numx, numy, 
                             epsilon, Nna, Nk, gna, vna, gk, vk, gl, 
                             vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2, 
                             x, y,
                             iadj, jadj, igamma_range, jgamma_range, ijconsidered,
                             nadj, ngamma_range, nconsidered,
                             phic, phi, p, q, dx, dy,
                             accepted_front);
                             """
        ijconsidered = (self.I[inds].flatten()*self.nx + self.J[inds].flatten()).astype(double)
        nconsidered = ijconsidered.size
        iadj,  jadj = self.adj[:, 0].astype(double), self.adj[:, 1].astype(double)
        igamma_range,  jgamma_range = self.gamma_range[:, 0].astype(double), self.gamma_range[:, 1].astype(double)
        nadj, ngamma_range = iadj.size, igamma_range.size
        epsilon, Nna, Nk, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2 = self.parameters
        numx, numy, x, y = self.nx, self.ny, self.x.astype(double), self.y.astype(double)
        accepted_front = self.accepted_front.astype(bool)
        phic, phi, p, q, dx, dy = self.phic.astype(double), self.phi.astype(double), self.p.astype(double), self.q.astype(double), self.dx.astype(double), self.dy.astype(double)
        inline(code, 
                           ['numx', 'numy', 'epsilon', 'Nna', 'Nk', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'x', 'y', 'iadj',  'jadj',  'igamma_range',  'jgamma_range', 'ijconsidered', 'nadj', 'ngamma_range', 'nconsidered', 'phic', 'phi', 'p', 'q', 'dx', 'dy', 'accepted_front'], 
                                   support_code = self.scode, 
                                   libraries=['gsl'], 
                                   include_dirs=['/opt/local/include'],
                                   library_dirs=['/opt/local/lib'],
                                   auto_downcast = 1)
        self.accepted_front = accepted_front
        self.phic, self.p, self.q, self.dx, self.dy = phic, p, q, dx, dy
    def advance_solution(self, nups):
        ## weave wrapper for compute_considered
        code = self._id + """
                //double *psub = p, *qsub = q, *npsub = nump, *nqsub = numq, *phisub = phi, *phicsub = phic;
                //bool *unconsidered_sub = unconsidered, *considered_sub = considered, *accepted_front_sub = accepted_front, *accepted_sub = accepted_front;
                advance_soln(nups, numx, numy, epsilon, Nna, Nk, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2, x, y, iadj,  jadj, inn,  jnn,  igamma_range,  jgamma_range, nadj, nnn, ngamma_range, phic, phi, p, q, dx, dy, unconsidered, considered, accepted_front, accepted);
                """
        iadj,  jadj = self.adj[:, 0].astype(double), self.adj[:, 1].astype(double)
        inn,  jnn = self.nn[:, 0].astype(double), self.nn[:, 1].astype(double)
        igamma_range,  jgamma_range = self.gamma_range[:, 0].astype(double), self.gamma_range[:, 1].astype(double)
        nadj, nnn, ngamma_range = iadj.size, inn.size, igamma_range.size
        epsilon, Nna, Nk, gna, vna, gk, vk, gl, vl, Iapp, betak, thetak1, thetak2, thetana1, thetana2 = self.parameters
        numx, numy, x, y = self.nx, self.ny, self.x.astype(double), self.y.astype(double)
        unconsidered, considered, accepted_front, accepted = self.unconsidered.astype(bool), self.considered.astype(bool), self.accepted_front.astype(bool), self.accepted.astype(bool)
        phic, phi, p, q, dx, dy = self.phic.astype(double), self.phi.astype(double), self.p.astype(double), self.q.astype(double), self.dx.astype(double), self.dy.astype(double)
        inline(code, 
                           ['nups', 'numx', 'numy', 'epsilon', 'Nna', 'Nk', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'x', 'y', 'iadj',  'jadj', 'inn',  'jnn',  'igamma_range',  'jgamma_range', 'nadj', 'nnn', 'ngamma_range', 'phic', 'phi', 'p', 'q', 'dx', 'dy', 'unconsidered', 'considered', 'accepted_front', 'accepted'], 
                                   support_code = self.scode, 
                                   libraries=['gsl'], 
                                   include_dirs=['/opt/local/include'],
                                   library_dirs=['/opt/local/lib'],
                                   auto_downcast = 1)
        self.unconsidered, self.considered, self.accepted_front, self.accepted = unconsidered, considered, accepted_front, accepted
        self.phic, self.phi, self.p, self.q, self.dx, self.dy = phic, phi, p, q, dx, dy
        self.data = {'phi': self.phi, 'p': self.p, 'q': self.q, 'dx': self.dx, 'dy': self.dy, 'x': self.x, 'y': self.y}
        np.save('data/'+self.filename+'_data', array([self.data]))
    def phicset(self, gamma):
        # reset the value of gamma, recompute current considered points that are inf
        self.set_gamma(gamma)
        inds = (self.phic == inf)*self.considered
        self.compute_considered(inds)
    def plot_points(self, sfn, fn = 1, fy = 2, fx = 3):
        # plots a grid of symbols to track what population each grid point is in (i.e., unconsidered, considered, accepted, or AF)
        plot(self.x0, self.y0, 'bo', ms = 3)
        plot(self.X[self.accepted], self.Y[self.accepted], linestyle = 'none', mec='orange', mfc = 'orange', marker = 'x', ms = 5)
        plot(self.X[self.accepted_front], self.Y[self.accepted_front], linestyle = 'none', mec = 'r', mfc = 'r', marker = 'o', ms = 5)
        inds = (self.phic == inf)*self.considered
        inds2 = (self.phic != inf)*self.considered
        plot(self.X[inds2], self.Y[inds2], linestyle = 'none', mec = 'c', mfc = 'c', marker = 'x', ms = 3)
        plot(self.X[inds], self.Y[inds], linestyle = 'none', mec = 'b', mfc = 'b', marker = 'x', ms=3)
        axis([self.x[0], self.x[-1], self.y[0], self.y[-1]])

###################################### script ####################################################

# if __name__ == '__main__':
# dZ = dt*(-Z*D*Z - Z*C - C^T*Z - G)
# 
# Kdot = gradK . xdot = gradK . Hp
# Kdot = Hp/(l^T*F*r)*l^T*(F*r)_x
