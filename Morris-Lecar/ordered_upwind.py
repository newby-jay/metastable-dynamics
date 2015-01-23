# Written by Jay Newby 2013
# Send questions and comments to newby@math.utah.edu.
from pylab import *
import time
from string import Template
import scipy.weave
from itertools import product

class OUM:
    ## needs Z(), x0, y0, Hx(), Hp(), Hy(), Hq() from MLexcite
    def __init__(self, model, x, y, args, gamma0 = 1, sx = 1, sy = 1, debug = False):
        """
        Input
        -----
        model : class containing methods for Z() (the Hessian of Phi), the fixed point, and the four partial derivatives of H.
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
        self.Z = model.Z
        self.Hx, self.Hy, self.Hp, self.Hq = model.Hx, model.Hy, model.Hp, model.Hq
        if debug:
            self._id = '//' + str(rand(4))
        else:
            self._id = '//newV0.15'
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
        self.initialize()
        self.set_structure()
    def set_gamma(self, gamma, sx = 1., sy = 1.):
        "Set the max range for FD updates.  Optional variables sx and sy scale the x and y relative distance"
        self.gamma = gamma
        if gamma == 1:
            G = array([xi for xi in product(arange(-gamma, gamma+1), arange(-gamma, gamma+1))])
        else:
            G = array([xi for xi in product(arange(-gamma, gamma+1), arange(-gamma, gamma+1)) if sqrt((xi[0]/sx)**2 + (xi[1]/sy)**2) <= gamma])
        gamma_range = array([[i, j]  for i, j in G if i != 0 or j != 0])
        d = array([sqrt((xi[0]/sx)**2 + (xi[1]/sy)**2) for xi in gamma_range])
        self.gamma_range = gamma_range[d.argsort()]
    def initialize(self):
        "Set everything up"
        Z = self.Z(self.x0, self.y0)
        def PHI0(x, y):
            return 0.5*((self.x0 - x)**2*Z[0, 0] + 2*(self.x0 - x)*(self.y0 - y)*Z[0, 1] + (self.y0 - y)**2*Z[1, 1])
        def p0(x, y):
            return -Z[0, 0]*(self.x0 - x) - Z[0, 1]*(self.y0 - y)
        def q0(x, y):
            return -Z[1, 0]*(self.x0 - x) - Z[1, 1]*(self.y0 - y)
        self.phi, self.p, self.q = nan*ones(self.meshdim), nan*ones(self.meshdim), nan*ones(self.meshdim)
        Phi1 = PHI0(self.X, self.Y)
        inds = Phi1 <= 0.5*self.delta
        self.phi[inds] = Phi1[inds]
        self.p[inds] = p0(self.X[inds], self.Y[inds])
        self.q[inds] = q0(self.X[inds], self.Y[inds])
        self.accepted_front = zeros(self.meshdim).astype(bool)
        self.accepted = zeros(self.meshdim).astype(bool)
        self.accepted[inds] = True
    def set_structure(self):
        "Determine which grid points are accepted, accepted front, considered, and unconsidered."
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
        "Compute a tentative value of Phi for every considered point"
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
        scipy.weave.inline(code, 
                           ['numx', 'numy', 'epsilon', 'Nna', 'Nk', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'x', 'y', 'iadj',  'jadj',  'igamma_range',  'jgamma_range', 'ijconsidered', 'nadj', 'ngamma_range', 'nconsidered', 'phic', 'phi', 'p', 'q', 'dx', 'dy', 'accepted_front'], 
                                   support_code = self.scode, 
                                   libraries=['gsl'], 
                                   include_dirs=['/opt/local/include'],
                                   library_dirs=['/opt/local/lib'],
                                   auto_downcast = 1)
        self.accepted_front = accepted_front
        self.phic, self.p, self.q, self.dx, self.dy = phic, p, q, dx, dy
    def advance_solution(self, nups):
        """Main OUM function.
           Parameters
           ------------
           nups: the number of iterations to perform. Each iteration determines one new accepted front point from among the considered point.
           ------------
           The result is stored in the data atribute.
           """
        code = self._id + """
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
        scipy.weave.inline(code, 
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
    def plot_points(self, sfn, fn = 1, fy = 2, fx = 3):
        "plots a grid of symbols to track what population each grid point is in (i.e., unconsidered, considered, accepted, or AF)"
        plot(self.x0, self.y0, 'bo', ms = 3)
        # plot(self.X[self.unconsidered], self.Y[self.unconsidered], 'kx', ms = 2)
        plot(self.X[self.accepted], self.Y[self.accepted], linestyle = 'none', mec='orange', mfc = 'orange', marker = 'x', ms = 5)
        plot(self.X[self.accepted_front], self.Y[self.accepted_front], linestyle = 'none', mec = 'r', mfc = 'r', marker = 'o', ms = 5)
        inds = (self.phic == inf)*self.considered
        inds2 = (self.phic != inf)*self.considered
        plot(self.X[inds2], self.Y[inds2], linestyle = 'none', mec = 'c', mfc = 'c', marker = 'x', ms = 3)
        plot(self.X[inds], self.Y[inds], linestyle = 'none', mec = 'b', mfc = 'b', marker = 'x', ms=3)
        axis([self.x[0], self.x[-1], self.y[0], self.y[-1]])
