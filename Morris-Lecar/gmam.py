# Written by Jay Newby 2014
# Send questions and comments to newby@math.utah.edu.

from pylab import *
from string import Template
import weave
from scipy.interpolate import interp1d

class GMAM:
    def __init__(self, model, debug = False):
        """
        Input
        ------------
        model: an instance of MLexcite.
        debug = False: if True then code is recompiled before each execution
        """
        self.model = model
        for k, v in model.__dict__.iteritems():
            setattr(self, k, v)
        self.Z0 = model.Z(self.x0, self.y0)
        if debug:
            self._id = '//' + str(rand(4))
        else:
            self._id = '//V0.0'
        cfile = open('gmam.c')
        code_template = Template(cfile.read())
        code = code_template.substitute(self.func_string)
        self.scode = code
        cfile.close()
    def __call__(self, Npnts, Niter, dt, x, y, tol = 1e-13):
        """
        Input
        -------------
        Npnts: number of grid points
        Niter: number of iterations to run
        dt: time step. Should be small enough for stability
        x: initial guess for the x coordinate
        y: initial guess for the y coordinate
        tol = 1e-13: error tolerance
        -------------
        Output
        -------------
        dictionary containing the solution
        """
        code = self._id + """
        gmam(t, s, x, y, p, q, phi, Npnts, Niter, dt, tol, pars);"""
        pars = self.parameters.astype(double)
        x, y = x.copy().astype(double), y.copy().astype(double)
        t, s = zeros(Npnts+1).astype(double), zeros(Npnts+1).astype(double)
        p, q = zeros(Npnts+1).astype(double), zeros(Npnts+1).astype(double)
        phi = zeros(Npnts+1).astype(double)
        weave.inline(code, ['t', 's', 'x', 'y', 'p', 'q', 'phi', 'Npnts', 'Niter', 'dt', 'pars', 'tol'],
                                   support_code = self.scode,
                                   libraries=['gsl'],
                                   include_dirs=['/opt/local/include'],
                                   library_dirs=['/opt/local/lib'],
                                   auto_downcast = 1)
        return {'t': t, 's': s, 'x': x, 'y': y, 'p': p, 'q': q, 'phi': phi}
