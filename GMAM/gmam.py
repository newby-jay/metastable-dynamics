# Written by Jay Newby 2014
# Send questions and comments to newby@math.utah.edu.

from pylab import *
import weave

class GMAM:
    def __init__(self, debug = False):
        """
        Input
        ------------
        debug = False: if True then code is recompiled before each execution
        """
        if debug:
            self._id = '//' + str(rand(4))
        else:
            self._id = '//V0.05'
        with open('gmam.c', 'r') as cfile:
            code = cfile.read()
            self.scode = code
    def __call__(self, Npnts, Niter, dt, x, y, pars, tol = 1e-13):
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
        pars = float64(pars)
        dt = double(dt)
        tol = double(tol)
        Npnts = int(Npnts)
        Niter = int(Niter)
        assert Npnts > 0
        assert Niter > 0
        assert dt > 0
        assert tol > 0
        assert x.size == Npnts + 1
        assert x.size == y.size
        x, y = float64(x), float64(y)
        t, s = zeros(Npnts+1), zeros(Npnts+1)
        p, q = zeros(Npnts+1), zeros(Npnts+1)
        phi = zeros(Npnts+1)
        error = float64([0.0])
        code = self._id + """
        gmam(t, s, x, y, p, q, phi, Npnts, Niter, dt, tol, error, pars);"""
        cvars = [
            't', 's', 'x', 'y', 'p', 'q', 'phi',
            'Npnts', 'Niter', 'dt', 'pars', 'tol', 'error']
        weave.inline(
            code,
            cvars,
            support_code = self.scode,
            auto_downcast = 1)
        return {
            't': t, 's': s, 'x': x, 'y': y,
            'p': p, 'q': q, 'phi': phi, 'error': error}
