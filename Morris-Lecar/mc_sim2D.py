## Jay Newby, 2013
## newby.23@mbi.osu.edu
from pylab import *
import time
import scipy.weave
from scipy.interpolate import interp1d


def simulate_traj(x0, n0, m0, Npoints,
                  epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betana, betak, thetak1, thetak2, thetana1, thetana2,
                  **kwargs):
    cfile = open('montecarlo2D.c')
    scode = cfile.read()
    cfile.close()
    t0 = 0.
    solt, solv, soln, solm = (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double)
    instid =  '// traj2D  V0.00' #+ str(rand(4))
    code =    instid + """
              monte_carlo_2D_traj(t0, n0, m0, x0, Npoints, epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betana, betak, thetak1, thetak2, thetana1, thetana2, solt, solv, soln, solm);
              """
    tm1 = time.time()
    scipy.weave.inline(code,
                       ['t0', 'x0', 'n0', 'm0', 'epsilon', 'N', 'M', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'solt', 'solv', 'soln', 'solm', 'Npoints', 'betana'],
                       support_code = scode,
                       libraries=['gsl'],
                       include_dirs=['/opt/local/include'],
                       library_dirs=['/opt/local/lib'])
    tm2 = time.time()
    # print 'compute time = ', (tm2-tm1)
    inds = ~isnan(solt)
    traj = {'t': solt[inds], 'v': solv[inds], 'n': soln[inds], 'w': solm[inds]/M}
    inds = traj['t'].argsort()
    for k in traj:
        traj[k] = traj[k][inds]
    return traj

def simulate_traj_sep(x0, n0, m0, Npoints, sep,
                      epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betana, betak, thetak1, thetak2, thetana1, thetana2,
                      **kwargs):
    cfile = open('montecarlo2D.c')
    scode = cfile.read()
    cfile.close()
    t0 = 0.
    solt, solv, soln, solm = (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double), (nan*zeros(Npoints)).astype(double)
    sep = sep.astype(double)
    instid =  '// traj 2D_sep V0.00' #+ str(rand(4))
    code =    instid + """
              monte_carlo_2D_traj_sep(t0, n0, m0, x0, Npoints, sep, epsilon, N, M, gna, vna, gk, vk, gl, vl, Iapp, betana, betak, thetak1, thetak2, thetana1, thetana2, solt, solv, soln, solm);
              """
    tm1 = time.time()
    scipy.weave.inline(code,
                       ['t0', 'x0', 'n0', 'm0', 'epsilon', 'N', 'M', 'gna', 'vna', 'gk', 'vk', 'gl', 'vl', 'Iapp', 'betana', 'betak', 'thetak1', 'thetak2', 'thetana1', 'thetana2', 'solt', 'solv', 'soln', 'solm', 'Npoints', 'sep'],
                       support_code = scode,
                       libraries=['gsl'],
                       include_dirs=['/opt/local/include'],
                       library_dirs=['/opt/local/lib'])
    tm2 = time.time()
    outstr = 'compute time = ' +  str(tm2-tm1)
    print outstr
    sys.stdout.write("\033[F")
    time.sleep(.001)
    inds = ~isnan(solt)
    traj = {'t': solt[inds], 'v': solv[inds], 'n': soln[inds], 'w': solm[inds]/M}
    inds = traj['t'].argsort()
    for k in traj:
        traj[k] = traj[k][inds]
    return traj
