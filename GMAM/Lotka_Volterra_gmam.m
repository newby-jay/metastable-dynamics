function sol = Lotka_Volterra_gmam(x, y, p, q, Npnts, Niter, dt, abs_err, pars)
    [sol0, err] = Lotka_Volterra_gmam_mex(x, y, p, q, Npnts, Niter, dt, abs_err, pars);
    sol.t = sol0(:, 1);
    sol.s = sol0(:, 2);
    sol.x = sol0(:, 3);
    sol.y = sol0(:, 4);
    sol.p = sol0(:, 5);
    sol.q = sol0(:, 6);
    sol.phi = sol0(:, 7);
    sol.error = err;
    sol.dt = dt;
    sol.abs_err = abs_err;
    sol.pars = pars;
end