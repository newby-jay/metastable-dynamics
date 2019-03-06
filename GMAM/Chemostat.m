%%%% mex -output Chemostat_gmam_mex gmam.c Chemostat.c Chemostat_gmam_mex.c

%% Model parameters
alpha1 = 2.0;
alpha2 = 3.0;
ey = 1e-3;
pars = [alpha1, alpha2, ey];

%% fixed points
xfp = alpha1*(alpha2 - 1/(alpha1 - 1));
yfp = 1/(alpha1 - 1);
xS1 = 0;
yS1 = alpha2;
%% GMAM parameters
Npnts = 500; % number of (interior) points along trajectory
Niter = 20000; % Max number of iterations for the gmam method
% time step (not physical time) for iterations (make smaller if convergence is an issue)
dt = 1e-3;
% error tolerance for convergence (max distance between nodes of two successive iterations)
% Note that iterations will stop if either tolerance is reached or max
% number of iterations is reached
abs_err = 1e-13;



%% GMAM
% Initial guess for solution
x = linspace(xfp, xS1, Npnts+1);
y = linspace(yfp, yS1, Npnts+1);
p = linspace(0, 0, Npnts+1);
q = linspace(0, 0, Npnts+1);
ray1 = Chemostat_gmam(x, y, p, q, Npnts, Niter, dt, abs_err, pars);
ray1.error

%% Plotting
figure(1);
clf;
hold on

[X, Y]  = meshgrid(linspace(0, 12, 30), linspace(0, 4, 30));
U = alpha1*(Y./(1 + Y)).*X - X;
V = -(Y./(1 + Y)).*X - Y + alpha2;
S = sqrt(U.^2 + V.^2);
Q = quiver(X, Y, U./S, V./S);

plt = plot(ray1.x, ray1.y, 'r');
plt.LineWidth = 2;


plot(xfp, yfp, 'ko')
plot(xS1, yS1, 'ko')

axis([0 12 0 4]);
hold off

figure(2);
clf;
hold on;
plt = plot(ray1.y, ray1.phi, 'r');
plt.LineWidth = 2;


x = linspace(xfp, 0, Npnts+1);
y = linspace(yfp, 1.5, Npnts+1);
ray = Chemostat_gmam(x, y, p, q, Npnts, Niter, dt, abs_err, pars);


hold off;
