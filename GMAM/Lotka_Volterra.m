%%%% To compile, you need to install the mex toolboxes (I think)
%%%% The command is "mex -output Lotka_Volterra_mex *.c"

%% Model parameters
alpha = 4.0;
beta = 1.0;
delta = 0.5;
gamma = 2.5;
rho = 0.2;
pars = [alpha, beta, delta, gamma, rho];

%% fixed points
xfp = gamma/delta;
yfp = (alpha - rho*gamma/delta)/beta;
xS1 = 0;
yS1 = 0;
xS2 = alpha/rho;
yS2 = 0;
%% GMAM parameters
Npnts = 500; % number of (interior) points along trajectory
Niter = 20000; % Max number of iterations for the gmam method
% time step (not physical time) for iterations (make smaller if convergence is an issue)
dt = 1e-2; 
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
ray1 = Lotka_Volterra_gmam(x, y, p, q, Npnts, Niter, dt, abs_err, pars);
ray1.error

x = linspace(xfp, xS2, Npnts+1);
y = linspace(yfp, yS2, Npnts+1);
ray2 = Lotka_Volterra_gmam(x, y, p, q, Npnts, Niter, dt, abs_err, pars);
ray2.error

%% Plotting
figure(1);
clf;
hold on

[X, Y]  = meshgrid(linspace(0, xS2+1, 30), linspace(0, 15, 30));
U = alpha*X - beta*X.*Y - rho*X.^2;
V = delta*X.*Y - gamma*Y;
S = sqrt(U.^2 + V.^2);
Q = quiver(X, Y, U./S, V./S);

plt = plot(ray1.x, ray1.y, 'r');
plt.LineWidth = 2;

plt = plot(ray2.x, ray2.y, 'r');
plt.LineWidth = 2;

plot(xfp, yfp, 'ko')
plot(xS1, yS1, 'ko')
plot(xS2, yS2, 'ko')
hold off
