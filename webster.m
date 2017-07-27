clear all;
close all;

data = load('Areafunktio.mat');

% mesh
p = data.ts;
%p = 0:0.005:0.2;
N = size(p,2);

% material parameters
rho = 1.0;
c = 2;

% initialize solution vector
W = poisson_hermite(data.As/max(data.As), p);
M = mass_hermite(1+0*p, p);

N = length(p);
delta = 1e-5;
plotevery = 1e-3;
eps = 1e-8;

% boundary conditions
speed = 20;
dirichlet = @(t)sin(2*pi*t*speed)*(t<0.5/speed);
ddirichlet = @(t)2*pi*speed*cos(2*pi*t*speed)*(t<0.5/speed);
dddirichlet = @(t)-4*pi^2*speed^2*sin(2*pi*t*speed)*(t<0.5/speed);
neumann = @(t)0*cos(t);
dneumann = @(t)-0*sin(t);
ddneumann = @(t)-0*sin(t);

% boundary nodes
d = [1 2*N];
i = setdiff(1:(2*N), d);

% reduction to first-order problem in time
Z = 0*M;
EN = [M(i,i) Z(i,i);
     Z(i,i) M(i,i)];
AN = [Z(i,i) M(i,i)/rho;
     rho*c^2*W(i,i) Z(i,i)];
 
B1 = [-M(i,d)  Z(i,d); 
       Z(i,d) -M(i,d)];

B2 = [Z(i,d), M(i,d)/rho;
      rho*c^2*W(i,d) Z(i,d)];

% initialize
u = zeros(size(EN, 1), 1);
%u(floor(N/2)) = 1;
T = 0;

D = [d d+2*N];
I = setdiff(1:(4*N), D);

figure(1);
clf;
subplot(2,1,1);
plot(data.ts,data.As);
subplot(2,1,2);
U = zeros(4*N, 1);
U(I) = u;
U(D) = [dirichlet(T); neumann(T); ddirichlet(T)/rho; dneumann(T)/rho];
values = U(1:(2*N));
plot_hermite(values,p);
xlim([min(p),max(p)]);
ylim([-1,1]);
drawnow;
pause;

prevplot = 0;

for itr=1:1e7
    T = (itr-1)*delta;
    f = @(t) B2*[dirichlet(t); neumann(t); ddirichlet(t)/rho; dneumann(t)/rho] + B1*[ddirichlet(t); dneumann(t); dddirichlet(t)/rho; ddneumann(t)/rho];
    u = (EN/delta - AN/2)\((EN/delta + AN/2)*u + (f(T+delta) + f(T))/2);
    prevplot = prevplot + delta;
    
    if prevplot>plotevery
        figure(1);
        clf;
        subplot(2,1,1);
        plot(data.ts,data.As);
        subplot(2,1,2);
        U = zeros(4*N, 1);
        U(I) = u;
        U(D) = [dirichlet(T); neumann(T); ddirichlet(T)/rho; dneumann(T)/rho];
        values = U(1:(2*N));
        plot_hermite(values,p);
        xlim([min(p),max(p)]);
        ylim([-1,1]);
        drawnow;
        pause(0.02);
        prevplot = 0;
    end
end

U = zeros(4*N, 1);
U(I) = u;
U(D) = [dirichlet(T); neumann(T); ddirichlet(T)/rho; dneumann(T)/rho];
values = U(1:(2*N));
plot_hermite(values,p);
xlim([min(p),max(p)]);
ylim([min(values(1:N)),max(values(1:N))]);
