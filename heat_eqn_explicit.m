clear;clc;
%% Set up the domain data structures
n = 4000; % 4000 initially
x_a = 0; x_b = 2*pi; % The beginning and end of the Cartesian grid
h = (x_b - x_a)/(n-1);

m =  4; % Gram Polynomial max order
d =  5; % Number of orthogonal polynomial matching points
C = 25; % Number of continuation points

% The grid we have data on and wish to interpolate/compute derivatives/solve DEs
x = (x_a:h:x_b).';
% The grid of continuation points (for plotting purposes)
x_cont = x_b + h*(1:C).';

% Period of our FC interpolating Fourier Series -- note the extended period
fourPts = n + C;
prd = fourPts*h;
if (mod(fourPts, 2) == 0)
    k = transpose([0:fourPts/2, -fourPts/2+1:-1]);
else
    k = transpose([0:(fourPts-1)/2, -(fourPts-1)/2:-1]);
end

%% Set up the initial condition of the explicit ODE
uexact = @(x, t) exp(-4 * t) * sin(2*x);
u0x = uexact(x, 0);

tic;
deltat = 10^-8; % initially 1e-8
maxit = 2500;
u = zeros(length(x), maxit);
u(:, 1) = u0x;
for it=2:maxit
    % Forward Euler method
    %uder = fc_der(u(:, it-1), [1, 2], 1, prd, k, d, C);
    %u(:, it) = u(:, it-1) + deltat*uder(:, 2);
    %u(1, it) = 0;
    %u(end, it) = 0;

    % Classical RK4 method
    k1 = deltat * fc_der(u(:, it-1), ...
                    2, 1, prd, k, d, C);
    k2 = deltat * fc_der(u(:, it-1) + 1/2 * k1, ...
                    2, 1, prd, k, d, C);
    k3 = deltat * fc_der(u(:, it-1) + 1/2 * k2, ...
                    2, 1, prd, k, d, C);
    k4 = deltat * fc_der(u(:, it-1) + k3, ...
                    2, 1, prd, k, d, C);
    u(:, it) = u(:, it-1) + 1/6*k1 + 1/3*k2 + 1/3*k3 + 1/6*k4;
    u(1, it) = 0;
    u(end, it) = 0;
end
toc;

plottime = maxit*deltat;

figure(1)
plot(x, u(:, maxit));
title(['Solution at t = ', num2str(plottime)]);

figure(2)
err = u(:, maxit) - uexact(x, plottime);
semilogy(x, abs(err));
title(['Error at t = ', num2str(plottime)]);
fprintf('Maximum error: %1.3e\n', max(abs(err)));
