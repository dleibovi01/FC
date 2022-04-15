clear all
close all

n = 1000; % initially 4000
x_a = 0; x_b = 1; % The beginning and end of the Cartesian grid
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
K = 1e3;
% u = @(t) K*t.^3;
% u_der = @(t) K*3*t.^2;
% u_der_2 = @(t) K*6*t;
u = 1000*(x < 0.5) + 0*(x >= 0.5);

[~, ~, fcont] = fc_der(u, 0, 0, prd, k, d, C)
[~, ~, fcont_der] = fc_der(u, 1, 0, prd, k, d, C)
[~, ~, fcont_der_2] = fc_der(u, 2, 0, prd, k, d, C)
s = [x', x_cont'];

figure
plot(s, fcont)

figure
plot(s, fcont_der)

% figure
% plot(s, abs(u_der_2(s) - fcont_der_2))
