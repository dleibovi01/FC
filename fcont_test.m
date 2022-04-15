clear;clc;
%% Set up the domain data structures
n = 4000;
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

%% Evaluate the function and perform continuations
%f = @(y) y.^2 - 1*y + 1;
%fder = @(y) 2*y - 1;
%f2der = @(y) 2;

a = 5.4*pi;
b = 2.7*pi;
c = 2*pi;
% f = @(y) exp(sin(a*y - b) - cos(c*y));
% fder = @(y) (a*cos(a*y - b) + c*sin(c*y)).*exp(sin(a*y - b) - cos(c*y));
% f2der = @(y) -(a^2*sin(a*y - b) - c^2*cos(c*y) - (a*cos(a*y - b) + c*sin(c*y)).^2) .* exp(sin(a*y - b) - cos(c*y));
f = @(y) y.^3;
fder = @(y) 3*y.^2;
f2der = @(y) 6*y;
tic;

fx = f(x);
[fx_cont_coeffs fc fc_l fc_r] = fcont_gram_blend(fx, d, C);
[fx_der fx_der_coeffs] = fc_der(fx, 1, prd, k, d, C);
[fx_2der fx_2der_coeffs] = fc_der(fx, 2, prd, k, d, C);
toc;

figure(1);
plot(x, fx, 'b*', x_cont, fc(n+1:end), 'rd');
xlabel('x');
ylabel('f(x)');

z = linspace(x_a, x_b, 8000);
fc_fine_grid = real(exp(2*pi*i*(z - x_a).'*k.'/prd) * ...
        fx_cont_coeffs(:));
ferr = abs(fc_fine_grid - f(z).');
fprintf('||f(x)||: %1.3e\n', max(ferr));

figure(2)
plot(z, log10(ferr));
xlabel('x', 'Fontsize', 15);
ylabel('log_1_0(e)', 'Fontsize', 15);


figure(3)
plot(1:length(fc_l),fc_l, 1:length(fc_r),fc_r, 1:length(fc_l), fc_l+fc_r);
xlabel('blend region (pts)');
ylabel('f(x)');
title('Blend to zero (zoomed)');


figure(4);
err = abs((fder(x) - fx_der));
fprintf('||f^(1)(x)||: %1.3e\n', max(err));
plot(x, log10(err));
xlabel('x');
ylabel('log_1_0(err)');
title('Error in Derivative');

figure(5);
err = abs((f2der(x) - fx_2der));
fprintf('||f^(2)(x)||: %1.3e\n', max(err));
plot(x, log10(err));
xlabel('x');
ylabel('log_1_0(err)');
title('Error in 2nd Derivative');
