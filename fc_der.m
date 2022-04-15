function [fx_der, fx_der_coeffs, fx_der_cont] = fc_der(fx, der_order, filter, prd, k, d, C)

% FC_DER - Routine to compute the derivative of a function using Fourier
% continuation.
%   [fx_der fx_der_coeffs] = fc_der(fx, der_order, prd, k, d, C)
%   returns the derivative (fx_der) and its coefficients (fx_der_coeffs) of
%   the function f, given its values fx on an equispaced grid
%
%   Inputs:
%       f_x - value of the function f at the nodes x_j (real n-vector)
%       der_order - the order of the derivative desired (integer array)
%       prd - structure containing parameters of the domain, including the
%       k - wavenumber vector
%       d - number of FC matching points
%       C - number of FC continuation points
%
%   Outputs:
%       fx_der - the (der_order)^th derivative of the function computed at
%       the grid points x_j  (real, n-vector)
%       fx_der_coeffs - the (continuation) Fourier coefficients of the
%       derivative (complex, n+C length vector, where n = length(x_j))
%
%   Author(s):
%       Thomas G. Anderson
%       Email: tanderson@caltech.edu
%

% Compute the Fourier continuation coefficients (fx_cont_coeffs) of the
% function from its grid point values (fx) using the FC(Gram) method.
[fx_cont_coeffs f_dp] = fcont_gram_blend(fx, d, C);
fourPts = length(fx) + C;
n = length(fx);
if (filter ~= 0)
    fx_cont_coeffs = fx_cont_coeffs .* specFilter(2*k/fourPts);
end

% Compute the coefficients of the (der_order)^th derivative of the function
fx_der_coeffs = zeros(length(der_order), fourPts);
fx_der_cont = zeros(fourPts, 1);
fx_der = zeros(n, length(der_order));
for i=1:length(der_order)
    fx_der_coeffs(i,:) = fx_cont_coeffs.*(2*pi*1j*k/prd).^der_order(i);

    % Compute the (der_order)^th derivative of the function in physical-space
    fx_der_cont = fourPts*real(ifft(fx_der_coeffs(i,:)));
    fx_der(:, i) = fx_der_cont(1:n);
end

return
