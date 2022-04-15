% My version of blend_to_zero

function [Q, A] =  my_blend_to_zero(d, C, Z, E, n_over, num_digits, svd_order)


tol = 1e-16;
digits(num_digits);
A = zeros(C, d);
N_coarse = d + C + Z + E;
if (nargin == 6)
    svd_order = N_coarse;
end

%Forming the grids we will use
h = 1/N_coarse; 
interp_coarse = sym(h*(0 : d - 1).'); 
interp_fine = sym(h*linspace(0, d - 1, n_over*(d - 1) + 1).'); 
zero_coarse = sym(h*(d + C : d - 1 + C + Z).');
zero_fine = sym(h*linspace(d + C, d - 1 + C + Z, n_over*(Z - 1) + 1).'); 
continuation = sym(h*(d : d - 1 + C).');


% Forming Q by using the modified Gram-Schmidt algorithm
fprintf( 'Performing QR decomposition... '); tic;
P_coarse = sym(zeros(d, d));
for i = 0 : d - 1
    P_coarse(:, i + 1) = interp_coarse.^i;
end
[Q, R] = mgs(P_coarse);
Q = double(Q);

% Forming Q_fine
P_fine = sym(zeros(n_over*(d - 1) + 1, d));
for i = 0 : d - 1
   P_fine(:, i + 1) = interp_fine.^i; 
end
Q_fine = P_fine*inv(R);
fprintf( 'Done. (%1.3fs)\n', toc );

% Approximate the Gram polynomials by trigonometric polynomials using an
% SVD
fprintf( 'Performing SVD... '); tic;
if( mod(svd_order,2) == 0 )
    k_cos = sym(0 : svd_order/2 - 1);
    k_sin = sym(1 : svd_order/2 - 1);
    k = sym(-svd_order/2 : svd_order/2-1);
else
    k_cos = sym(0 : (svd_order - 1)/2);
    k_sin = sym(1 : (svd_order - 1)/2);
    k = sym(-(svd_order-1)/2 : (svd_order - 1)/2);
end
X = [interp_fine; zero_fine];
C = [cos(2*sym('pi')*X*k_cos), sin(2*sym('pi')*X*k_sin)];
% C = exp(2 * 1i * sym('pi') * X * k);
[U, S, V] = svd(C, 'econ');
fprintf( 'Done. (%1.3fs)\n', toc );
Coeffs = sym(zeros(size(C, 2), d));
index = 1;
delta = diag(S);
delta_inv = 1./delta;


for i = 1 : d
    b = [Q_fine(:, i); sym(zeros(n_over*(Z - 1) + 1, 1))];
    Coeffs(:, i) = V*(delta_inv.*(U'*b)); % inverting the SVD
    r = max( abs( double( C*Coeffs(:, i) - b ) ) );
    fprintf( '\t%d\tresidual = %e\n', d, r );    
end

% Evaluating the trigonometric polynomials at the continuation points
fprintf( 'Evaluating at continuation points... '); tic;
A = double([cos(2*pi*continuation*k_cos), sin(2*pi*continuation*k_sin)]*Coeffs);
fprintf( 'Done. (%1.3fs)\n', toc );
end