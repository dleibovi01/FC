function [A,Q] = fcprecomp(N,d,C,Z,E,nos)
% PRECOMP - Routine to compute the matrices A,Q in (12) in [1].
%   Inputs:
%   N - number of point values provided
%   d-1 - degree at the Gram polynomial at the boundary section
%   C - number of continuation points
% 
%   Outputs:
%   A,Q - A * Q.' the matrix for continuation (see (12) in [1])
%   Note: We can separately return all the four matrices Al, Ar, Ql and Qr.
% 
%   Ref:
%   [1] "An FC-based spectral solver for elastodynamics problems in general
%        three-dimensional domains" by Faisal Amlani and Oscar Bruno, Journal
%        of Computational Physics 307 (2016) 333-354.
% 
% Author(s):
%   Jagabandhu Paul
%   jpaul@caltech.edu
% 

%% Set some comon parameters

% We will use variable precision arithmatic (via sym(.) and vpa(.)) to
% compute the continuation matrices.
fprintf('Precomputations: ... \n');
% set the number of digits of precision
digits(256);

x = linspace(0,1,N);
delta = x(end)-x(end-d+1);
h = 1/(N-1);
% continuation length
c = 1 + (C+1)*h;
% points on the boundary section
xr = x(end-d+1:end);

% Modified Gram-Schmidt with (full) re-orthogonalization
% Computes Q amd R (required to compute he Gram basis on a oversampled grid)
Q = sym(zeros(d));
R = sym(zeros(d));
fprintf('\t QR-decomposition: ... \n');
% we start with the monomial basis {1, x, x^2, ... , x^(d-1)}
std_basis = sym(ones(d, d));
for ii = 2:d
    std_basis(:,ii) = (xr.^(ii-1)).';
end

for ii = 1:d
    Q(:,ii) = std_basis(:,ii);    
    for jj = 1:ii-1
        proj_subspace = Q(:,jj).' * Q(:,ii);        
        R(jj, ii) = proj_subspace;       
        Q(:,ii) = Q(:,ii) - proj_subspace * Q(:,jj);
    end    
    % (Full) re-orthogonalization; do the above procedure again to
    % safuguard against roundoff errors and ensure orthoganality
    for jj = 1:ii-1
        proj_subspace = Q(:,jj).' * Q(:,ii); 
        R(jj, ii) = R(jj, ii) + proj_subspace;        
        Q(:,ii) = Q(:,ii) - proj_subspace * Q(:,jj);
    end    
    % normalize for orthonormality
    p_nrm = sqrt(Q(:,ii).' * Q(:,ii));
    Q(:,ii) = Q(:,ii) / p_nrm;
    R(ii, ii) = p_nrm;
end

% oversampling
xro = sym(linspace(x(end-d+1), x(end), nos*(d-1) + 1));
std_basis_over = sym(zeros(nos * (d-1) + 1, d));
for ii = 1:d
    std_basis_over(:,ii) = sym((xro.^(ii-1))).';
end

% Q - oversampled
Qover = std_basis_over * inv(R);

g = sym(d + C + Z + E);
if( mod(g,2) == 0 )
    k = sym(-g/2:g/2-1);
else
    k = sym(-(g-1)/2:(g-1)/2);
end

Imatch = xro;% sym(linspace(x(end-dr+1), x(end), nos * (dr-1)+1));
% Iblend = sym(linspace(N*h, (N+C-1)*h, nos * (C-1)+1));
% Izero = sym(linspace(c, c + x(Z)-x(1), nos * (Z-1) + 1));
% Izero = sym(linspace(c+h-hl, c + x(Z)-x(1), nos * (Z-1) + 1));
Izero = sym(linspace(c, c + (Z-1)*h, nos * (Z-1) + 1));

y = [Imatch Izero].';
Bover = vpa(exp(2 * 1i * sym('pi') * y * k / ((g-1)*h)));

fprintf('\t SVD: ... \n');
[U, S, V] = svd(Bover, 'econ');
S_inv = diag(1./diag(S));

Iblend = sym(linspace(N*h, 1+C*h, C));
y = [Iblend].'; 
A = sym(zeros(C, d));
for jj = 1:d
    qj_over = [Qover(:,jj); zeros((Z-1)*nos +1,1)];
    D = U' * qj_over; 
    a_svd = V * S_inv * D;
    r = max( abs( double( (Bover*a_svd - qj_over ) ) ));
    fprintf( '\t%d\tresidual = %e\n', d, r );        
    A(:,jj) = real(exp(2 * 1i * sym('pi') * y * k / ((g-1)*h)) * a_svd);
end
fprintf('Precomputations: Complete! \n');

return;