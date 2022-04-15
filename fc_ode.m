function [u_j eta v_x c1 c2] = fc_ode(alpha, f_xj, dom, B_l, B_r)

% FC_ODE - Routine to compute the operator S_{alpha^2 x} which approximates
% the inverse of a constant coefficient differential operator
% S_{alpha^2 x} ~ ( 1 - alpha^2 d^2/dx^2 )^{-1}
% i.e., we solve the ODE:
% -alpha^2 v''(x) + v(x) = f(x_j); v(x_l) = B_l, v(x_r) = B_r  ----- (*)
% (see (58), pp. 2023 in [1])
%   [u_j eta v_x c1 c2] = fc_ode(alpha, f_xj, dom, B_l, B_r, F) returns the
%   ODE solution v(x) and the constants associated with the two homogeneous
%   solutions c_1 and c_2, as well as finite difference correction term eta
%   and the final corrected solution u_j
%
%   Inputs:
%       alpha - (constant) coefficient of the second-order ODE term see (*)
%               above (real)
%       f_xj - value of the function f at the nodes x_j (real n-vector)
%       dom - structure containing parameters of the domain, including the
%       continuation parameters. The members of the structure are given
%       below:
%       dom
%       |-> n - number of grid points in the domain (integer)
%       |-> x_a - first grid point (real)
%       |-> x_b - last grid point (real) (x_j=x_a+j(x_b-x_a)/(n-1),j=[0,n-1])
%       |-> x_j - the equispaced grid on which ODE is solved (real n-vector)
%       |-> h - spatial discretization size (real) (Note: h = 1/(n-1))
%       |-> x_l - left boundary point (boundary values given at these pts)
%       |-> x_r - right boundary point
%       |-> cntr - structure containing continuation data and parameters
%           |-> n_delta - no of boundary section points (integer)
%           |-> n_d - no. of continuation points (integer)
%           |-> m - maximum order of Gram polynomial (integer)
%           |-> prd - period of the FC(Gram) approximation  (real)
%           |-> k - Fourier modes of FC(Gram) approximation (integer vec)
%       B_l - left boundary value of the solution (real)
%       B_r - right boundary value of the solution (real)
%       F - matrix of continuation values containing point values of the
%           Gram polynomials at the boundary sections and their
%           continuation
%           (real, matrix size (2n_d+2n_delta-4)x(2m+2))
%
%   Outputs:
%       u_j - final solution of the ODE incorporating stability corrections
%             (real, n-vector)
%       eta - finite difference correction term (see (63), pp. 2023 in [1])
%             (real, n-vector)
%       v_x - approximate solution of the ODE  (real, n-vector)
%       c_1 - constant associated with the left boundary homogeneous
%       solution (real)
%       c_2 - constant associated with the right boundary homogeneous
%       solution (real)
%
%   Ref:
%       [1] High-order unconditionally stable FC-AD solvers for general
%       smooth domains I. Basic elements
%       O. P. Bruno, M. Lyon
%       J. Comp. Phy. 229 (2010) pp. 2009-2033
%
%   Author(s):
%       Aditya Viswanathan
%       Email: aditya.v@caltech.edu
%       Last revised: Sep 23 2010
%
%       Thomas Anderson
%       Email: tanderson@caltech.edu
%       Last revised: June 16, 2015

%% Compute the Fourier continuation of f from its data values f(x_j)

% First, access data structure members so we do not need to reference
% the structure as much, making things simpler to read
n = dom.n;
x_j = dom.x_j;
h = dom.h;
x_l = dom.x_l;
x_r = dom.x_r;
d = dom.cntr.d;
C = dom.cntr.C;
m = dom.cntr.m;
prd = dom.cntr.prd;
k = dom.cntr.k;
fourPts = dom.cntr.fourPts;

% Compute the Fourier continuation coefficients of the function using the
% FC(Gram) algorithm
% denote these continuation coeffs by fx_cont_coeffs
%n_t = n + C;
[fx_cont_coeffs,~] = fcont_gram_blend(f_xj, dom);

% Use a spectral filter on f_xj
fx_cont_coeffs = fx_cont_coeffs.*specFilter(k / ((fourPts - 1)/2));


%% Compute coefficients of the inverse differential operator

% Compute the coefficients of the inverse of the constant coefficient
% differential operator
% (59), pp. 2023 in [1]

% coefficients of the inverse operator
v_hat = fx_cont_coeffs./(1+(4*pi^2*alpha^2*k.^2)/(prd^2));

% the particular solution of the ODE is the ID(F)FT of these coefficients
v_tilde_x = fourPts*real(ifft(v_hat));

%% Compute the homogeneous solutions
% define the homogeneous solutions
h1_def = @(v) exp((x_l - v)/abs(alpha));
h2_def = @(v) exp((v - x_r)/abs(alpha));

% compute these homogeneous solutions in the interior of the domain
h1 = h1_def(x_j);
h2 = h2_def(x_j);

% and their values at the boundary points
h1_xl = h1_def(x_l); h1_xr = h1_def(x_r);
h2_xl = h2_def(x_l); h2_xr = h2_def(x_r);

% compute the values of the continued function at the boundaries
% denote these by v^c_x(x_l) and v^c_x(x_r)
v_c_xl = real(exp(2*pi*1i*(x_l-x_j(1))*...
            fftshift(transpose(k))/prd)*fftshift(v_hat));
v_c_xr = real(exp(2*pi*1i*(x_r-x_j(1))*...
            fftshift(transpose(k))/prd)*fftshift(v_hat));

% Compute the constants associated with the homogeneous solutions
% Note: We solve the system Ac=b.
% A = [ h1_xl h2_xl ;...
%       h1_xr h2_xr  ];
% c = [c1 c2]^T

b = [B_l - v_c_xl; B_r - v_c_xr];
A = [ h1_xl h2_xl ;...
      h1_xr h2_xr  ];
c = A\b;
c1 = c(1); c2 = c(2);

% the homogeneous solution
homg_sol = c1*h1 + c2*h2;

% the homogeneous solution at the boundaries
homg_sol_bdry_xl = c1*h1_xl + c2*h2_xl;
homg_sol_bdry_xr = c1*h1_xr + c2*h2_xr;

%% Solution of the ODE (*)
% Add the particular and homogeneous solutions
v_x = v_tilde_x(1:n) + transpose(homg_sol);

% Also compute solution at the boundaries x_l and x_r (this is required
% later when computing additional stability-related corrections to the
% solution
v_xl = v_c_xl + homg_sol_bdry_xl;
v_xr = v_c_xr + homg_sol_bdry_xr;

%% Additional components for stability of the FC-ODE solver
% (see pp. 2023-2024 in [1])

%% Low-order finite difference corrections
% (see (63)-(64), pp. 2023 in [1])

% this is the first component of the system in (63)
eta_1 = toeplitz([1; zeros(2*d-1,1)], [1 -2 1 zeros(1,2*d-1)]);
% boundary conditions
eta_1(:,1) = []; eta_1(:,end) = [];
eta_1(d,d+1) = 0; eta_1(d+1, d) = 0;

% this is the second component of the system in (63)
eta_2 = eye(2*d);

% the complete system matrix
eta_mat = (-alpha^2/h^2)*eta_1 + eta_2;

% perform LU decomposition of the above system matrix
[L U] = lu(eta_mat);

% the FC(Gram) approximation
f_c_xj = fourPts*real(ifft(fx_cont_coeffs));
% denote the residual error in the FC(Gram) approximation, f_j - f^c(x_j)
% by rhs_res
rhs_res = transpose([f_xj(1:d) f_xj(n-d+1:end)]) - ...
                            [f_c_xj(1:d); f_c_xj(n-d+1:n)];

% solve the system (63) to obtain eta
% Ax = LUx = b; first solve Ly = b, then Ux = y
y = L\rhs_res;
eta_bdry = U\y;

% Note that the above system only gives the values of eta at the boundaries
% sections. Since eta_j = 0 elsewhere, we have
eta = [eta_bdry(1:d); ...
                zeros(dom.n-2*d, 1); ...
                        eta_bdry(d+1:end)];

%% Final discrete FC-ODE operator, incorporating the stability corrections
% (see pp. 2024 in [1])

% Step 1 - the open boundary projection of the approximate solution v(x)
v_p = open_proj(v_x, dom);

% Step 2 - the closed boundary projection of the approximate solution v(x)
v_b = closed_proj(v_x, v_xl, v_xr, dom);

% Step 3 - open boundary projection of the finite-difference solution
eta_p = open_proj(eta, dom);

% Step 4 - the final corrected solution
% (see (68), pp. 2024 in [1])
chi = min(25*alpha^2/h^2, 1);
u_j = eta - eta_p + (1-chi)*v_p + chi*v_b(2:end-1);

return

function f_p_j = open_proj(f_xj, dom)
% OPEN_PROJ - Routine to compute the open boundary projection f^p of a
% function f as defined in (65), pp. 2023 of [1]
%   f_p_j = open_proj(f_xj, n_delta, m, F) returns the open boundary
%   projection f_p_j of the function with data values f_xj at the
%   equispaced grid points x_j = j/(n-1), j = 0,1,..., n-1.
%   Boundary projections of size n_delta are computed using the precomputed
%   Gram polynomial values stored in the continuation matrix F.
%
%   Inputs:
%       f_xj - value of the function f at the nodes x_j (real n-vector)
%       n_delta - no. of boundary section points (integer)
%       m - order of Gram polynomials to use (integer)
%       F - matrix of continuation values containing point values of the
%           Gram polynomials at the boundary sections and their
%           continuation
%           (real, matrix size (2n_d+2n_delta-4)x(2m+2))
%
%   Outputs:
%       f_p_j - open boundary projection of the function f (real, n-vector)
%
%   Ref:
%       [1] High-order unconditionally stable FC-AD solvers for general
%       smooth domains I. Basic elements
%       O. P. Bruno, M. Lyon
%       J. Comp. Phy. 229 (2010) pp. 2009-2033

% First, access data structure members so we do not need to reference
% the structure as much, making things simpler to read
n = dom.n;
h = dom.h;
x_l = dom.x_l;
x_r = dom.x_r;
d = dom.cntr.d;
C = dom.cntr.C;
m = dom.cntr.m;
prd = dom.cntr.prd;
k = dom.cntr.k;

%% Gram Polynomial basis
% Read in the Gram polynomial basis used for the open boundary projection
% denote these by P_open_left and P_open_right
% Note: For equal boundary section lengths, P_open_right = P_open_left.
P_open_left = dom.cntr.Q(1:d, 1:m+1);
% P_open_right = F(dom.cntr.n_delta+dom.cntr.n_d-1:2*dom.cntr.n_delta+dom.cntr.n_d-2, 11:11+dom.cntr.m);
P_open_right = P_open_left;

%% Compute projections onto these basis
% Compute the (open) Gram polynomial projections of the function at the
% boundary sections; denote these by a_left and a_right respectively
% the boundary section function values
fx_right = f_xj(1:d);
fx_left = f_xj(n-d+1:end);
% the projections
a_left = P_open_left.'*fx_left;
a_right = P_open_right.'*fx_right;

%% Build the open boundary projection f^p
% Compute the open boundary projection as defined in (65), pp. 2023 in [1]
f_p_j = [P_open_right*a_right;  ...
            f_xj(d+1:n-d); ...
                        P_open_left*a_left];

return

function f_b_j = closed_proj(f_xj, f_xl, f_xr, dom)
% CLOSED_PROJ - Routine to compute the closed boundary projection f^b of a
% function f as defined in (65), pp. 2023 of [1]
%   f_b_j = closed_proj(f_xj, x_j, f_xl, f_xr, x_l, x_r, n_delta, n_d, m) returns
%   the closed boundary projection f_b_j of the function with data values
%   [f_xl f_xj f_xr] at the equispaced grid points
%   x_j = x_l U {j/(n-1)} U x_r, j = 0,1,..., n-1.
%   Boundary projections of size n_delta are computed using the precomputed
%
%   Inputs:
%       f_xj - value of the function f at the nodes x_j (real n-vector)
%       x_j - equispaced grid in the interior of the domain (real n-vector)
%       f_xl - vale of f at the left boundary x_l (real)
%       f_xr - vale of f at the right boundary x_r (real)
%       x_l - left boundary (real)
%       x_r - right boundary (real)
%       n_delta - no. of boundary section points (integer)
%       n_d - no. of continuation points (integer)
%       m - order of Gram polynomials to use (integer)
%
%   Outputs:
%       f_b_j - closed boundary projection of the function f (real, n-vector)
%
%   Ref:
%       [1] High-order unconditionally stable FC-AD solvers for general
%       smooth domains I. Basic elements
%       O. P. Bruno, M. Lyon
%       J. Comp. Phy. 229 (2010) pp. 2009-2033

% First, access data structure members so we do not need to reference
% the structure as much, making things simpler to read
n = dom.n;
x_j = dom.x_j;
h = dom.h;
x_l = dom.x_l;
x_r = dom.x_r;
d = dom.cntr.d;
C = dom.cntr.C;
m = dom.cntr.m;
prd = dom.cntr.prd;
k = dom.cntr.k;

%% Gram Polynomial basis
% Compute a Gram polynomial basis for the closed boundary projections
% denote these by P_closed_left and P_closed_right

% (Modified) Gram-Schmidt(MGS) on the stencil S_{left}, i.e., in the region
% [1-delta, 1] U x_r

x_left = [x_j(n-d+1:n) x_r];
P_closed_left = zeros(d+1, m+1);

% we start with the monomial basis {1, x, x^2, x^3, ..., x^m}
stand_basis = zeros(d+1, m+1);
for ind = 1:m+1
    stand_basis(:,ind) = (x_left.^(ind-1)).';
end

% MGS with (full) re-orthogonalization
for ind_a = 1:m+1
    P_closed_left(:,ind_a) = stand_basis(:,ind_a);

    % project onto space spanned by previous (ind<ind_a) basis vectors
    % use discrete norm as defined in Def. 4.2, pp. 2024, [1]
    for ind_b = 1:ind_a-1
        proj_subspace = P_closed_left(:,ind_b).'*P_closed_left(:,ind_a);

        % subtract projection
        P_closed_left(:,ind_a) = P_closed_left(:,ind_a) - ...
                            proj_subspace*P_closed_left(:,ind_b);
    end

    % (Full) re-orthogonalization; do the above procedure again to
    % safuguard against roundoff errors and ensure orthoganality
    for ind_b = 1:ind_a-1
        proj_subspace = P_closed_left(:,ind_b).'*P_closed_left(:,ind_a);

        % subtract projection
        P_closed_left(:,ind_a) = P_closed_left(:,ind_a) - ...
                                proj_subspace*P_closed_left(:,ind_b);
    end

    % normalize for orthonormality
    p_nrm = sqrt(P_closed_left(:,ind_a).'*P_closed_left(:,ind_a));
    P_closed_left(:,ind_a) = P_closed_left(:,ind_a)/p_nrm;
end

%% Create Gram polynomial basis for the right boundary section
% (Modified) Gram-Schmidt(MGS) on the stencil S_{right}, i.e.,in the region
% x_l U [1+d, 1+d+delta]

x_right = [x_l x_j(1:d)]+prd;
P_closed_right = zeros(d+1, m+1);

% we start with the standard basis {1, x, x^2, x^3, ..., x^m}
stand_basis = zeros(d+1, m+1);
for ind = 1:m+1
    stand_basis(:,ind) = (x_right.^(ind-1)).';
end

% MGS with (full) re-orthogonalization
for ind_a = 1:m+1
    P_closed_right(:,ind_a) = stand_basis(:,ind_a);

    % project onto space spanned by previous (ind<ind_a) basis vectors
    % use discrete norm as defined in [1]
    for ind_b = 1:ind_a-1
        proj_subspace = P_closed_right(:,ind_b).'*P_closed_right(:,ind_a);

        % subtract projection
        P_closed_right(:,ind_a) = P_closed_right(:,ind_a) - ...
                    proj_subspace*P_closed_right(:,ind_b);
    end

    for ind_b = 1:ind_a-1
        proj_subspace = P_closed_right(:,ind_b).'*P_closed_right(:,ind_a);

        % subtract projection
        P_closed_right(:,ind_a) = P_closed_right(:,ind_a) - ...
                            proj_subspace*P_closed_right(:,ind_b);
    end

    % normalize for orthonormality
    p_nrm = sqrt(P_closed_right(:,ind_a).'*P_closed_right(:,ind_a));
    P_closed_right(:,ind_a) = P_closed_right(:,ind_a)/p_nrm;
end

%% Compute projections onto these basis
% Compute the (closed) Gram polynomial projections of the function at the
% boundary sections; denote these by a_left and a_right respectively
% the boundary section function values
fx_right = [f_xl; f_xj(1:d)];
fx_left = [f_xj(n-d+1:end); f_xr];
% the projections
a_left = P_closed_left.'*fx_left;
a_right = P_closed_right.'*fx_right;

%% Build the closed boundary projection f^b
% Compute the closed boundary projection as defined in (66), pp. 2023 in [1]
f_b_j = [P_closed_right*a_right;  ...
            f_xj(d+1:n-d); ...
                        P_closed_left*a_left];

return
