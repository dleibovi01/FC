%
% function dummy = blend_to_zero( num_original_points, ...
%     num_continuation_points, num_zero_points, num_extra_points, ...
%     oversampling_factor, modes_to_reduce, num_digits )
%
% Creates a blending-to-zero for use in the accelerated FC method.
% Continuations are computed so that the first FP points are taken from the
% Gram polynomial basis.  The next CP points are continuation points used
% to form the continuation.  ZP points follow where the continuation
% function is forced to zero.  Finally, EP points follow this which will
% allow for a smooth, periodic function.  Note: The EP points are never
% actually used in the continuation.  All computations are performed at OV
% times as many point values (the oversampling part of the FC method).
%
% num_original_points     : The number of function values to sample.
% num_continuation_points : The number of continuation points to add.
% num_zero_points         : The number of points in the "zero" region.
% num_extra_points        : The number of points following the zero region.
% oversampling_factor     : The oversampling factor.
% modes_to_reduce         : The number of Fourier modes to subtract from
%                           the least squares problem.
% num_digits              : The number of digits to compute with.
%

function [Q,A] = blend_to_zero( num_original_points, ...
    num_continuation_points, num_zero_points, num_extra_points, ...
    oversampling_factor, modes_to_reduce, num_digits )

fprintf( '------------------------------------------------------------\n' );
fprintf( '%d POINT CONTINUATION\n\n', num_original_points );
start_time = tic;

fprintf( 'Setting up... '); tic;

digits(num_digits);

% the total number of points in the coarse domain
coarse_N = num_original_points + num_continuation_points + ...
    num_zero_points + num_extra_points;

% the bandwidth
bandwidth = floor(coarse_N/2)-modes_to_reduce;

% the total number of oversampled points in the sampled interval
num_oversampled_original = oversampling_factor*(num_original_points-1) + 1;

% the total number of oversampled points in the zero interval
num_oversampled_zero = oversampling_factor*(num_zero_points-1) + 1;

% the total number of points in the oversampled domain
oversampled_N = coarse_N*oversampling_factor;

% the coarse grid for the long interval [0,1)
coarse_full_grid = sym( (0:coarse_N-1).' )/coarse_N;

% the coarse grid for the sampled interval
coarse_short_grid = sym( (0:num_original_points-1).' )/coarse_N ;

% the fine grid for the sampled interval
fine_short_grid = sym( (0:num_oversampled_original-1).' )/oversampled_N;

% the 'midpoint' where zeroing starts
midpoint = (num_original_points + num_continuation_points)/coarse_N;

% the coarse grid for the zero interval
coarse_zero_grid = sym( coarse_N*midpoint + ...
    (0:num_zero_points-1).' )/coarse_N;

% the fine grid for the zero interval
fine_zero_grid = sym( oversampled_N*midpoint + ...
    (0:num_oversampled_zero-1).' )/oversampled_N;

% the Fourier modes
k  = (-bandwidth:bandwidth).';

% the number of Fourier modes
num_modes = 2*bandwidth+1;

fprintf( 'Done. (%1.3fs)\n', toc );

%% generate the polynomial basis

fprintf( 'Generating polynomial basis... '); tic;

% these variables will hold the polynomial basis on the coarse and fine
% grids respectively
coarse_p  = sym( zeros(num_original_points,num_original_points) );
fine_p    = sym( zeros(num_oversampled_original,num_original_points) );

% generate monomials
for d = 0:num_original_points-1
    coarse_p(:,d+1)  = coarse_short_grid.^d;
    fine_p(:,d+1)    = fine_short_grid.^d;
end

% use QR to construct coarse basis
[coarse_Q,coarse_R] = sym_qr(coarse_p);

% use R to construct fine basis
fine_Q = fine_p*inv(coarse_R);

fprintf( 'Done. (%1.3fs)\n', toc );

%% generate the continuations

fprintf( 'Generating reconstruction matrices... '); tic;

% generate reconstruction matrix on fine mesh
A1 = [ vpa( cos( 2*sym('pi')*fine_short_grid*k.') );
        vpa( cos( 2*sym('pi')*fine_zero_grid*k.' ) ) ];
A2 = [ vpa( sin( 2*sym('pi')*fine_short_grid*k.') );
        vpa( sin( 2*sym('pi')*fine_zero_grid*k.' ) ) ];
A  = [A1,A2];

% generate reconstruction matrix on total coarse mesh
C1 = vpa( cos( 2*sym('pi')*coarse_full_grid*k.' ) );
C2 = vpa( sin( 2*sym('pi')*coarse_full_grid*k.' ) );
C  = [C1,C2];

fprintf( 'Done. (%1.3fs)\n', toc );

fprintf( 'Performing SVD... '); tic;

% invert the SVD
[U,S,V] = svd(A);
U = U(:,1:size(S,1));
s = diag(S);
z = double( s(1:end-1)./s(2:end) );
ind = min( [find( z > 1e16 ).', num_modes] );
t = 1./s(1:ind);
V = V(:,1:ind);
U = U(:,1:ind);

fprintf( 'Done. (%1.3fs)\n', toc );

fprintf( 'Computing continuations...\n'); tic;

% create the continuations
cont_data = zeros(coarse_N,num_original_points);
for d = 1:num_original_points
    b = [fine_Q(:,d);zeros(num_oversampled_zero,1)];
    a = V*(t.*(U'*b));
    r = max( abs( double( A*a - b ) ) );
    fprintf( '\t%d\tresidual = %e\n', d, r );

    cont_data(:,d) = double( C*a(1:2*num_modes) );
end

fprintf( 'Done. (%1.3fs)\n', toc );

Q = real( double( coarse_Q ) );
fprintf( 'Done. Error = %1.3e. (%1.3fs)\n', ...
    max( max( real( double( abs( Q'*Q - eye(num_original_points) ) ) ) ) ), toc );

clear A;
A = cont_data((num_original_points+1):num_original_points+num_continuation_points, ...
              1:num_original_points);
size(A)

fprintf( 'Saving continuations... '); tic;
fprintf( '\n\nDONE. (%1.3fs)\n', toc(start_time) );
fprintf( '------------------------------------------------------------\n\n' );
