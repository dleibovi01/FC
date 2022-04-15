%% This is a driver script which computes continuation data for FC methods

clear;clc;

d = 5; m = d - 1; % d is number of matching points, m is max polynomial order
C = 37;
E = 25;
Z = 12;
n_ovr = 10;
num_digits = 256;

tic;
Q = zeros(d, d);
A = zeros(C, d);
fprintf('Performing precomputations...\n');
% [Q, A] = blend_to_zero(d, C, Z, E, n_ovr, 4, num_digits);
[Q, A] = blend_to_zero(d, C, Z, E, n_ovr, num_digits, d+Z+5);
save(['FC_data_d',num2str(d),'_C_', num2str(C), '.mat'], 'Q', 'A');
toc;
