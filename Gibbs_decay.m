clear all
close all


N = 10;
S = 0;
L = 100;
x = (-50 : 0.1 : 50);
% for k = 1 : N
%     S = S + 4/pi*sin(pi/L*(2*k - 1)*x)/(2*k-1);
% end
S = Fourier_sum(x, N, L);
plot(x, S);
y = (x >=0).*x;

figure
S_plus = Fourier_sum(y, N, L);
plot(y, S_plus - 1)