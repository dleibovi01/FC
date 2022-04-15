function S = Fourier_sum(x, N, L)

S = 0;
for k = 1 : N
    S = S + 4/pi*sin(pi/L*(2*k - 1)*x)/(2*k-1);
end


end