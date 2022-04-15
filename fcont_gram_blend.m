function [fx_cont_coeffs f_dp fc_l fc_r] = fcont_gram_blend(fx, d, C)

% load(['FC_data_d', num2str(d), '.mat']);
load('FC_data');
fourPts = length(fx) + C;

fr = fx((length(fx) - (d-1)):length(fx));
fl = fx(1:d);

fc_r = A*(Q.'*fr);
fc_l = flipud(A*(Q.'*flipud(fl)));
f_dp = [fx; fc_l + fc_r];

fx_cont_coeffs = fft(f_dp)/fourPts;

return
