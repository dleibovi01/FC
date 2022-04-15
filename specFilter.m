function sigma = specFilter(eta)

% eta_c is the cutoff  wavenumber
alpha = 35; % alpha = -log(\eps_m), with \eps_m the machine epsilon
p = 5; % order of filter

sigma = exp(-alpha*abs(eta).^(2*p));
return
