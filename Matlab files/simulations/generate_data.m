%% Function to generate data from a network

function [w,r_out,e] = generate_data(T,N,L,Z,K,var_e)

if K > 0
    r = randn(K,N);
else
    r = zeros(K,N);   % disabling input r by filling it with zeros
end

warning('external excitations are disabled in generate_data!')
e = var_e * randn(Z,N);

w = lsim(T,[r; e]);
r_out = r';
end