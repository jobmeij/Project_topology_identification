%% Function to check the AIC criterion 
function AICc = check_AICc(eps,M)
    
    N = size(eps.y,1); % number of data
    
    V = trace(eps.y' * eps.y) / N; % PEM cost for the ARX model
    
    n = size(M.Report.Parameters.Free,1); % number of free parameters
    
    AIC = 0.5 * log(V) + n / N; % AIC criterion
    AICc = AIC + (2*n^2 + 2*n) / (N - n - 1); % AICc correction for small sample sizes

end