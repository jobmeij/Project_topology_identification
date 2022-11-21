%% PEM cost function
function [J, Lambda] = PEM_cost(data,G,R,H)
    
    % Preliminaries
    L = size(data.y,2);
    N = size(data.y,1);
    transfer_w_to_eps = inv(H)*(eye(L)-G);
    transfer_r_to_eps = -inv(H)*R;
    
    % Compute prediction error for model
    eps = lsim([transfer_w_to_eps transfer_r_to_eps],[data.y data.u]);
    
    % Compute cost
    J = trace(eps'*eps) / N;
    
    % Compute covariance
    Lambda = cov(eps);

end