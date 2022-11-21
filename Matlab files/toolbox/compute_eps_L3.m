%% Function to compute error epsilon_L2, to be used to optimize parameters
function [eps_L3] = compute_eps_L3(Hi,eps_L2)

    %% Compute first the regression of eps_L2
    
    % eps_L2 structure is
    % 1st-dimension: node
    % 2nd-dimension: time

    % Regressor structure is
    % 1st-dimension: node, or excitation
    % 2nd-dimension: time
    % 3rd-dimension: delays due to regression
    
    max_delay = size(Hi,3); % number of delays to take into account
    L = size(eps_L2,1); % number of nodes
    N = size(eps_L2,2); % number of data


    % L2 toeplitz structure        
    for node = 1:1:L
        zos = 0 * eps_L2(1,1:max_delay);
        zos(1) = eps_L2(node,1);
        L2_delay(node,:,:) = toeplitz(eps_L2(node,:)',zos);
    end

    
    %% Sum contribution of each filter over the delay coefficients
    
    % initialize first part
    eps_L3 = Hi(:,:,1) * L2_delay(:,:,1);
    for delay = 2:1:max_delay
        eps_L3 = eps_L3 + Hi(:,:,delay) * L2_delay(:,:,delay);
    end





end