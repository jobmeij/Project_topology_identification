%% Function to compute error epsilon_L2, to be used to optimize parameters
function [eps_L2] = compute_eps_L2(theta,netw,data)

    % Extract the current estimated model from theta
    [D,NG,NR,NH] = parameters_to_D3model(theta,netw);

    % Compensate for known dynamics
    w_minus_known = data.y -lsim(netw.fixedG,data.y) -lsim([netw.fixedR netw.fixedH],data.u);
    arx_data = iddata(w_minus_known,[data.y data.u],data.Ts);
    
    % Create data matrices with Toeplitz structures
    [wout_delay,win_delay,r_delay,eps_delay] = data_to_regression(arx_data,netw);
    
    
    
    %% Sum contribution of each filter over the delay coefficients
    
    % D
    max_delay_D = 1 + netw.orders.D; % number of delays to take into account
    for delay = 2:1:max_delay_D
        eD(:,:,delay) = D(:,:,delay) * wout_delay(:,:,delay);
    end
    eD(:,:,1) = wout_delay(:,:,1); % This is a manual addition of the feedthrough of the D model, because the D does not contain this term
    
    % N_G
    max_delay_NG = 1 + netw.orders.NG; % number of delays to take into account
    for delay = 1:1:max_delay_NG
        eNG(:,:,delay) = NG(:,:,delay) * win_delay(:,:,delay);
    end
    
    % N_R
    max_delay_NR = 1 + netw.orders.NR; % number of delays to take into account
    for delay = 1:1:max_delay_NR
        eNR(:,:,delay) = NR(:,:,delay) * r_delay(:,:,delay);
    end
    
    % N_H
    max_delay_NH = 1 + netw.orders.NH; % number of delays to take into account
    for delay = 2:1:max_delay_NH 
        eNH(:,:,delay) = NH(:,:,delay) * eps_delay(:,:,delay);
    end
    eNH(:,:,1) = eps_delay(:,:,1); % This is a manual addition of the feedthrough of the noise model, because the NH does not contain this term
    
    % Sum all contributions
    eps_L2 = sum(eD,3) - sum(eNG,3) - sum(eNR,3) - sum(eNH,3);
    
end