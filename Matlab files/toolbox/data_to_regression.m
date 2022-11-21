%% Function to compute regressors from the provided data

% Data structure is
% 1st-dimension: node, or excitation
% 2nd-dimension: time

% Regressor structure is
% 1st-dimension: node, or excitation
% 2nd-dimension: time
% 3rd-dimension: delays due to regression

function [wout_delay,win_delay,r_delay,eps_delay] = data_to_regression(data,netw)

    % Build regressors containing delayed data
    N = size(data.y,1);
    
    % nodes w output side
    max_delay_wout = 1 + netw.orders.D; % contains maximum delay used by regression
    
    wout_delay = zeros(netw.L,N,max_delay_wout); % pre-allocate for speed
    
    for node = 1:1:netw.L
        wout_delay(node,:,:) = toeplitz(data.y(:,node),[data.y(1,node) zeros(1,max_delay_wout-1)]);
    end
    
    
    
    % nodes w input side
    max_delay_win = 1 + netw.orders.NG; % contains maximum delay used by regression
    
    win_delay = zeros(netw.L,N,max_delay_win); % pre-allocate for speed
    
    for node = 1:1:netw.L
        win_delay(node,:,:) = toeplitz(data.u(:,node),[data.u(1,node) zeros(1,max_delay_win-1)]);
    end
    
    
    
    % external excitations
    max_delay_r = 1 + netw.orders.NR; % contains maximum delay used by regression
    
    r_delay = zeros(netw.K,N,max_delay_r); % pre-allocate for speed
    
    for excite = 1:1:netw.K
        r_delay(excite,:,:) = toeplitz(data.u(:,netw.L+excite)',[data.u(1,netw.L+excite) zeros(1,max_delay_r-1)]);
    end


    % estimated innovations
    max_delay_eps = 1 + netw.orders.NH; % contains maximum delay used by regression
    
    eps_delay = zeros(netw.L,N,max_delay_eps); % pre-allocate for speed
    
    for node = 1:1:netw.L
        eps_delay(node,:,:) = toeplitz(data.u(:,netw.L+netw.K+node)',[data.u(1,netw.L+netw.K+node) zeros(1,max_delay_eps-1)]);
    end
end