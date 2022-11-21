%% Function to compute 3-D matrices from parameter vector theta

% The structure of theta is [theta_D; theta_NG; theta_NR; theta_NH]

% The structure of the 3-D matrices is:
% 1st-dimension=outputs;
% 2nd-dimension=inputs; 
% 3rd-dimension=time delay

function [D,NG,NR,NH] = parameters_to_D3model(theta,netw)
    
    
    used_parameters = 0; % counter for tracking how many from parameters from the theta vector have been entered into the model so far
    
    
    

    % D polynomial
    for ii = 1:1:netw.L
        wide = size(netw.D{ii},2);
        D(ii,ii,1:wide) = [0; theta(used_parameters+(1:wide-1))];
        used_parameters = used_parameters + wide-1;
    end
    %% BEWARE!! The D polynomial did not receive 1 on the diagonal as is required for identification, this is fixed in the function where D is used %%

    % N_G polynomial    
    for colcol = 1:1:netw.L
        for rowrow = 1:1:netw.L
            wide = size(netw.NG{colcol,rowrow},2);
            if wide > 0
                NG(colcol,rowrow,1:wide) = [0; theta(used_parameters+(1:wide-1))];
                used_parameters = used_parameters + wide-1;
            end
        end
    end

    % The code above may result in an NG of wrong dimension, so correct it
    max_delay_NG = 1 + netw.orders.NG;
    if ~exist('NG')
        NG = zeros(netw.L,netw.L,max_delay_NG);% when it was empty, then make it
    elseif size(NG,1) < netw.L || size(NG,2) < netw.L
        NG(netw.L,netw.L,:) = zeros(1,1,max_delay_NG);% when the dimension was too small, increase size
    end

    % N_R polynomial
    for colcol = 1:1:netw.L
        for rowrow = 1:1:netw.K
            wide = size(netw.NR{colcol,rowrow},2);
            if wide > 0
                NR(colcol,rowrow,1:wide) = theta(used_parameters+(1:wide));
                used_parameters = used_parameters + wide;
            end
        end
    end

    % The code above may result in an NR of wrong dimension, so correct it
    max_delay_NR = 1 + netw.orders.NR;
    if ~exist('NR')
        NR = zeros(netw.L,netw.K,max_delay_NR);% when it was empty
    elseif size(NR,1) < netw.L || size(NR,2) < netw.K
        NR(netw.L,netw.K,:) = zeros(1,1,max_delay_NR);% when the dimension was too small
    end
    
    % N_H polynomial
    for colcol = 1:1:netw.L
        for rowrow = 1:1:netw.L
            wide = size(netw.NH{colcol,rowrow},2);
            if wide > 0
                NH(colcol,rowrow,1:wide) = [0; theta(used_parameters+(1:wide-1))];
                used_parameters = used_parameters + wide-1;
            end
        end
    end

    % The code above may result in an NH of wrong dimension, so correct it
    max_delay_NH = 1 + netw.orders.NH;
    if ~exist('NH')
        NH = zeros(netw.L,netw.L,max_delay_NH); % when it was empty
    elseif size(NH,1) < netw.L || size(NH,2) < netw.L
        NH(netw.L,netw.L,:) = zeros(1,1,max_delay_NH); % when the dimension was too small
    end
    %% BEWARE!! The NH polynomial did not receive 1 on the diagonal as is required for identification, this is fixed in the function where NH is used %%
end