%% Run simulation for identification with SLS

clear all;close all;clc % Start clean


%% Parameters & Network

% Simulation patameters
N = 1000; % number of datapoints 
var_e = 0.3; % power of noise
Ts = 1; % sample time

L = 5; % number of nodes
K = 5; % number of external excitations
Z = L; % Added function for rank-reduced noise
original_order = 2; % order of original network modules

addpath(genpath('../toolbox')) % access the path with the functions


% Generate network
while(1)
    [G0,R0,H0,T0,unstable] = generate_physical_network(L,K,Z,Ts,0);
    if unstable + sum(abs(eig(G0)) > 0.95) == 0 && norm(G0) > 0
        break;
    end
end

%bodemag(G0)



% Original network topology
G0_imp = impulse(G0);
original_structure.NG = squeeze(G0_imp(2,:,:)) ~= 0; 

R0_imp = impulse(R0);
original_structure.NR = squeeze(R0_imp(1,:,:)) ~= 0; 

H0_imp = impulse(H0);
original_structure.NH = squeeze(H0_imp(2,:,:)) ~= 0; 

% Fixed dynamics
% fixed.R = R0;


orders.D = original_order;
orders.NG = original_order;
orders.NR = original_order;
orders.NH = original_order;


for mc = 1:1:5 % number of Monte-Carlo data sets
    mc

    % Generate data
    [w,r,e] = generate_data(T0,N,L,Z,K,var_e);

    % make the true noise model monic, and adjust the noise source accordingly
    % for ii = 1:1:L
    %     scale = H0.num{ii,ii}(2);
    %     H0.num{ii,ii} = [H0.num{ii,ii}(2:end) 0] ./ scale;
    %     e(ii,:) = [0 e(ii,1:end-1)] .* scale;
    % end


    %% Identification with SLS
    % Network object
    %netw = idnetw('diagH',L,K,original_order,original_structure.G);
    netw = idnetw(original_structure,orders);

    % Data object
    data = iddata(w,r,Ts);

    % Option set
    opt = slsOptions();

    % Estimation
    netw_estimated = SLS(data,netw,opt);


    %% Validation
    est{mc} = netw_estimated;

end


[fit_G_fast,~] = validate_module(est,G0);
boxplot([nonzeros(fit_G_fast)])


