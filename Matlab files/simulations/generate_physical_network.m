%% Generate a dynamic network based on a mass-spring-damper system
% Nodes are the position of each mass
% External inputs are forces added externally

function [Gd,Rd,Hd,Td,unstable] = generate_physical_network(L,K,Z,Ts,Reye)
    
    % Parameters
    pd = 0.3; % chance of nodes connected by damper
    ps = 0.3; % chance of nodes connected by spring

    %% Generate random values
    % Masses
    M = diag(10.^random('Uniform',-1,1,L,1)); % Random mass between 0.1 and 10 

    % Dampers
    damper_connections = tril(binornd(ones(L),pd*ones(L)),-1); % connections between masses
    damper_initial = 10.^random('Uniform',-1,0,L,L); % Random damper between 0.1 and 1 
    weights_d = damper_connections.*damper_initial + damper_connections'.*damper_initial'; % Weights of damper graph
    dissipation_d = 10.^random('Uniform',-2,-1,L,1); % Dissipation of each mass as damper between 0.01 and 0.1 
    damper_diag = dissipation_d + sum(weights_d,2); % diagonal of D matrix
    D = diag(damper_diag) - weights_d; % Damper matrix

    % Springs
    spring_connections = tril(binornd(ones(L),ps*ones(L)),-1); % connections between masses
    spring_initial = 10.^random('Uniform',-1,1,L,L); % Random spring between 0.1 and 10
    weights_s = spring_connections.*spring_initial + spring_connections'.*spring_initial'; % Weights of spring graph
    dissipation_s = 10.^random('Uniform',-2,-1,L,1); % Dissipation of each mass as spring between 0.01 and 0.1 
    spring_diag = dissipation_s + sum(weights_s,2); % diagonal of S matrix
    S = diag(spring_diag) - weights_s; % Spring matrix

    % Input matrix
    F = eye(L,K);

    % Noise matrix
    C = eye(L);
    
    if (L ~= Z)
        C = C(:,1:Z);   % Removing some columns to make H0 rank-reduced 
        C(1,Z) = 1;     % Making H positive definite
    end

    %% Make transfer functions
    s = tf('s');

    Q = inv(s^2 * diag(diag(M)) + s * diag(diag(D)) + diag(diag(S))); % denominators
    G2 = s * (D - diag(diag(D))) + (S - diag(diag(S))); % numerators of modules

    Gc = Q * G2; % modules
    Rc = Q * F; % external modules
    Hc = Q * C; % noise modules 

    % Discrete time system
    Gd = c2d(Gc,Ts,'zoh');
    
    Hd = c2d(Hc,Ts,'zoh');
    
    if Reye % this toggles whether the R0 is identity or contains dynamics
        Rd = tf(ss([],[],[],F,Ts));
    else
        Rd = c2d(Rc,Ts,'zoh');
    end
    
    Td = (eye(L)-Gd) \ [Rd Hd];

    %% Tests
    % The network may become unstable due to discretization...
    unstable = sum(abs(eig(Td)) >= 1);
end


