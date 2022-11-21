%% Splitting nodes test network
% Goal: generate network data using the 3 node system, then identify it
% with splitting nodes.
%
% Generated by Job Meijer with aid of Shengling Shi, January 2020
%
clear all; close all; clc


%% Settings
N = 1000;           % Data length (samples)
Ts = 1;             % Sample time
var_e = 0.1;        % Variance of noise signals e1 and e2
var_r = 1;          % Variance of reference signal r1

% Full random G or preset G system
randomG = false;
plotEstG = false;

% Plotting
plotSignals = false;

% ARMAX fitting
removePoorFits = false;
fitPercentThreshold = 20;

% Converting dynamics to topology identification threshold
Treshold2norm = 5;
plotPositiveRates = false;

% Determine if G is identified using a MISO approach
MISO = true;

% Use bayesian algorithm to compare results
bayesianCheck = false;

% Iterations
nIterations = 100;
% TPRvect = zeros(nIterations,1);
% FPRvect = TPRvect;

% Start for loop
for i = 1:nIterations
    %% Determined network structure of 3 node system
    % Creating needed symbols
    syms w1 w2 w3 'real';
    syms e1 e2 'real';
    syms r1 'real';
    syms G12 G13 G23 G21 G31 G32 'real';
    syms H11 H22 H21 H13 'real';
    
    % Creating needed dynamic network representation vectors and matrices
    w = [w1 w2 w3]';
    e = [e1 e2]';
    r = [r1 0 0]';
    G = [0 0 0; G21 0 0; G31 0 0];
    H = [H11 H21; 0 H22; H13 0];
    R = zeros(3,3);
    R(1,1) = 1;
    
    % Combine into dynamic network representation
    w_ = G*w + H*e + R*r;
    
    % Determined transfer functions of G
    z = tf('z',Ts);
    G12 = 0;
    G13 = 0;
    G21 = (-0.2264*z - 0.1613)/(z^2 - 0.9692*z + 0.3673);
    G31 = (-0.2264*z - 0.1613)/(z^2 - 0.9692*z + 0.3673);
    G32 = 0;
    G23 = 0;
    % Combine
    G = [0, G12, G13;
        G21 0 G23;
        G31 G32 0];
    
    % Determine transfer functions of H
    H11 = 1/(z^2 - 0.9692*z + 0.3673);
    H22 = 1;
    % H21 must be strictly proper, otherwise there is correlation between H and
    % G resulting in a not identifiable network!
    H21 = 1/(z^2 - 0.9692*z + 0.3673);
    H13 = 1;
    H = [H11 H21; 0 H22; H13 0];
    
    % Set original G topology (adjust manually!)
    G_original = [0 0 0; 1 0 0; 1 0 0];
    
    %% Generate a random network (copied, adjust where needed)
    if randomG
        L = 3;      % Number of nodes
        K = 3;      % Number of external excitations
        Z = 2;      % Number of noise sources (columns of H)
        % Generate network
        while(1)
            [G0,R0,H0,T0,unstable] = generate_physical_network(L,K,Z,Ts,1);
            if unstable + sum(abs(eig(G0)) > 0.95) == 0 && norm(G0) > 0
                break;
            end
        end
        
        % Overwrite preset G
        G12 = G0(2,1);
        G13 = G0(3,1);
        G21 = G0(1,2);
        G23 = G0(3,2);
        G31 = G0(1,3);
        G32 = G0(2,3);
    end
    
    %% Create measurement data for determined network
    % Create noise data (NOT available as measurement for identification)
    e1 = sqrt(var_e).*randn(N,1);
    e2 = sqrt(var_e).*randn(N,1);
    
    % Create input signals (available as measurement for identification)
    r1 = sqrt(var_r).*randn(N,1);
    
    % Generate noise signals (NOT available as measurement for identification)
    v1 = lsim(H11,e1) + lsim(H21,e2);
    v2 = e2;
    v3 = e1;
    
    % Excite system with generated data (running simulation)
    if randomG
        w1 = r1 + v1 + lsim(G12,w2) + lsim(G13,w3);
        w2 = lsim(G21,w1) + lsim(G23,w3) + v2;
        w3 = lsim(G31,w1) + lsim(G32,w2) + v3;
    else
        w1 = r1 + v1;
        w2 = lsim(G21,w1) + v2;
        w3 = lsim(G31,w1) + v3;
    end
    
    % Plot simulation inputs and outputs
    if plotSignals
        figure()
        subplot 211
        hold on
        plot(e1);
        plot(e2);
        plot(r1);
        hold off
        grid on
        title('Simulation inputs')
        xlabel('Sample [n]')
        ylabel('Value [-]')
        legend('e1','e2','r1','Location','Best')
        subplot 212
        hold on
        plot(w1);
        plot(w2);
        plot(w3);
        hold off
        grid on
        xlabel('Sample [n]')
        ylabel('Value [-]')
        title('Simulation output')
        legend('w1','w2','w3','Location','Best')
    end
    
    %% Use generated data to identify a model using split node theory - MISO approach
    opt = armaxOptions;
    opt.Display = 'off';
    
    % Goal: estimate all G's and v's (H*e)
    % w1 = r1 + v1 + v2
    % w2 = G21*w1 + v3
    % w3 = G31*w1 + v4
    
    % Estimate transfers between nodes using SISO approach
    G12_id = iddata(w1,w2,Ts);
    G13_id = iddata(w1,w3,Ts);
    G21_id = iddata(w2,w1,Ts);
    G23_id = iddata(w2,w3,Ts);
    G31_id = iddata(w3,w1,Ts);
    G32_id = iddata(w3,w2,Ts);
    %
    G12_est = armax(G12_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    G13_est = armax(G13_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    G21_est = armax(G21_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    G23_est = armax(G23_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    G31_est = armax(G31_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    G32_est = armax(G32_id,'na',4,'nb',[4],'nc',2,'nk',[1],opt)
    
    % Estinate transfers between nodes using MISO approach
    G21_id_M = iddata(w2,[w1 w3],Ts);
    G31_id_M = iddata(w3,[w1 w2],Ts);
    %
    
    % estimate w2
    G21_est_M = armax(G21_id_M,'na',4,'nb',[4 4],'nc',2,'nk',[1 1],opt);
    G23_est_M = G21_est_M(1,2);
    G21_est_M = G21_est_M(1,1);
    
    % estimate noise v2: H22*e2
    e2_hat = w2 - lsim(G21_est_M,w1) - lsim(G23_est_M,w3);
    
    % Estimate w3
    G31_est_M = armax(G31_id_M,'na',4,'nb',[4 4],'nc',2,'nk',[1 1],opt);
    G32_est_M = G31_est_M(1,2);
    G31_est_M = G31_est_M(1,1);
    
    % estimate noise v3
    e1_hat = w3 - lsim(G31_est_M,w1) - lsim(G32_est_M,w2);
    
    
    % Left with non diagonal part
    w1x = w1 - r1 - lsim(H11,e1_hat) - lsim(H21,e2_hat);
    G12_id_M = iddata(w1x,[w2 w3],Ts);
    
    % G12_est_M = armax(G12_id_M,'na',4,'nb',[4 4],'nc',2,'nk',[1 1],opt);
    G12_est_M = armax(G12_id_M,'na',4,'nb',[4 4],'nc',0,'nk',[1 1],opt);
    G13_est_M = G12_est_M(1,2);     % Copy
    % v2_H_est = G12_est_M(1,3);
    % v3_H_est = G12_est_M(1,4);
    G12_est_M = G12_est_M(1,1);     % Overwrite
    
    % Bayesian
    if bayesianCheck
        addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\3. BSalgorithmPublic')
        W = [w1x w2 w3];
        order = 20;
        [Graph ,Score] = ForwardBackwardMIMO( W, order)
    end
    
    
    % Determined transfer function SISO approach
    G_est = [0, G12_est, G13_est;
        G21_est, 0, G23_est;
        G31_est, G32_est, 0];
    
    % Determined transfer function MISO approach
    G_est_M = [0, G12_est_M, G13_est_M;
        G21_est_M, 0, G23_est_M;
        G31_est_M, G32_est_M, 0];
    
    % Select MISO network if enabled
    if MISO
        G_est = G_est_M;
    end
    
    %% Estimating noise sources v1, v2, v3 and v4
    % v3_hat: w2 - G21*w1
    v3_hat = w2 - lsim(G21_est,w1);
    
    % v2_hat assumed to be the same as v3_hat, since it is originating from the
    % same noise source e2
    v2_hat = v3_hat;
    
    % v1_hat: w1 - r1 - v2_hat
    v1_hat = w1 - r1 - v2_hat;
    
    % v4_hat: w3 - G31*w1
    v4_hat = w3 - lsim(G31_est,w1);
    
    % Some more tests, check if still needed
    testdata = iddata(w1,[v1_hat v2_hat r1])
    testfit = armax(testdata,'na',4,'nb',[4 4 4],'nc',2,'nk',[1 1 1])
    
    r1_hat = w1 - r1;
    testdata2 = iddata(r1_hat, v2_hat,Ts);
    testfit2 = armax(testdata2,'na',4,'nb',[4],'nc',2,'nk',[1])
    
    w3_hat = iddata(w3,w1)
    testfit3 = armax(w3_hat,'na',4,'nb',[4],'nc',2,'nk',[1])
    
    %% Thresholding
    [m,n] = size(G);
    estimatedTopology = zeros(m,n);
    
    for i = 1:m
        for j = 1:n
            [vect2,~,~] = bode(G_est(i,j));
            
            Vect2norm = vecnorm(vect2)
            
            if (Vect2norm > Treshold2norm)
                estimatedTopology(i,j) = 1;
            end
        end
    end
    
    estimatedTopology = estimatedTopology
    
    % Determine TPR and FPR ratios
    L = m;
    graphvector = reshape(estimatedTopology,[L^2,1]);
    truth = reshape(G_original,[L^2,1]);
    CP = classperf(truth, graphvector,'Positive', 1, 'Negative', 0);
    TP = CP.DiagnosticTable(1,1);
    P = sum(CP.DiagnosticTable(:,1));
    FP = CP.DiagnosticTable(1,2);
    Ne = sum(CP.DiagnosticTable(:,2));
    TPR = TP/P;
    FPR = FP/Ne;
    
    % Plot TPR and FPR
    if plotPositiveRates
        figure(5)
        % set(gcf,'Position',[100 100 900 500])
        % subplot(2,2,[1 3])
        % hold all
        plot(FPR,TPR,'blackx')
        axis([0 1 0 1])
        ylabel('True Positive Rate')
        xlabel('False Rositive Rate')
        grid on
        title(['ROC curve - 3 node example network'])
    end
    %% Create boxplots of results over different realizations
    % TODO
    
    
    %% Plot transfer functions G
    if plotEstG
        figure()
        bodemag(G_est)
        hold on
        bodemag(G)
        hold off
        grid on
    end
    
    % Only plot G21 and G31
    if false
        figure()
        hold on
        bodemag(G21);
        bodemag(G21_est);
        hold off
        grid on
        legend('G21 original','G21 estimated')
        
        figure()
        hold on
        bodemag(G31);
        bodemag(G31_est);
        hold off
        grid on
        legend('G31 original','G31 estimated')
    end
    
    TPRvect(i) = TPR;
    FPRvect(i) = FPR;
end
