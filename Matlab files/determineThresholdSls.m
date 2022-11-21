%% Comparing different MIMO identification algorithms with rank reduced noise - determining threshold SLS
clear all; close all; clc
%
nTrials = 50; % was 50
trial = 1;
counter = 1;
num_loop = 1000;
%
thresholdVector = [1 5 10 25 50 100];

for threshold = thresholdVector
    disp(['>>> Running threshold value ',num2str(threshold),' <<<'])
    counter = 1;
    
    disp(['Starting loop, counter is ', num2str(counter)])
    while counter <= num_loop
        disp(['Running loop, counter is ',num2str(counter)])
        
        try     % Used to catch and skip error in case one occurs
            %% restart this section in case of error!
            for trial = trial:nTrials
                disp(['Running iteration ',num2str(trial),' of ',num2str(nTrials),', so far so good...'])
                
                %% Generate a MIMO model
                tic
                % Simulation patameters
                N = 1000;               % number of datapoints
                var_e = 1; % 3;          % power of noise
                Ts = 1;                 % sample time
                L = 5;                  % number of nodes
                K = 5;                  % number of external excitations (disabled since r is always set to zero)
                Z = 5;                  % number of columns for H, must be =< L
                original_order = 2;     % order of original network modules
                
                disp(['Going to create a network with L=',num2str(L),', K=',num2str(K),', N=',num2str(N),'...'])
                
                % Adding Harm's toolbox to create the network
                addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\2. Harm rank reduced noise algorithm\toolbox') % access the path with the functions
                addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\2. Harm rank reduced noise algorithm\simulations')
                
                networkAccepted = false;
                while(networkAccepted == 0)
                    
                    % Generate network
                    while(1)
                        [G0,R0,H0,T0,unstable] = generate_physical_network(L,K,Z,Ts,1);
                        if unstable + sum(abs(eig(G0)) > 0.95) == 0 && norm(G0) > 0
                            break;
                        end
                    end
                    
                    % Original network topology
                    G0_imp = impulse(G0);
                    original_structure.NG = squeeze(G0_imp(2,:,:)) ~= 0;
                    
                    R0_imp = impulse(R0);
                    original_structure.NR = squeeze(R0_imp(1,:,:)) ~= 0;
                    
                    H0_imp = impulse(H0);
                    original_structure.NH = squeeze(H0_imp(2,:,:)) ~= 0;
                    
                    % Fixed dynamics
                    fixed.R = R0;
                    
                    orders.D = original_order;
                    orders.NG = original_order;
                    orders.NR = original_order;
                    orders.NH = original_order;
                    
                    % Generate data
                    [w,r,e] = generate_data(T0,N,L,Z,K,var_e);        % r is set to zeros
                    
                    % make the true noise model monic, and adjust the noise source accordingly
                    for ii = 1:1:Z
                        scale = H0.num{ii,ii}(2);
                        H0.num{ii,ii} = [H0.num{ii,ii}(2:end) 0] ./ scale;
                        e(ii,:) = [0 e(ii,1:end-1)] .* scale;
                    end
                    
                    % Check if there are nodes with absolutely no signals (all zeros) and fill
                    % it with a random signal
                    networkAccepted = true;
                    for i = 1:L
                        if nnz(w(:,i)) == 0
                            warning(['Filling node ',num2str(i),' with new random noise since the data consisted of all zeros!'])
                            networkAccepted = false;
                        end
                    end
                end
                
                disp(['Done creating network and data, this took ',num2str(toc),' seconds.'])
                
                % Overview of parameters:
                % G: connection between nodes
                % H: noise entering nodes
                % R: input signals
                % w: measured data from nodes
                % e: generated noise inserted through H
                % r: generated noise inserted through R
                
                %     figure(2)
                %     plot(w)
                %     grid on
                %     xlabel('Sample [n]')
                %     ylabel('Amplitude [-]')
                %     title('generated data from nodes')
                
                %% Estimate model from data with Granger causality algorithm
                addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\7. Granger toolbox\mvgc_v1.0')
                startup;    % Start toolbox
                
                disp('Starting (1) Granger causality algorithm...')
                momax = 100;        % model order for estimation
                icregmode = 'LWR';
                verb = 'true';
                regmode = 'OLS';    % TODO check what this does
                acmaxlags = 2000;   % TODO check what this does
                [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(w',momax,icregmode,verb);      % does not work with moAIC and moBIC since w is not multi trial?
                
                amo = 1;    % Only one order?
                
                % Select model order.
                morder = 'AIC'; % initial setting
                if     strcmpi(morder,'actual')
                    morder = amo;
                    fprintf('\nusing actual model order = %d\n',morder);
                elseif strcmpi(morder,'AIC')
                    morder = moAIC;
                    fprintf('\nusing AIC best model order = %d\n',morder);
                elseif strcmpi(morder,'BIC')
                    morder = moBIC;
                    fprintf('\nusing BIC best model order = %d\n',morder);
                else
                    fprintf('\nusing specified model order = %d\n',morder);
                end
                
                % Estimate VAR model of selected order from data.
                [A,SIG] = tsdata_to_var(w',morder,regmode);
                
                % Check for failed regression
                assert(~isbad(A),'VAR estimation failed');
                
                % NOTE: at this point we have a model and are finished with the data! - all
                % subsequent calculations work from the estimated VAR parameters A and SIG.
                
                % Autocovariance calculation
                [G,info] = var_to_autocov(A,SIG,acmaxlags);
                
                % Plot results
                var_info(info,false); % report results (and bail out on error)
                
                % Determine the estimated network 1's and 0's from Granger algorithm
                if (info.error == 0)
                    F = autocov_to_pwcgc(G);
                    nobs = N;       % Number of observations
                    ntrials = 1;   % Copied from demo
                    nvars = L;
                    tstat = '';     % Copied from demo
                    pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
                    alpha = 0.05;   % Copied from demo
                    mhtc = 'FDR';   % Copied from demo
                    GrangerEstimate  = significance(pval,alpha,mhtc);
                    for i = 1:L
                        GrangerEstimate(i,i) = 0;
                    end
                else
                    % In case of an error, set estimated model to all zeros
                    warning('Granger estimate not complete, setting it to all zeros!')
                    GrangerEstimate = zeros(L,L);
                end
                
                % Determine TPR and FPR for Granger algorithm
                G_graphvector = reshape(GrangerEstimate,[L^2,1]);
                G_truth = reshape(original_structure.NG,[L^2,1]);
                G_CP = classperf(G_truth, G_graphvector,'Positive', 1, 'Negative', 0);
                G_TP = G_CP.DiagnosticTable(1,1);
                G_P = sum(G_CP.DiagnosticTable(:,1));
                G_FP = G_CP.DiagnosticTable(1,2);
                G_Ne = sum(G_CP.DiagnosticTable(:,2));
                G_TPR = G_TP/G_P;
                G_FPR = G_FP/G_Ne;
                
                disp('Done running (1) Granger causality algorithm.')
                
                %% Estimate model from data with Bayesan algorithm
                disp('Starting (2) Bayesian algorithm identification...')
                tic
                bayesian_order = 100 %original_order;
                addpath('C:\Users\Job\Dropbox\school 2019 2020\Blok 2 Internship\1. Received files\3. BSalgorithmPublic')
                [Graph_bayesian, Score_bayesian] = ForwardBackwardMIMO(w, bayesian_order,[],[],[],[],0); % use "help ForwardBackwardMIMO" for more information
                
                % Plot TPR vs FPR for Bayesian algorithm
                B_graphvector = reshape(Graph_bayesian,[L^2,1]);
                B_truth = reshape(original_structure.NG,[L^2,1]);
                B_CP = classperf(B_truth, B_graphvector,'Positive', 1, 'Negative', 0);
                B_TP = B_CP.DiagnosticTable(1,1);
                B_P = sum(B_CP.DiagnosticTable(:,1));
                B_FP = B_CP.DiagnosticTable(1,2);
                B_Ne = sum(B_CP.DiagnosticTable(:,2));
                B_TPR = B_TP/B_P;
                B_FPR = B_FP/B_Ne;
                
                disp(['Done running (2) Bayesian algorithm, this took ', num2str(toc),' seconds.'])
                
                %% Estimate model from data with Rank-reduced process noise algorithm (SLS)
                disp('Starting (3) SLS algorithm...')
                tic
                % Network object
                %netw = idnetw('diagH',L,K,original_order,original_structure.G);
                netw = idnetw(original_structure,orders,fixed);
                
                % Data object
                data = iddata(w,r,Ts);
                
                % Option set
                opt = slsOptions();
                
                % Estimation
                SLS_netw_estimated = SLS(data,netw,opt);
                
                % Validation
                mc = 1;     % Monte carlo
                est{mc} = SLS_netw_estimated;
                
                % opt_dc = slsOptions('regularization','dc');
                % netw_dc = SLS(data,netw,opt_dc);
                %
                % opt_tc = slsOptions('regularization','tc');
                % netw_tc = SLS(data,netw,opt_tc);
                %
                % opt_ss = slsOptions('regularization','ss');
                % netw_ss = SLS(data,netw,opt_ss);
                
                [fit_G_fast, SLS_err(:,:,trial)] = validate_module(est,G0);
                
                [Gw,Gh] = size(G0);
                disp('Filling SLS matrices with zeros')
                SLSabsDiff = zeros(Gw,Gh,nTrials);
                SLSpeakDiff = SLSabsDiff;
                SLSestimate = zeros(Gw,Gh);
                
                % Compare SLS by checking the norm of every estimated transfer
                % function, if the 2-norm is above a certain treshold set it to 1.
                %             Treshold2norm = 1e-1; % was 1e-1 before for mean(vect2)
                Treshold2norm = threshold;
                omega = 0.01:0.001:pi;
                for i = 1:Gw;
                    for j = 1:Gh;
                        if (isempty(pole(SLS_netw_estimated.G(i,j))) == 0)
                            disp(['Comparing original and estimated transfer function at ',num2str(i),'x',num2str(j)])
                            [vect2,~,~] = bode(SLS_netw_estimated.G(i,j),omega);
                            
                            Vect2norm = vecnorm(vect2);
                            %                         Vect2norm = mean(vect2)
                            
                            SLSabsDiff(i,j,trial) = Vect2norm;
                            
                            if (Vect2norm > Treshold2norm)
                                SLSestimate(i,j) = 1;
                            end
                        else
                            disp(['Estimated transfer function for ',num2str(i),'x',num2str(j),' is empty, skipping...'])
                        end
                    end
                end
                
                % Plot TPR vs FPR for Bayesian algorithm
                S_graphvector = reshape(SLSestimate,[L^2,1]);
                S_truth = reshape(original_structure.NG,[L^2,1]);
                S_CP = classperf(S_truth, S_graphvector,'Positive', 1, 'Negative', 0);
                S_TP = S_CP.DiagnosticTable(1,1);
                S_P = sum(S_CP.DiagnosticTable(:,1));
                S_FP = S_CP.DiagnosticTable(1,2);
                S_Ne = sum(S_CP.DiagnosticTable(:,2));
                S_TPR = S_TP/S_P;
                S_FPR = S_FP/S_Ne;
                
                disp(['Done running SLS algorithm (3), this took ', num2str(toc),' seconds.'])
                
                %% Compare estimated models with eachother and the generated/imported model
                % save TPR and FPR values over iterations in a vector
                TPR(trial,:) = [G_TPR, B_TPR, S_TPR];
                FPR(trial,:) = [G_FPR, B_FPR, S_FPR];
                
                disp('...')
                disp('Checking if SLS algorithm model compares correctly to true model:')
                SLSmatch = zeros(L,L);
                SLSmatch = (original_structure.NG == SLSestimate);
                SlsCorrectPercent(trial) = (sum(SLSmatch(:))/numel(SLSmatch))*100;
                Sls_TPR = SlsCorrectPercent;
                Sls_FPR = 0;
                for i = 1:size(SLSestimate,1)
                    for j = 1:size(SLSestimate,2)
                        if (SLSestimate(i,j) == 1)
                            if (original_structure.NG(i,j) == 0)
                                Sls_FPR = Sls_FPR + 1;
                            end
                        end
                    end
                end
                Sls_FPR = (Sls_FPR/size(SLSestimate,1)^2) * 100;
                disp('...')
                
                %
                disp('...')
                disp('Checking if Bayesian algorithm model compares correctly to true model:')
                BayesianMatch = (original_structure.NG == Graph_bayesian);
                BayesianCorrectPercent(trial) = (sum(BayesianMatch(:))/numel(BayesianMatch))*100;
                Bayesian_TPR = BayesianCorrectPercent;
                Bayesian_FPR = 0; % Set initially to zero
                % Add one if algorithm thinks there is a link while original structure
                % has no link
                for i = 1:size(Graph_bayesian,1)
                    for j = 1:size(Graph_bayesian,2)
                        if (Graph_bayesian(i,j) == 1)
                            if (original_structure.NG(i,j) == 0)
                                Bayesian_FPR = Bayesian_FPR + 1;
                            end
                        end
                    end
                end
                Bayesian_FPR = (Bayesian_FPR/size(Graph_bayesian,1)^2) * 100;
                disp('...')
                
                disp('Comparing results using transfer functions, going to plot them now...')
                tic
                
                if false
                    figure(1)
                    hold all
                    bodemag(G0)
                    grid on
                    
                    figure(2)
                    hold all
                    bodemag(H0)
                    grid on
                end
                
                if false
                    figure(4)
                    set(gcf,'Position',[100 100 900 500])
                    subplot(2,2,[1 3])
                    hold all
                    plot(FPR(:,1),TPR(:,1),'bluex')
                    plot(FPR(:,2),TPR(:,2),'ro')
                    plot(FPR(:,3),TPR(:,3),'blackp')
                    legend('Granger algorithm','Bayesian algorithm','SLS algorithm')
                    axis([0 1 0 1])
                    ylabel('True Positive Rate')
                    xlabel('False Rositive Rate')
                    grid on
                    sgtitle(['ROC curve - ',num2str(L),' node network with rank (',num2str(Z),') reduced noise over ',num2str(nTrials),' realizations'])
                    % TPR rate
                    subplot(2,2,2)
                    hold on
                    plot(1:size(TPR,1),TPR(:,1),'bluex')
                    plot(1:size(TPR,1),TPR(:,2),'ro')
                    plot(1:size(TPR,1),TPR(:,3),'blackp')
                    hold off
                    grid on
                    xlabel('Iteration [n]')
                    ylabel('TPR rate')
                    axis([1 nTrials -inf 1])
                    % FPR rate
                    subplot(2,2,4)
                    hold on
                    plot(1:size(FPR,1),FPR(:,1),'bluex')
                    plot(1:size(FPR,1),FPR(:,2),'ro')
                    plot(1:size(FPR,1),FPR(:,3),'blackp')
                    hold off
                    grid on
                    xlabel('Iteration [n]')
                    ylabel('FPR rate')
                    axis([1 nTrials 0 inf])
                end
                
                disp(['Done plotting transfer function matrices, this took ', num2str(toc),' seconds.'])
                
                if trial == nTrials
                    counter = num_loop + 1;      % stop while loop
                end
                
            end
            
        catch err
            if strcmp(err.identifier, 'whatever identifier is thrown by GMM')
                disp('error detected')
            else
                disp('no error detected')
                if trial == nTrials
                    counter = num_loop + 1;     % stop while loop
                end
            end
        end
        counter = counter + 1;
        
    end % while end
    
    % Compute average TPR/FPR rates
    avg.avgTPRgranger = mean(TPR(:,1));
    avg.avgFPRgranger = mean(FPR(:,1));
    avg.avgTPRbayesian = mean(TPR(:,2));
    avg.avgFPRbayesian = mean(FPR(:,2));
    avg.avgTPRsls = mean(TPR(:,3));
    avg.avgFPRsls = mean(FPR(:,3));
    
    dist.avgDistGranger = mean(sqrt(FPR(:,1).^2+(ones(length(TPR(:,1)),1)-TPR(:,1)).^2));
    dist.avgDistBayesian = mean(sqrt(FPR(:,2).^2+(ones(length(TPR(:,2)),1)-TPR(:,2)).^2));
    dist.avgDistSls = mean(sqrt(FPR(:,3).^2+(ones(length(TPR(:,3)),1)-TPR(:,3)).^2));
    
    % Create bars with height TPR
    X_TPR = categorical({'Granger','Bayesian','SLS'});
    X_TPR = reordercats(X_TPR,{'Granger','Bayesian','SLS'});
    Y_TPR = [avg.avgTPRgranger avg.avgTPRbayesian avg.avgTPRsls];
    %
    if false
        figure(5)
        set(gcf,'Position',[100 100 900 500])
%         sgtitle(['Average TPR/FPR - ',num2str(L),' node network with rank (',num2str(Z),') reduced noise over ',num2str(nTrials),' realizations'])
        title(['Average TPR/FPR - ',num2str(L),' node network with rank (',num2str(Z),') reduced noise over ',num2str(nTrials),' realizations'])
        subplot 121
        hold on
        bar(X_TPR,Y_TPR)
        hold off
        ylabel('Average TPR')
        
        % FPR
        X_FPR = categorical({'Granger','Bayesian','SLS'});
        X_FPR = reordercats(X_FPR,{'Granger','Bayesian','SLS'});
        Y_FPR = [avg.avgFPRgranger avg.avgFPRbayesian avg.avgFPRsls];
        %
        subplot 122
        hold on
        bar(X_FPR,Y_FPR)
        hold off
        ylabel('Average FPR')
    end
    
    
    
    % Compute average TPR and FPR and add to figure    
    figure(14)
    hold all
    plot(avg.avgFPRsls,avg.avgTPRsls,'*')
%     hold off
    grid on
    xlabel('FPR')
    ylabel('TPR')
    sgtitle(['Average TPR/FPR - ',num2str(L),' node network with rank (',num2str(Z),') reduced noise over ',num2str(nTrials),' realizations'])
end

% legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','Location','Best')
