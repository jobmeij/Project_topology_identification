%% Easy access to SLS algorithm by a 'frontend'


function [varargout] = SLS(varargin)
    % This function applies the SLS network identification algorithm
    %
    % Use it as
    % [G_est,R_est,H_est] = SLS(data,netw,opt)
    %
    % Where data is an IDDATA structure
    % Where netw is an IDNETW structure, (optional)
    % Where opt is generated by slsOptions(), (optional)
    
    
    
    
    
    % Define some shorthand notation variables
    data = varargin{1};
    L = size(data.y,2);
    K = size(data.u,2);
    N = size(data.y,1);
    
    if nargin > 1
        netw = varargin{2};
    else
        netw = idnetw('diagH',L,K,2); % Default model structure
    end
    
    if nargin > 2
        opt = varargin{3};
    else
        opt = slsOptions(); % Default options
    end
    
    
    ARX_orders = opt.ARX_orders;
    
    
    
    % Rank reduced noise is implemented by modifying the model structure
    if netw.p < netw.L % only execute if rank reduced noise
        netw_original = netw; % save the original network model
        
        % edit the model to equivalent form by adding unity columns
        netw.p = L; % because it is a 'fake' full rank noise model
        netw.adjacencyH(:,netw_original.p+1:netw_original.L) = [zeros(netw_original.p,netw_original.L-netw_original.p);eye(netw_original.L-netw_original.p)]; % columns are added, because must have monic square noise model
        netw.fixedH(:,netw_original.p+1:netw_original.L) = tf(ss([],[],[],double(netw.adjacencyH(:,netw_original.p+1:netw_original.L)),-1)); % the extra columns have no parameters
        netw.H(:,netw_original.p+1:netw_original.L) = netw.fixedH(:,netw_original.p+1:netw_original.L); % update the H transfer function
        
        % update the polynomial coefficients
        for colcol = 1:1:netw_original.L
            for rowrow = netw_original.p+1:1:netw_original.L
                netw.NH{colcol,rowrow} = zeros(1,0);% remain empty when there is no link
            end
        end
        
        % update the covariance matrix
        netw.Lambda(netw_original.p+1:1:netw_original.L,netw_original.p+1:1:netw_original.L) = eye(netw_original.L-netw_original.p);
    end


    %% Quality choice is made here
    switch opt.quality
        case 'best'
            % The best ARX order must be found, estimate for many orders
            for ii = 1:1:length(ARX_orders)
                n_arx = ARX_orders(ii);
                
 
                [G_step2{ii},R_step2{ii},H_step2{ii},M2{ii},G_step3{ii},R_step3{ii},H_step3{ii},M3{ii},epsi{ii}] = SLS_core(data,netw,n_arx,opt);
            end
            
            % Compute performance of all estimates based on PEM criterion
            for ii = 1:1:length(ARX_orders)
                % cost of step 2: consistency for each ARX order
                [J_PEM(1,ii),~] = PEM_cost(data,G_step2{ii},R_step2{ii},H_step2{ii});
                
                % cost of step 3: efficiency for each ARX order and iteration
                for iter = 1:1:size(G_step3{ii},2)
                    [J_PEM(1+iter,ii),~] = PEM_cost(data,G_step3{ii}{iter},R_step3{ii}{iter},H_step3{ii}{iter});
                end
                J_PEM(1+size(G_step3{ii},2)+1:opt.max_iter+1,ii) = inf; % Iterations that are skipped should not be selected, therefore cost is infinite
            end
            
            % Select best estimate, both: ARX order & iteration
            [A,B] = min(J_PEM,[],1);
            [~,opt_ARX] = min(A);
            opt_iter = B(opt_ARX)-1; 
            
            
            % Output
            if opt_iter == 0
                netw.G = G_step2{opt_ARX};
                netw.R = R_step2{opt_ARX};
                netw.H = H_step2{opt_ARX};
                netw.D = M2{opt_ARX}.D;
                netw.NG = M2{opt_ARX}.NG;
                netw.NR = M2{opt_ARX}.NR;
                netw.NH = M2{opt_ARX}.NH;
            else
                netw.G = G_step3{opt_ARX}{opt_iter};
                netw.R = R_step3{opt_ARX}{opt_iter};
                netw.H = H_step3{opt_ARX}{opt_iter};    
                netw.D = M3{opt_ARX}{opt_iter}.D;
                netw.NG = M3{opt_ARX}{opt_iter}.NG;
                netw.NR = M3{opt_ARX}{opt_iter}.NR;
                netw.NH = M3{opt_ARX}{opt_iter}.NH;         
            end
            
            
            
        case 'fast'
            % A good ARX order must be found, estimate for many orders
            for ii = 1:1:length(ARX_orders)
                n_arx = ARX_orders(ii);  
                
                optio = opt;
                optio.max_iter = 0;
                
                [G_step2{ii},R_step2{ii},H_step2{ii},~,~,~,~,~,epsi{ii}] = SLS_core(data,netw,n_arx,optio);

            end

            

            % Compute performance of all estimates based on PEM criterion
            for ii = 1:1:length(ARX_orders)
                % cost of step 2: consistency for each ARX order
                J_PEM(1,ii) = PEM_cost(data,G_step2{ii},R_step2{ii},H_step2{ii});
            end
            
            % Select optimal ARX order based on PEM criterion
            [~,opt_ARX] = min(J_PEM);
            
            
            % Now do iterations
            n_arx = ARX_orders(opt_ARX);

            [G_step2,R_step2,H_step2,M2,G_step3,R_step3,H_step3,M3,~] = SLS_core(data,netw,n_arx,opt);

            
            

            % Select best iteration based on PEM criterion
            J_PEM = PEM_cost(data,G_step2,R_step2,H_step2);
            for iter = 1:1:size(G_step3,2)
                J_PEM(1+iter) = PEM_cost(data,G_step3{iter},R_step3{iter},H_step3{iter});
            end
            [~,opt_iter] = min(J_PEM);
            
            
            % Output
            if opt_iter == 1
                netw.G = G_step2;
                netw.R = R_step2;
                netw.H = H_step2;
                netw.D = M2.D;
                netw.NG = M2.NG;
                netw.NR = M2.NR;
                netw.NH = M2.NH;
            else
                netw.G = G_step3{opt_iter-1};
                netw.R = R_step3{opt_iter-1};
                netw.H = H_step3{opt_iter-1};     
                netw.D = M3{opt_iter-1}.D;
                netw.NG = M3{opt_iter-1}.NG;
                netw.NR = M3{opt_iter-1}.NR;
                netw.NH = M3{opt_iter-1}.NH;               
            end

            
        case 'very fast'
            % A good ARX order must be found by AIC criterion
            for ii = 1:1:length(ARX_orders)
                n_arx = ARX_orders(ii);
                
                [eps_hat,M_ARX] = SLS_step1(data,netw,n_arx);
                %AIC(ii) = check_AIC(eps_hat,M_ARX); % Uses AIC 
                AICc(ii) = check_AICc(eps_hat,M_ARX); % Uses AIC corrected for small sample sizes 
            end
            [~,opt_ARX] = min(AICc);
            
            
            % Now do step 2 and 3 of SLS without additional iterations
            [~,~,~,~,G_step3,R_step3,H_step3,M3,~] = SLS_core(data,netw,ARX_orders(opt_ARX),opt);
            netw.G = G_step3{1};
            netw.R = R_step3{1};
            netw.H = H_step3{1};
            netw.D = M3.D{1};
            netw.NG = M3.NG{1};
            netw.NR = M3.NR{1};
            netw.NH = M3.NH{1};
            
    end
    
    % Estimate the covariance matrix
    [~,Lambda] = PEM_cost(data,netw.G,netw.R,netw.H);
    netw.Lambda = Lambda;
    
    
    % In case of rank reduced noise, convert back to the non-square H format
    if exist('netw_original')
    if netw_original.p < netw_original.L 
        
        netw.p = netw_original.p; % revert to original rank reduced value
        netw.NH(:,netw_original.p+1:netw_original.L) = []; % remove excess zero columns
        netw.adjacencyH(:,netw_original.p+1:netw_original.L) = []; % remove excess columns
        netw.fixedH(:,netw_original.p+1:netw_original.L) = []; % remove excess columns
        netw.H(:,netw_original.p+1:netw_original.L) = []; % remove excess columns
        
        % obtain the rank reduced covariance matrix
        small_Lambda = netw.Lambda(1:netw_original.p,1:netw_original.p);
        netw.Gamma = netw.Lambda(netw_original.p+1:end,1:netw_original.p) / small_Lambda;
        netw.Lambda = small_Lambda;
        
    end
    end
     
     
% Produce requested output format
if nargout == 1
    varargout{1} = netw;
elseif nargout >= 3
    varargout{1} = netw.G;
    varargout{2} = netw.R;
    varargout{3} = netw.H;
    varargout{4} = netw.Lambda;
end


end