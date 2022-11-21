%% SLS step 1

function [eps_hat,M_ARX] = SLS_step1(data,netw,n_arx,optio)

    % Specify ARX model orders 
    order.na = ones(netw.L) * n_arx;
    order.nb = ones(netw.L,netw.K) * n_arx;
    order.nk = zeros(netw.L,netw.K);
    
    ARX_orders = [order.na order.nb order.nk];

    
    % Possible regularization
    if strcmp(optio.regularization,'none')
        arxopt = arxOptions('OutputWeight',optio.WLS);  % Add Weighted Least Squares possibility
    else        
        % regularization kernel
        reg_option = arxRegulOptions('RegularizationKernel',optio.regularization);
        [lambda,Kern] = arxRegul(data,ARX_orders,reg_option);
        
        % place kernel in options
        arxopt = arxOptions('OutputWeight',optio.WLS);  % Add Weighted Least Squares possibility
        arxopt.Regularization.Lambda = lambda;
        arxopt.Regularization.R = Kern;
    end   
    
    % Eestimate the model
    M_ARX = arx(data,ARX_orders,arxopt);
    
    % Recover innovation sequence
    [eps_hat,~] = resid(data,M_ARX);
    
    
end