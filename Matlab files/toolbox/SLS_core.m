%% Sequential Least Squares core function

function [G_step2,R_step2,H_step2,M2,G_step3,R_step3,H_step3,M3,eps_hat] = SLS_core(data,netw,arx_ord,optio)



    %% Step 1: ARX
    [eps_hat,~] = SLS_step1(data,netw,arx_ord,optio);
    
    
    %% Step 2: Consistency
    % estimated innovation is an input in the new data set
    est_data = iddata(data.y,[data.u eps_hat.OutputData],data.Ts);
    

    [G_step2,R_step2,H_step2,M2] = SLS_step2_arx(est_data,netw);

    
    %% Step 3: Efficiency and iterations
    
    % Check if there is an efficiency step
    if optio.max_iter >= 1
        % Estimate step 3 and iterations
        [G_step3,R_step3,H_step3,M3] = SLS_step3_cvx(est_data,netw,M2.NH,optio);
    else
        G_step3 = [];
        R_step3 = [];
        H_step3 = [];
        M3 = [];
    end

end