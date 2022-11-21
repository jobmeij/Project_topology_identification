%% SLS step 3
% Solving with use of cvx to translate easily to a solver

function [G_step3,R_step3,H_step3,M3] = SLS_step3_cvx(est_data,netw,NH_step2,optio)


	% Determine size of parameter vector based on network model structure
	theta_size = number_of_parameters(netw);
        
    % If the noise model is rank-reduced, the numerator is not invertible, and this is circumvented here
    NH_mod = NH_step2;
    for ii = 1:1:netw.L
        NH_mod{ii,ii}(1) = 1;
    end
    
    
    % Initialize NH inverse for step 3
    den = NH_mod;
    for colcol = 1:1:netw.L
        for rowrow = 1:1:netw.L
            den{colcol,rowrow} = 0*den{colcol,rowrow};
            den{colcol,rowrow}(1) = 1;
        end
    end
    NH_step_k = tf(NH_mod,den,est_data.Ts); 


    for iter = 1:1:optio.max_iter
        % Get inverse
        NH_inverse = inv(NH_step_k);
                
        NH_inverse_impulse = impulse(NH_inverse);
        NH_inverse_impulse = permute(NH_inverse_impulse,[2,3,1]);

        
        % if impulse response of noise model is long then code breaks...
        if size(NH_inverse_impulse,3) >= est_data.N
            break;
        end


        % Compute step 3 using cvx to translate to solver
        cvx_begin quiet
            variable theta(theta_size,1);
            J = cost_step3(theta,NH_inverse_impulse,netw,est_data);
            minimize J
        cvx_end % Solved!        
        theta_step3(:,iter) = theta; % save the variable

        

        % Convert to network model
        [D3d,NG3d,NR3d,NH3d] = parameters_to_D3model(theta,netw);
        [M3{iter}.D,M3{iter}.NG,M3{iter}.NR,M3{iter}.NH] = D3model_to_polynomial(D3d,NG3d,NR3d,NH3d,netw);
        [G_step3{iter},R_step3{iter},H_step3{iter}] = model_to_transfer(M3{iter},est_data.Ts,netw);

        
        % Check stability of inverse, stop iterations if unstable
        if sum(abs(eig(NH_inverse)) > 1)
            break;
        end
        

        % Check whether estimated noise model has converged
        den = NH_mod;
        for colcol = 1:1:netw.L
            for rowrow = 1:1:netw.L
                den{colcol,rowrow} = 0*den{colcol,rowrow};
                den{colcol,rowrow}(1) = 1;
            end
        end
        NH_new = tf(H_step3{iter}.num,den,est_data.Ts);
        
        
        if norm(NH_step_k - NH_new) > optio.T_conv
            % New NH and initial parameter
            NH_step_k = NH_new;
            theta_init_k = theta_step3(:,iter);
        else
            % Converged
            break;
        end
    
    end

end