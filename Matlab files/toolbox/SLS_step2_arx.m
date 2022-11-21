%% SLS step 2

function [G_step2,R_step2,H_step2,M2] = SLS_step2_arx(est_data,netw)

    % shorthand notations
    L = netw.L;
    K = netw.K;
    
    % Network model translated to ARX model
    for ii = 1:1:L
        A{ii,ii} = netw.D{ii,1};
    end
    
     
    B = [netw.NG netw.NR netw.NH];
    arx_poly = idpoly(A,B);
    
    % Monic noise model enforced
    for rowrow = L+K + (1:1:L)
        for colcol = 1:1:L
            arx_poly.Structure.B(colcol,rowrow).Free(1) = 0;
        end
    end
    
    % Modify data object appropriately with known dynamics
    w_minus_known = est_data.y -lsim(netw.fixedG,est_data.y) -lsim([netw.fixedR netw.fixedH],est_data.u); 
    arx_data = iddata(w_minus_known,[est_data.y est_data.u],est_data.Ts);
    
    % estimate with monic N_H
    Marx = arx(arx_data,arx_poly);

    

    % output is as transfer functions of the correct format
    [M2.D,M2.NG,M2.NR,M2.NH] = ARX_to_model(Marx,netw);
    [G_step2,R_step2,H_step2] = model_to_transfer(M2,est_data.Ts,netw);

end