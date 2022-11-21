%% 
% Function to compute network polynomial matrices D,N_G,N_R,N_H from the 
% ARX estimate A,B via these relations
% A = D - N_G
% B = [R  H]




function [D,NG,NR,NH] = ARX_to_model(M,netw)
    
    
    % D polynomial
    for ii = 1:1:netw.L
        D{ii,1} = M.A{ii,ii};
    end
    
    NG = M.B(1:netw.L,1:netw.L);
    NR = M.B(1:netw.L,netw.L + (1:netw.K));
    NH = M.B(1:netw.L,netw.L+netw.K + (1:netw.L));

    
end