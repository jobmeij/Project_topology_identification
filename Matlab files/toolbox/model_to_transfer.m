%% Function to compute transfer functions from polynomials D, N_G, N_R, N_H



% The transfer functions are G,R,H

function [G,R,H] = model_to_transfer(polyno,Ts,netw)

    L = netw.L;
    K = netw.K;
    p = netw.p;
    
    % Create appropriate denominator cell
    D_wideL = repmat(polyno.D,1,L);
    D_wideK = repmat(polyno.D,1,K);
    D_widep = repmat(polyno.D,1,p);
    
    % assign transfer functions
    G = tf(polyno.NG,D_wideL,Ts) + netw.fixedG;
    R = tf(polyno.NR,D_wideK,Ts) + netw.fixedR;
    H = tf(polyno.NH,D_widep,Ts) + netw.fixedH;
    
    
end