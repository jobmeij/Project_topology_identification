%% Function that computes the cost function in step 2 of SLS
function J = cost_step3(theta,Hi,netw,data)

    eps_L2 = compute_eps_L2(theta,netw,data);
    eps_L3 = compute_eps_L3(Hi,eps_L2);

    

    for ii = 1:1:netw.L
        Jvec(ii) = eps_L3(ii,:) * eps_L3(ii,:)';
    end
    J = sum(Jvec);

end