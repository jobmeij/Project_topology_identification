function net = network_polynomial(L,K,p,order,topo)

    net.L = L;
    net.K = K;
    net.p = p;
    net.orders = order;

    % Assign D polynomial with normalized denominator
    tempD = [ones(L,1) zeros(net.L,order.D)];
    net.D = mat2cell(tempD,ones(net.L,1)); 

    % Assign N_G polynomial with 0
    for colcol = 1:1:net.L
        for rowrow = 1:1:net.L
            if topo.NG(colcol,rowrow) 
                net.NG{colcol,rowrow} = zeros(1,order.NG+1); % default parameters if there is a link
                net.NG{colcol,rowrow}(2) = 1;
            else
                net.NG{colcol,rowrow} = zeros(1,0); % remain empty when there is no link
            end
        end
    end
 
    % Assign N_R polynomial with 0
    for colcol = 1:1:net.L
        for rowrow = 1:1:net.K
            if topo.NR(colcol,rowrow) 
                net.NR{colcol,rowrow} = [1 zeros(1,order.NR)];% default parameters if there is a link
            else
                net.NR{colcol,rowrow} = zeros(1,0);% remain empty when there is no link
            end
        end
    end
 
    % Assign N_H polynomial with 0
    for colcol = 1:1:net.L
        for rowrow = 1:1:net.p
            if topo.NH(colcol,rowrow) 
                net.NH{colcol,rowrow} = zeros(1,order.NH+1);% default parameters if there is a link
            else
                net.NH{colcol,rowrow} = zeros(1,0);% remain empty when there is no link
            end
        end
        if colcol <= net.p
            net.NH{colcol,colcol}(1) = 1; % Make it monic!
        end
    end
end