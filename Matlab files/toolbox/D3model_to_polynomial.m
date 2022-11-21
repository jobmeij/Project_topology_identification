%% Function to compute polynomial representation from parameter vector theta

% The polynomials are D, N_G, N_R, N_H

function [D,NG,NR,NH] = D3model_to_polynomial(D3d,NG3d,NR3d,NH3d,netw)
    
% Make D monic
D3d(:,:,1) = eye(size(D3d,1));

% In H make the entries with parameters monic

NH3d(:,:,1) = diag(diag(double(sum(NH3d.^2,3) ~= 0)));
    

% D polynomial
for ii = 1:1:netw.L
    D{ii,1} = squeeze(D3d(ii,ii,:))';
end


% N_G polynomial
for colcol = 1:1:netw.L
    for rowrow = 1:1:netw.L
        NG{colcol,rowrow} = squeeze(NG3d(colcol,rowrow,:))';
    end
end

% N_R polynomial
for colcol = 1:1:netw.L
    for rowrow = 1:1:netw.K
        NR{colcol,rowrow} = squeeze(NR3d(colcol,rowrow,:))';
    end
end

% N_H polynomial
for colcol = 1:1:netw.L
    for rowrow = 1:1:netw.L
        NH{colcol,rowrow} = squeeze(NH3d(colcol,rowrow,:))';
    end
end

end