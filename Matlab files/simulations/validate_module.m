%% Validation of transfer function

function [fit,err] = validate_module(netw_estimated,G0)

    L = size(G0,1);
    K = size(G0,2);

    for mc = 1:1:size(netw_estimated,2)
        for colcol = 1:1:L
            for rowrow = 1:1:K
                err(colcol,rowrow,mc) = norm(netw_estimated{mc}.G(colcol,rowrow)-G0(colcol,rowrow));
                if norm(G0(colcol,rowrow)) > 0
                    fit(colcol,rowrow,mc) = 1 - (err(colcol,rowrow,mc) / norm(G0(colcol,rowrow)));
                else
                    fit(colcol,rowrow,mc) = 0;
                end
            end
        end
    end


end