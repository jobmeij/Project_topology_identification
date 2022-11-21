%% Function to retrieve the parameter vector from a network model

function count = number_of_parameters(netw)
    
    count = 0;

    
    for colcol = 1:1:netw.L
        % Count parameters in D
        count = count + size(netw.D{colcol},2)-1;
        
        % Count parameters in G
        for rowrow = 1:1:netw.L
            if size(netw.NG{colcol,rowrow},2) > 0
                count = count + size(netw.NG{colcol,rowrow},2)-1;
            end
        end
        
        % Count parameters in R
        for rowrow = 1:1:netw.K
            count = count + size(netw.NR{colcol,rowrow},2);
        end
        
        % Count parameters in H
        for rowrow = 1:1:netw.L
            if size(netw.NH{colcol,rowrow},2) > 0
                count = count + size(netw.NH{colcol,rowrow},2)-1;
            end
        end
    end


    


end