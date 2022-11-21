function [G_est_fixed] = checkEstimation(G_est,fitPercentThreshold)

G_est_fixed = G_est;

[m,n] = size(G_est_fixed);

for i = 1:m
    for j = 1:n
        if (isempty(G_est_fixed(i,j).Report.Fit.FitPercent) || G_est_fixed(i,j).Report.Fit.FitPercent < fitPercentThreshold)
            G_est_fixed(i,j) = 0;
        end
    end
end

end