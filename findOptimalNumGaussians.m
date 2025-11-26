function [bestNumPeaks, results] = findOptimalNumGaussians(x, y, maxPeaks, MinPeakProminence, MinPeakHeight,doPlot)
% Find optimal number of Gaussians using information criteria

    if nargin < 6
        doPlot = true;
    end

    results = struct();
    aic_values = zeros(maxPeaks, 1);
    bic_values = zeros(maxPeaks, 1);
    r2_values = zeros(maxPeaks, 1);
    
    for n = 1:maxPeaks
        try
            [fitresult, gof, ~] = fitGaussianMix(x, y, n, MinPeakProminence, MinPeakHeight, false);
            n_params = 3 * n; % parameters per Gaussian
            
            % Calculate AIC and BIC
            n_points = length(y);
            rss = sum((y - feval(fitresult, x)).^2);
            aic_values(n) = n_points * log(rss/n_points) + 2 * n_params;
            bic_values(n) = n_points * log(rss/n_points) + n_params * log(n_points);
            r2_values(n) = gof.rsquare;
            
            results(n).numPeaks = n;
            results(n).aic = aic_values(n);
            results(n).bic = bic_values(n);
            results(n).r2 = r2_values(n);
        catch
            aic_values(n) = Inf;
            bic_values(n) = Inf;
            r2_values(n) = 0;
        end
    end
    
    % Find optimal number (minimize AIC/BIC)
    [~, bestAIC] = min(aic_values);
    [~, bestBIC] = min(bic_values);
    
    bestNumPeaks = min(bestAIC, bestBIC);
    
    % Plot comparison
    if doPlot
        figure;
        subplot(2,1,1);
        plot(1:maxPeaks, aic_values, 'bo-', 'LineWidth', 2);
        hold on;
        plot(1:maxPeaks, bic_values, 'ro-', 'LineWidth', 2);
        xlabel('Number of Gaussians');
        ylabel('Information Criterion');
        legend('AIC', 'BIC');
        grid on;

        subplot(2,1,2);
        plot(1:maxPeaks, r2_values, 'go-', 'LineWidth', 2);
        xlabel('Number of Gaussians');
        ylabel('R²');
        grid on;
    end
end