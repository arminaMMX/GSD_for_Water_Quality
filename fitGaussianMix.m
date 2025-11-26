function [fitresult, gof, FitParam, coeffTable] = fitGaussianMix(x, y, numGauss, MinPeakProminence, MinPeakHeight, showPlot, Title)
% Fit spectrum y(x) with a sum of numGauss Gaussian peaks
%
% Inputs:
%   x        - wavelength or x-data (Nx1)
%   y        - spectrum or signal (Nx1)
%   numGauss - number of Gaussians (1..8 supported by fittype)
%   MinPeakProminence - findpeaks option
%   MinPeakHeight - findpeaks option
%   showPlot - boolean to show plot (default: true)
%
% Outputs:
%   fitresult - cfit object from MATLAB fit()
%   gof       - goodness-of-fit struct (R2, RMSE, etc.)
%   FitParam  - structure with Gaussian parameters
%   coeffTable - table of Gaussian parameters

    % Default to show plot if not specified
    if nargin < 6
        showPlot = true;
    end

    
    % Ensure column vectors
    x = x(:);
    y = y(:);

    % Remove NaN values
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);

    % Prepare data
    [xData, yData] = prepareCurveData(x, y);

    % Find initial peaks for better starting guesses
    [pks, locs, widths] = findpeaks(yData, xData, ...
        'MinPeakProminence', MinPeakProminence * max(yData), ...
        'MinPeakHeight', MinPeakHeight * max(yData), ...
        'SortStr', 'descend');
    
    % If we found fewer peaks than requested, adjust
    numFound = min(numGauss, length(pks));
    
    % Create fit options with better starting guesses
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'off';
    
    % Set starting points based on detected peaks
    if numFound > 0
        startPoint = [];
        lowerBounds = [];
        upperBounds = [];
        
        for i = 1:numGauss
            if i <= numFound
                % Use detected peak parameters
                A_start = pks(i);
                mu_start = locs(i);
                sigma_start = widths(i) / (2*sqrt(2*log(2))); % Convert FWHM to sigma
                c_start = sigma_start * sqrt(2); % Convert to c parameter
            else
                % Use reasonable defaults for additional peaks
                A_start = 0.1 * max(yData);
                mu_start = min(xData) + (max(xData)-min(xData)) * i/(numGauss+1);
                c_start = (max(xData)-min(xData)) / (2*numGauss);
            end
            
            startPoint = [startPoint, A_start, mu_start, c_start];
            lowerBounds = [lowerBounds, 0, min(xData), 0];
            upperBounds = [upperBounds, Inf, max(xData), Inf];
        end
        
        opts.StartPoint = startPoint;
        opts.Lower = lowerBounds;
        opts.Upper = upperBounds;
    end

    % Choose Gaussian model
    modelName = sprintf('gauss%d', numGauss);
    ft = fittype(modelName);

    % Fit the model
    try
        [fitresult, gof] = fit(xData, yData, ft, opts);
    catch ME
        warning('Fit failed with initial guesses. Trying with default parameters.');
        opts = fitoptions('Method', 'NonlinearLeastSquares', 'Display', 'off');
        [fitresult, gof] = fit(xData, yData, ft, opts);
    end

    % Extract coefficients
    coeffs  = coeffvalues(fitresult);
    confInt = confint(fitresult);
    
    nPeaks = numGauss;
    A      = coeffs(1:3:end); % amplitudes
    mu     = coeffs(2:3:end); % centers
    c      = coeffs(3:3:end); % width parameter
    sigma  = c ./ sqrt(2); % true σ
    FWHM   = 2.355 .* sigma;
    
    % Calculate confidence intervals
    A_ci     = confInt(:,1:3:end);
    mu_ci    = confInt(:,2:3:end);
    sigma_ci = confInt(:,3:3:end) ./ sqrt(2);

    % Calculate integrated area under each Gaussian
    areas = A .* sigma .* sqrt(2*pi);

    % Store parameters
    FitParam.Amplitude    = A;
    FitParam.Center       = mu;
    FitParam.Width        = c;
    FitParam.Sigma        = sigma;
    FitParam.FWHM         = FWHM;
    FitParam.Area         = areas;
    FitParam.Amplitude_CI = A_ci;
    FitParam.Center_CI    = mu_ci;
    FitParam.Sigma_CI     = sigma_ci;

    % Create output table
    coeffTable = table((1:nPeaks)', A', mu', sigma', FWHM', areas', ...
        'VariableNames', {'Peak#', 'Amplitude', 'Center', 'Sigma', 'FWHM', 'Area'});

    if showPlot
        disp(coeffTable);
        
        % Display goodness of fit
        fprintf('Goodness of fit:\n');
        fprintf('  R² = %.2f\n', gof.rsquare);
        fprintf('  Adjusted R² = %.2f\n', gof.adjrsquare);
        fprintf('  RMSE = %.4f\n', gof.rmse);
        
        % Plot results
        figure;
        hold on;
        
        % Plot individual Gaussians
        colors = lines(nPeaks);
        y_fit_total = feval(fitresult, xData);
        
        for i = 1:nPeaks
            % Calculate individual Gaussian component
            individual_gauss = A(i) * exp(-((xData - mu(i)) / c(i)).^2);
            
            % Plot individual Gaussian
            plot(xData, individual_gauss, '-', 'Color', colors(i,:), ...
                'LineWidth', 1.5, 'HandleVisibility','off');
            
            % Mark peak center
            plot(mu(i), A(i), 'o', 'Color', colors(i,:), ...
                'MarkerFaceColor', colors(i,:), 'MarkerSize', 8, ...
                'HandleVisibility', 'off');
            
            % Show FWHM
            half_max = A(i)/2;
            fwhm_x = [mu(i)-FWHM(i)/2, mu(i)+FWHM(i)/2];
            plot(fwhm_x, [half_max, half_max], '-', 'Color', colors(i,:), ...
                'LineWidth', 2, 'HandleVisibility', 'off');
        end
        
        % Plot original data and total fit
        plot(xData, yData, 'k-', 'LineWidth', 3.0, 'DisplayName', 'Original Data');
        plot(xData, y_fit_total, 'r--', 'LineWidth', 3.0, 'DisplayName', 'Total Fit');

        hold off;
        xlim([min(x)*0.95 max(x)*1.05])

        % Add text box with fit statistics
        text(0.7, 0.98, sprintf('R^2 = %.2f\nRMSE = %.4f', gof.rsquare, gof.rmse), ...
            'Units', 'normalized', 'VerticalAlignment', 'top','FontSize',18);

        xlabel('Wavelength (nm)');
        ylabel('{\itRrs} (sr^{-1})');
        legend('Location', 'northwest');
        title(Title)
        grid on;
        box on;
        
    end
end