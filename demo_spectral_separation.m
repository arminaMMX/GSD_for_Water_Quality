%% Gaussian Spectral Separation
% Author: Masoud Moradi, 2025
% Purpose:
%   MATLAB code to separate mixed reflectance spectra
%   into Chlorophyll-a (Chla) and Total Suspended Matter (TSM) components using
%   Gaussian spectral decomposition, then estimate concentrations.
%
% Notes:
%   • Important bands are near ~650 nm (Chla dominated) and ~810 nm (TSM dominated).
%     These appear in the default centers below; tune to your sensor and water type.
%   • You can expand the number of Gaussians per component to better fit
%     your instrument and region.
%   • After separation, you can map peak heights to concentrations via
%     linear or log-linear regressions using your lab measurements.
%
% Inputs:
%   • wl  : wavelength vector `lambda` [nx1] in nm
%   • Rrs : [nx1]
%   • OPTIONAL: 
%       • Rrs_pure_Chl : pure component spectra for Chla-only
%       • Rrs_pure_TSM : pure component spectra for Chla-only
%       * Otherwise, defaults are used
%
% Model (analytical):
%   R(λ) = b0 + b1*λ + w_c * S_c(λ) + w_s * S_s(λ)
%   where the component spectra are each a sum of Gaussians
%   S_c(λ) = Σ_k A_{c,k} * exp(-(λ-μ_{c,k})^2 / (2σ_{c,k}^2))
%   S_s(λ) = Σ_k A_{s,k} * exp(-(λ-μ_{s,k})^2 / (2σ_{s,k}^2))
%   Parameters μ (centers) and σ (widths) are taken from pure spectra fits
%   or set to defaults; amplitudes (A), weights (w_c, w_s) and baseline (b0,b1)
%   are solved per mixed spectrum with non-negativity constraints on amplitudes.
%
% -------------------------------------------------------------------------
% Sample Rrs taken from the central part of the Dutch Wadden Sea.
% Chlta = 38.0 mg m-3, TSM = 53.3 g m-3, aCDOM(440) = 1.11 m-1
% Reference:
% Hommersom, A., Peters, S., Wernand, M. R., & de Boer, J. (2009). 
% Spatial and temporal variability in bio-optical properties of the Wadden Sea. 
% Estuarine, Coastal and Shelf Science, 83(3), 360-370
% -------------------------------------------------------------------------
clc;
clear;
% --- Load data (for demo only) -------------------------------------------
myFile = './Sample_Spectra.csv';
T      = readtable(myFile);

lambda = T.Wavelengths;                   % wavelengths [nm]
Rmix   = T.Rrs_sr_1;                      % Sample Rrs spectrum (sr-1)
MEMs   = double(table2array(T(:,3:11)));  % End-Members

% --- Gaussian decomposition of MEMs --------------------------------------
% 1. Evaluate the number of Gaussians
maxPeaks = 8;
MinPeakHeight = 0.01;
MinPeakProminence = 0.01;

numGauss = zeros(size(MEMs,2),1);

for i = 1:size(MEMs,2)
    [bestNumPeaks, results] = findOptimalNumGaussians(lambda, MEMs(:,i), ...
                            maxPeaks, MinPeakProminence, MinPeakHeight, false);
    numGauss(i) = bestNumPeaks;
end

% 2. Gaussian parameters of MEMs
showPlot = false;

GauParams = cell(1, size(MEMs,2));

for i=1:size(MEMs,2)
    [~, ~, ~, GauTable] = fitGaussianMix(lambda, MEMs(:,i), ...
                 numGauss(i), MinPeakProminence, MinPeakHeight, showPlot,sprintf('MEM #%d',i));

    GauParam{i} = double(table2array(GauTable(:,3:4)));
end

% --- Gaussian decomposition of Rmix --------------------------------------
% 1. Define default centers and widths
cfg = default_cfg(lambda, GauParam);

% 2. Spectral Separation with defaults
out = separate_spectrum_gauss(Rmix, cfg, MEMs);

fprintf('Fit RMSE: %.4f\nFit SAM: %3.f\n', out.rmse, out.sam);


% --- Visualize -----------------------------------------------------------
figure('Position',[38,205,835,640]);
hold on;
box on;
plot(lambda, Rmix          , 'b-' , 'linewidth', 2.5  , 'DisplayName','Mixed');
plot(lambda, out.Rfit      , 'r--', 'linewidth', 2.5  , 'DisplayName','Fit');
for i=1:size(MEMs,2)
    plot(lambda, out.hat_EMs(:,i), 'linewidth', 1.5, 'DisplayName', ...
                sprintf('Separated EM #%d (abundance: %.2f)',i, out.abundances(i)));
end
hold off
xlabel('\lambda (nm)');
ylabel('Reflectance (sr_{-1})'); 
legend show;
grid on
title('Gaussian Spectral Separation (demo)');

% Peak features (for downstream concentration models)
for i=1:size(MEMs,2)
    [h650, ~] = peak_height_at(lambda, out.hat_EMs(:,i), 650, 15);
    [h811, ~] = peak_height_at(lambda, out.hat_EMs(:,i),   811, 20);
    fprintf('Peak(Chla@650 nm) ≈ %.4f, Peak(SS@810 nm) ≈ %.4f\n', h650, h811);
end

