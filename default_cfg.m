function cfg = default_cfg(lambda, GauParam)
%   Default configuration for Gaussian spectral separation
%   Defines Gaussian centers (mu) and widths (sigma) for pure end-members
%
% Inputs:
%   lambda   : wavelength vector [nm]
%   GauParam : Sequence of Gaussian Fit parameters
% Output:
%   cfg : struct with fields for each endmember

    cfg.lambda = lambda(:);
    
    % Component Gaussians: e.g., Pure Wtar, CHl, TSM
    cfg.GauParam = GauParam;

    % Baseline terms included by default
    cfg.use_baseline = true;          % include b0 + b1*lambda

    % Optimization settings
    cfg.max_iter = 500;

    % must be 'none',  'off',  'iter',  'iter-detailed',  'final',  or 'final-detailed'
    cfg.display  = 'off';
end
