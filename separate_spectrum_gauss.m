function out = separate_spectrum_gauss(Rmix, cfg, EMs)
% Separate a mixed spectrum into contributions from end-members
% using Gaussian bases with fixed centers/widths.
%
% Inputs:
%   lambda : [N x 1] wavelength vector (nm)
%   Rmix   : [N x 1] measured reflectance spectrum
%   cfg    : configuration struct from default_cfg(lambda)
%   Ems    : end-members spectrs
%
% Outputs:
%   out : struct with reconstructed spectra, parameters, RMSE, etc.

%----- Combine all Gaussians ----------------------------------------------
p = size(cfg.GauParam,2);
Phi_all = [];

lambda = cfg.lambda;

for i=1:p
    Phi_ = build_basis(lambda, cfg.GauParam{i}(:,1), cfg.GauParam{i}(:,2));
    Phi_all = [Phi_all, Phi_];
end

%--------------------------------------------------------------------------

if cfg.use_baseline
    Phi_bl = [ones(size(lambda)) lambda];
else
    Phi_bl = [];
end

%----- Full design matrix -------------------------------------------------
Phi = [Phi_all, Phi_bl];

%----- Initial guesses for amplitudes and baseline ------------------------
nG = size(Phi_all,2);
p0 = [0.1*ones(nG,1)];

if ~isempty(Phi_bl)
    p0 = [p0; 0; 0];
end

%----- Use lsqnonlin with transformation to enforce nonnegativity on Gaussian amps
opts = optimoptions('lsqnonlin','MaxIterations',cfg.max_iter,'Display',cfg.display);
p_soft0 = as_softplus_inv(p0); % unconstrained initialization

fun = @(p_soft) residual(p_soft, Rmix, Phi_all, Phi_bl);
p_soft_hat = lsqnonlin(fun, p_soft0, [], [], opts);

%----- Retrieve physical parameters ---------------------------------------
p_hat = as_softplus(p_soft_hat);

%----- Extract amplitudes -------------------------------------------------
A = p_hat(1:nG);
if ~isempty(Phi_bl)
    b0 = p_hat(nG+1);
    b1 = p_hat(nG+2);
else
    b0 = 0;
    b1 = 0;
end

%----- Reconstruct --------------------------------------------------------
S_gauss = Phi_all .* A';
Rfit    = sum(S_gauss,2) + b0 + b1*lambda;

%----- Reconstruct separated spectra --------------------------------------
hat_EMs = nan(size(lambda,1),p);
nA = 0;
for i=1:p
    nB = size(cfg.GauParam{i},1);
    A_Hat = p_hat(nA+1:nA+nB);
    nA = nA + nB;

    Phi_ = build_basis(lambda, cfg.GauParam{i}(:,1), cfg.GauParam{i}(:,2));
    hat_EMs(:,i) = Phi_ * A_Hat;
end

%----- Abundances ---------------------------------------------------------
abundances = zeros(1,p);

opts = optimoptions('lsqlin','Display','off');
Aeq = ones(1,p);
beq = 1;
lb = zeros(p,1);

b = Rmix;
a = lsqlin(EMs, b, [], [], Aeq, beq, lb, [], [], opts);
abundances(1,:) = a';
recon = EMs * a;
rmse_sim = sqrt(mean((recon - b).^2));
sam_sim = acosd( (recon' * b) / (norm(recon)*norm(b)) );

%----- Pack results -------------------------------------------------------
out.Rfit       = Rfit;
out.rmse       = rmse_sim;
out.sam        = sam_sim;
out.params.b0  = b0;
out.params.b1  = b1;
out.cfg        = cfg;
out.abundances = abundances;
out.hat_EMs    = hat_EMs;

if int16(sum(abundances)) > 1
    error('An error rized, Sigma Abundances > 1');
end


end
