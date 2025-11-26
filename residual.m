function r = residual(p_soft, Rmix, Phi_all, Phi_bl)
    % Residual for lsqnonlin with softplus-constrained amplitudes
     p = as_softplus(p_soft);
    nG = size(Phi_all,2);

    A = p(1:nG);
    if ~isempty(Phi_bl)
        b0 = p(nG+1);
        b1 = p(nG+2);
    else
        b0 = 0; b1 = 0;
    end

    Rfit = Phi_all*A + b0 + b1*((1:length(Rmix))');
    r = Rfit - Rmix;
    
end
