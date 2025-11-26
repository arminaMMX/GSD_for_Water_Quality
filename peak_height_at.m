function [peakH, peakLam] = peak_height_at(lambda, S, center_nm, window_nm)
    % Extract local peak height around specified center using a narrow window
    idx = lambda >= (center_nm - window_nm) & lambda <= (center_nm + window_nm);
    [peakH, ii] = max(S(idx));
    lamLoc = lambda(idx);
    peakLam = lamLoc(ii);
end
