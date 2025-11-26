function x = as_softplus_inv(y)
    % Inverse softplus (approximate)
    x = log(exp(y) - 1 + eps);
end
