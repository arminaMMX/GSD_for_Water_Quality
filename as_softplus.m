function y = as_softplus(x)
    % Softplus to ensure positivity
    y = log1p(exp(x));
end
