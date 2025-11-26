function Phi = build_basis(lambda, mu, sigma)
    n = numel(mu);
    Phi = zeros(numel(lambda), n);
    for k = 1:n
        Phi(:,k) = exp(-0.5*((lambda - mu(k))/sigma(k)).^2);
    end
end
