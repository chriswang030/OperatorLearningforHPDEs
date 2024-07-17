function pass = test_qrw(tol)
% Test QR factorization with weighted inner product.

if nargin == 0
    tol = 1e-14;
end

rng(1);
X = randn(100,10);
[~, W] = legpts(100);
W = W';

% weighted
[Q, R] = qrw(X, W);
pass(1) = (norm(Q'*(W.*Q) - eye(10)) < tol); % orthonormality
pass(2) = (max(abs(X-Q*R), [], 'all') < tol); % correct product

% unweighted
[Q, R] = qrw(X);
pass(3) = (norm(Q'*Q - eye(10)) < tol); % orthonormality
pass(4) = (max(abs(X-Q*R), [], 'all') < tol); % correct product

end