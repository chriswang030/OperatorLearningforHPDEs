function [Q, R] = qrw(X, W, tol)
%% QRW
%   QR using modified Gram-Schmidt with weighted inner product. The 
%   returned matrices satisfy Q*R = X; Q has r orthonormal columns with 
%   respect to the weighted inner product, and zeroes in the remaining 
%   columns, where r = RANK(X,TOL).
%   (G. W. Stewart, Matrix Algorithms, Vol. 1, SIAM, 1998).
%
%   Weighted inner product sped up with discussion at 
%   https://stackoverflow.com/a/27126337.
% 
%   Input:
%   * X : n x m real full rank matrix, n >= m
%   * (optional) W   : n x 1 column vector of inner product weights
%                      = ones(n,1)
%   * (optional) tol : tolerance for numerical rank = 1e-12
% 
%   Output:
%   * Q : n x m matrix with r orthonormal columns with respect to W, where
%         r = RANK(X,TOL), and zeroes in the remaining columns
%   * R : n x n upper triangular matrix such that X = Q*R

arguments
    X   (:,:) double
    W   (:,1) double = ones(size(X,1),1)
    tol (1,1) double = 1e-12
end

[n,m] = size(X);
Q = zeros(n,m);
R = zeros(m,m);
S = zeros(n,m); % store some weighted values

for j = 1:m
    Q(:,j) = X(:,j);
    for i = 1:j-1
        R(i,j) = Q(:,j)'*S(:,i);
        Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
    end
    R(j,j) = sqrt(Q(:,j)'*(W.*Q(:,j)));
    if R(j,j) < tol
        Q(:,j) = zeros(n,1);
    else
        Q(:,j) = Q(:,j) / R(j,j);
    end
    S(:,j) = W .* Q(:,j);
end