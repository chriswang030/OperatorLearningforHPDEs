function [U, S, V, mvs] = rsvd(A, k, p, q, Wl, Wr, K)
%% RSVD
%   Randomized SVD (Algorithm 1). 
%   (N. Halko, P.-G. Martinsson, J. A. Tropp, SIAM Review, 2011)
%
%   Input:
%   * A   : n x m matrix
%   * k   : target rank
%   * p   : oversampling
%   * q   : iterations
%   * (optional) Wl : n x 1 column vector of inner product weights for
%                     left-multiplication with A. Default is uniform 
%                     weight Wl = ones(n,1).
%   * (optional) Wr : m x 1 column vector of inner product weights for
%                     right-multiplication with A. Default is uniform 
%                     weight Wr = ones(m,1).
%   * (optional) K  : n x n covariance kernel
% 
%   Output:
%   * U   : n x (k+p) orthonormal columns, left singular vectors
%   * S   : (k+p) x (k+p) diagonal matrix of singular values
%   * V   : n x (k+p) orthonormal columnsm, right singular vectors
%   * mvs : number of matvecs

arguments
    A  (:,:) double
    k  (1,1) double
    p  (1,1) double
    q  (1,1) double
    Wl (:,1) double = ones(size(A,1),1)
    Wr (:,1) double = ones(size(A,2),1)
    K  (:,:) double = []
end

m = size(A,2);

if isempty(K)
    X = randn(m,k+p);
else
    X = mvnrnd(zeros(1,m), K, k+p)';
end

if q > 0
    Y = 1./Wl .* ((Wl.*A)*(Wr.*A'))^q * (Wl.*A) * (Wr.*X);
else
    Y = A * (Wr.*X);
end
[Q, ~] = qrw(Y, Wl);

B = Q' * (Wl.*A);
[Q1, R1] = qrw(B', Wr);
[U1, S, V1] = svd(R1');
U = Q * U1;
V = Q1 * V1;

mvs = (2*q+2)*(k+p);

end