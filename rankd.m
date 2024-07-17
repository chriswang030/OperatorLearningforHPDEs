function [high, mvs] = rankd(A, k, p, q, tol, C, Wl, Wr, K)
%% RANKD
%   Decide whether a matrix A is sufficiently low-rank, numerically, using 
%   rSVD (Algorithm 3).
% 
%   Input:
%   * A    : n x m matrix
%   * k    : target rank
%   * p    : oversampling
%   * q    : power iterations
%   * tol  : error tolerance
%   * (optional) C  : positive scalar value for rank test
%   * (optional) Wl : n x 1 column vector of inner product weights for
%                     left-multiplication with A. Default is uniform 
%                     weight Wl = ones(n,1).
%   * (optional) Wr : m x 1 column vector of inner product weights for
%                     right-multiplication with A. Default is uniform 
%                     weight Wr = ones(m,1).
%   * (optional) K  : n x n covariance kernel
% 
%   Output:
%   * high : logical, true if high rank, else false
%   * mvs  : number of matvecs

arguments
    A   (:,:) double
    k   (1,1) double
    p   (1,1) double
    q   (1,1) double
    tol (1,1) double
    C   (1,1) double = 4
    Wl  (:,1) double = ones(size(A,1),1)
    Wr  (:,1) double = ones(size(A,2),1)
    K   (:,:) double = []
end

Aq = 1./Wl .* ((Wl.*A)*(Wr.*A'))^q * (Wl.*A);
[U, ~, ~, mvs] = rsvd(Aq, k, p, 0, Wl, Wr, K);
mvs = mvs * (2*q+1);

r = min([k size(A) size(U,2)]); % in case lower rank than we thought

U = U(:,1:r);
B = U' * A;
mvs = mvs + r;

S = svd(B);
high = (S(k) > C * tol * S(1));

end