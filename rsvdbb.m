function [U,S,V,mvs] = rsvdbb(solve,asolve,k,p,Nx,Nt,xtid,ysid,smooth)
%% RSVDBB
%   Randomized SVD (Algorithm 1), using a black-box numerical solver to 
%   generate input-output data.
%   (N. Halko, P.-G. Martinsson, J. A. Tropp, SIAM Review, 2011)
%
%   Input:
%   * solve  : @(f) -> double, forward numerical solver
%   * asolve : @(f) -> double, backward (adjoint) numerical solver
%   * k      : target rank
%   * p      : oversampling
%   * Nx     : number of spatial gridpoints (must be of form 2^n+1)
%   * Nt     : number of temporal gridpoints (must be of form 2^n+1)
%   * xtid   : [x(start) x(end) t(start) t(end)] grid indices for which
%              target subdomain to test
%   * ysid   : [y(start) y(end) s(start) s(end)] grid indices for which
%              source subdomain to test
%   * (optional) smooth : whether to artificially smooth random test
%                         functions = false
% 
%   Output:
%   * U   : n x k orthonormal columns, left singular vectors
%   * S   : k x k diagonal matrix of singular values
%   * V   : n x k orthonormal columns, right singular vectors
%   * mvs : number of matvecs

arguments
    solve  function_handle
    asolve function_handle
    k      (1,1) int16 
    p      (1,1) int16
    Nx     (1,1) double
    Nt     (1,1) double
    xtid   (1,4) double
    ysid   (1,4) double
    smooth logical = 0
end

% set basis for GP
v = @(x,k) sqrt(2) * sin(pi*k.*x); % ON basis of L2([0,1])
x = linspace(0,1,Nx);
t = linspace(0,1,Nt);
r1 = 10;
r2 = 10;
[xN,NNx] = meshgrid(x,1:r1);
[tN,NNt] = meshgrid(t,1:r2);

% setup
n1 = xtid(2)-xtid(1)+1;
n2 = xtid(4)-xtid(3)+1;
n3 = k+p;
assert(ysid(2)-ysid(1)+1 == n1);
assert(ysid(4)-ysid(3)+1 == n2);

% 2D trapezoidal rule weights
w = ones(n1,n2);
w(2:end-1,:) = w(2:end-1,:) + 1;
w(:,2:end-1) = w(:,2:end-1) + 1;
w(2:end-1,2:end-1) = w(2:end-1,2:end-1) + 1;
W = reshape(w,n1*n2,1)/(4*Nx*Nt);

% generate random inputs
M = randn(r1,r2,n3);
f = pagemtimes(pagemtimes(v(xN,NNx)',M),v(tN,NNt));

% restrict to ysdom
res = zeros(Nx,Nt,n3);
res(ysid(1):ysid(2),ysid(3):ysid(4),:) = 1;
f = f .* res;
if smooth
    f = smoothdata2(f,'gaussian',{(n1-1)/2, (n2-1)/2});
end

% solve and extract from xtdom: Y = A*X
u = solve(f);
Y = u(xtid(1):xtid(2),xtid(3):xtid(4),:);

% weighted QR: Y = Q*R
Y = reshape(Y,n1*n2,n3);
[Q,~] = qrw2(Y,W); % Q is (n1*n2)-by-n3

% extend by zero from xtdom, solve adjoint, extract from ysdom: B' = A'*Q
fq = zeros(size(f));
fq(xtid(1):xtid(2),xtid(3):xtid(4),:) = reshape(Q,n1,n2,n3);
if smooth
    fq = smoothdata2(fq,'gaussian',{(n1-1)/2, (n2-1)/2});
end
ub = asolve(fq);
Bt = ub(ysid(1):ysid(2),ysid(3):ysid(4),:);

% finish off rSVD, U*S*V' = Q*Q'*A
Bt = reshape(Bt,n1*n2,n3);
[Q1,R1] = qrw2(Bt,W); % Q1 is (n1*n2)-by-n3, R1 is n3-by-n3
[U1,S,V1] = svd(R1');
U = Q*U1;  % U is (n1*n2)-by-n3
V = Q1*V1; % V is (n1*n2)-by-n3

% truncate
S = S(1:k,1:k);
U = U(:,1:k);
V = V(:,1:k);

% return # of matvecs (PDE solves) used
mvs = 2*(k+p);
end

% simple weighted QR
function [Q,R] = qrw2(A,W)
V = sqrt(W);
[U,R] = qr(V.*A,'econ');
Q = (1./V) .* U;
end