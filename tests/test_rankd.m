function pass = test_rankd(tol)
% Test rank detection algorithm

if nargin == 0
    tol = 1e-16;
end

k = 20;
p = 5;
q = 1;

rng(1);
M1 = randn(100,100);
M2 = randn(100,100);
Slow  = diag([ones(1,k-1) tol*ones(1,100-k+1)]);    % k-1 to give some leeway
Shigh = diag(linspace(1,10*tol,100));

[U, ~] = qr(M1);
[V, ~] = qr(M2);
Xlow  = U*Slow*V';  % low-rank X
Xhigh = U*Shigh*V'; % high-rank X

% test low-rank
[high1, ~] = rankd(Xlow,k,p,q,tol);
pass(1) = ~high1;

% test high-rank
[high2, ~] = rankd(Xhigh,k,p,q,tol);
pass(2) = high2;

end