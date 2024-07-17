function pass = test_rsvd(tol)
% Test rSVD algorithm with weighted and unweighted inner product.

if nargin == 0
    tol = 1e-4;
end

k = 10;
p = 5;
q = 0;

M1 = randn(100,100);
M2 = randn(100,100);
[~, Wl] = legpts(100);
[~, Wr] = legpts(100, [0 0.5]);
Wl = Wl';
Wr = Wr';
s = 1./(1:100).^4;
S = diag(s);

% unweighted X
[U, ~] = qr(M1);
[V, ~] = qr(M2);
X = U*S*V';

% weighted X
[Uw, ~] = qrw(M1,Wl);
[Vw, ~] = qrw(M2,Wr);
Xw = Uw*S*Vw';

% unweighted test
[U1,S1,V1] = rsvd(X,k,p,q);
U1 = U1(:,1:k);
S1 = S1(1:k,1:k);
V1 = V1(:,1:k);
X1 = U1*S1*V1';

% check orthonormality
pass(1) = (norm(U1'*U1 - eye(k)) < tol);
pass(2) = (norm(V1'*V1 - eye(k)) < tol);

% spectral error
errspec = norm(X1-X);
pass(3) = (errspec - s(k+1) > 0); % Eckart-Young-Mirsky
pass(4) = ((errspec - s(k+1)) / norm(X) < tol);

% frobenius error
errfrob = norm(X1-X, 'fro');
pass(5) = (errfrob - norm(s(k+1:end)) > 0); % Eckart-Young-Mirsky
pass(6) = ((errfrob - norm(s(k+1:end))) / norm(X, 'fro') < tol);

% weighted test
[U2,S2,V2] = rsvd(Xw,k,p,q,Wl,Wr);
U2 = U2(:,1:k);
S2 = S2(1:k,1:k);
V2 = V2(:,1:k);
X2 = U2*S2*V2';

% check orthonormality
pass(7) = (norm(U2'*(Wl.*U2) - eye(k)) < tol);
pass(8) = (norm(V2'*(Wr.*V2) - eye(k)) < tol);

% spectral error
errspec  = norm(sqrt(Wl) .* (X2-Xw) .* sqrt(Wr)');
pass(9)  = (errspec - s(k+1) > 0); % Eckart-Young-Mirsky
pass(10) = ((errspec - s(k+1)) / norm(sqrt(Wl) .* Xw .* sqrt(Wr)') < tol);

% frobenius error
errfrob  = sqrt(Wl' * (X2-Xw).^2 * Wr);
pass(11) = (errfrob - norm(s(k+1:end)) > 0); % Eckart-Young-Mirsky
pass(12) = ((errfrob - norm(s(k+1:end))) / sqrt(Wl'*Xw.^2*Wr) < tol);
end