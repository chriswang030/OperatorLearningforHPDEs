function W = trap2(n1,n2,hx,ht)
%% TRAP2
%   Compute quadrature weights for the composite trapezoidal rule on a
%   2-dimensional domain.
%
%   Input:
%   * n1 : number of gridpoints in x-direction
%   * n2 : number of gridpoints in t-direction
%   * hx : mesh size in x direction
%   * ht : mesh size in t direction
%
%   Output:
%   * W : (n1*n2) column vector of weights

w = ones(n1,n2);
w(2:end-1,:) = w(2:end-1,:) + 1;
w(:,2:end-1) = w(:,2:end-1) + 1;
w(2:end-1,2:end-1) = w(2:end-1,2:end-1) + 1;
W = 0.25*hx*ht*reshape(w,n1*n2,1);
end