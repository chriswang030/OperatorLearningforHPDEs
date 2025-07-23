function u = solve_uw1(f,a,u0,du0)
%% SOLVE_UW1
%   Finite difference solver for the 1-dimensional hyperbolic PDE
%       u_tt - a(x)*u_xx = f
%       u(x,0) = u0(x)
%       u_t(x,0) = du0(x)
%   in the domain [0,1]x[0,1] with homogeneous boundary conditions using
%   the UW1 scheme of Banks & Henshaw.
%   (J. W. Banks, W. D. Henshaw, J. Computational Physics, 2012)
%
%   Input:
%   * f : Nx x Nt x m representing m forcing terms on Nx x Nt grid
%   * a : Nx x 1 time-independent wavespeed coefficient 
%   * (optional) u0   : initial value u(x,0) = 0
%   * (optional) du0  : initial value u_t(x,0) = 0
%
%   Output:
%   * u : Nx x Nt x m, solution(s) to Cauchy problem

arguments
    f    (:,:,:) double
    a    (:,1) double
    u0   (:,1) double = zeros(size(f,1),1)
    du0  (:,1) double = zeros(size(f,1),1)
end

% check that initial conditions satisfy homogeneous BCs
assert(abs(u0(1))    < 100*eps);
assert(abs(du0(1))   < 100*eps);
assert(abs(u0(end))  < 100*eps);
assert(abs(du0(end)) < 100*eps);

% setup
[Nx,Nt,m] = size(f); % allow parallel solves for multiple f
hx = 1/(Nx-1);
ht = 1/(Nt-1);
u0 = repmat(u0,1,1,m);
du0 = repmat(du0,1,1,m);
u = zeros(Nx,Nt,m);
v = zeros(Nx,Nt,m);
u(:,1,:) = u0;
v(:,1,:) = du0;
in = 2:Nx-1;

% iterate
for k = 1:Nt-1
    d2u = 1/hx^2 * (u(in+1,k,:) - 2*u(in,k,:) + u(in-1,k,:));
    d2v = 1/hx^2 * (v(in+1,k,:) - 2*v(in,k,:) + v(in-1,k,:));
    u(in,k+1,:) = u(in,k,:) + ht*v(in,k,:) + ht^2/2 * a(in) .* d2u ...
        + hx*ht^2/4 * sqrt(a(in)) .* d2v;
    v(in,k+1,:) = v(in,k,:) + ht * a(in) .* d2u ...
        + hx*ht/2 * sqrt(a(in)) .* d2v + ht * f(in,k,:);
end
end