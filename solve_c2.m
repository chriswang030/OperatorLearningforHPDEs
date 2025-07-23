function u = solve_c2(f,a,b,c,u0,du0,init)
%% SOLVE_C2
%   Finite difference solver for the 1-dimensional hyperbolic PDE
%       u_tt - a(x)*u_xx + b(x)*u_x + c(x)*u = f
%       u(x,0) = u0(x)
%       u_t(x,0) = du0(x)
%   in the domain [0,1]x[0,1] with homogeneous boundary conditions using 
%   centered-time/centered-space finite differences.
%
%   Input:
%   * f : Nx x Nt x m representing m forcing terms on Nx x Nt grid
%   * a : Nx x 1 time-independent wavespeed coefficient 
%   * (optional) b    : time-independent coefficient = 0
%   * (optional) c    : time-independent coefficient = 0
%   * (optional) u0   : initial value u(x,0) = 0
%   * (optional) du0  : initial value u_t(x,0) = 0
%   * (optional) init : whether to initiate time-step with one Euler step 
%                       or with a ghost point = "euler"
%
%   Output:
%   * u : Nx x Nt x m, solution(s) to Cauchy problem

arguments
    f    (:,:,:) double
    a    (:,1) double
    b    (:,1) double = zeros(size(f,1),1)
    c    (:,1) double = zeros(size(f,1),1)
    u0   (:,1) double = zeros(size(f,1),1)
    du0  (:,1) double = zeros(size(f,1),1)
    init string = "euler"
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
u(:,1,:) = u0;
in = 2:Nx-1;

% initialize with ghost point or one forward Euler step
if init == "ghost"
    u(in,2,:) = u(in,1,:) + ht*du0(in,1,:) + 0.5*(ht/hx*a(1))^2 ...
        * (u(in-1,1,:) - 2*u(in,1,:) + u(in+1,1,:));
elseif init == "euler"
    u(:,2,:) = u0 + ht*du0;
end

% iterate
for k = 2:Nt-1
    % divdiff = u(in-1,k,:) - 2*u(in,k,:) + u(in+1,k,:);
    % stencil = eno(divdiff);
    u(in,k+1,:) = 2*u(in,k,:) - u(in,k-1,:) + ht^2 * f(in,k,:) ...
        + (ht/hx)^2 * a(in) .* (u(in-1,k,:) - 2*u(in,k,:) + u(in+1,k,:)) ...
        - ht^2/(2*hx) * b(in) .* (u(in+1,k,:) - u(in-1,k,:)) ...
        - ht^2/2 * c(in) .* (u(in+1,k,:) + u(in-1,k,:));
end
end