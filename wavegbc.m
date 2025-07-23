function g = wavegbc(x, t, y, s, c)
%% WAVEGBC
%   Homogeneous Green's function for 1D wave equation in [0,1]x[0,1] domain
%   with constant wave speed and homogeneous initial/boundary values:
%       u_tt - c^2 * u_xx = f.
%
%   Input:
%   * x : array of x-values, from NDGRID
%   * t : array of t-values, from NDGRID
%   * y : array of y-values, from NDGRID
%   * s : array of s-values, from NDGRID
%   * (optional) c : positive wave speed = 1
%
%   Output:
%   * g : array of values for Green's function, of same size as inputs

arguments
    x double
    t double
    y double
    s double
    c (1,1) double = 1
end

assert(c > 0);

g = zeros(size(x));
dx1 = abs(y-x);
dx2 = abs(1-y-x);
dt = c*(t-s);

for j = 0:c/2
    g(dt > dx1 & dt <= 1 - dx2) = 1/(2*c);
    g(dt > dx2 + 1 & dt <= 2 - dx1) = -1/(2*c);
    dt = dt - 2;
end
end


