function Vq = hmtimes(tree,data,N,f,Xq,Tq,ivp)
%% HMTIMES
%   Computes operator-function product of the approximate Green's function
%   computed from CONSTRUCT or CONSTRUCTPAR against the forcing term f,
%   evaluated at the query points Xq, Tq in space and time, respectively.
%
%   Input:
%   * tree : integer array with tree structure
%   * data : cell array with stored SVDs of low-rank blocks
%   * N    : side length of the matrix representing each block, so that
%            cells of DATA store SVDs of N^2 x N^2 matrices
%   * f    : @(x,t) -> double, forcing term
%   * Xq   : vector of query points in space
%   * Tq   : vector of query points in time
%   * (optional) ivp : whether f represents an initial condition, for 
%                      speedup = false
%
%   Output:
%   * Vq   : matrix of values from the operator-function product 

arguments
    tree (:,1) double
    data (:,3) cell
    N    (1,1) double
    f    function_handle
    Xq   (:,1) double
    Tq   (:,1) double
    ivp  logical = 0
end

% setup
Nx = length(Xq);
Nt = length(Tq);
Vq = zeros(Nx,Nt);
[xxg,ttg] = meshgrid(Xq,Tq);

% compute base Legendre nodes
[xl,wb] = legpts(N,[0 1]);
[xb,tb] = ndgrid(xl,xl);
Wb = reshape(wb'*wb,N^2,1);

% initialize stack
stack = [0 1 0 1 0 1 0 1 1]; % store domain intervals + tree index
is = 1;                      % stack index

% stack loop
while is > 0
    % retrieve next item from stack
    dom = stack(is,1:8);    % domain intervals
    it  = stack(is,9);      % tree index
    len = dom(2) - dom(1);  % side length of block

    % update stack size
    is  = is - 1;
    
    % block is green or zero
    if it > length(tree) || tree(it) <= 0
        % retrieve block from f
        yy = xb * len + dom(5);
        ss = tb * len + dom(7);
        Fch = reshape(f(yy,ss),N^2,1);

        % if numerically zero, skip the computations
        if max(abs(Fch),[],'all') < 1e-8
            continue;
        end

        % get query points
        [xxq,ttq] = ndgrid(Xq(dom(1) <= Xq & Xq < dom(2)), ...
            Tq(dom(3) <= Tq & Tq < dom(4)));
    
        % if block is green, retrieve approximate Green's function
        if it <= length(tree) && tree(it) < 0
            id = -tree(it);    % green IDs were stored as negative integers
            U = data{id,1};
            S = data{id,2};
            V = data{id,3};
           
            % get Legendre nodes and weights in this domain
            xx = xb * len + dom(1);
            tt = tb * len + dom(3);
            W = len^2 * Wb;

            % do blockwise operator-function product
            Vch = reshape(U*(S.*(V'*(W.*Fch))),N,N);
            Vch = interpn(xx,tt,Vch,xxq,ttq,'linear');

            % add to global solution
            Vch = Vch'; % transpose for linear indexing
            ind = dom(1) <= xxg & xxg < dom(2) & dom(3) <= ttg & ttg < dom(4);
            Vq(ind) = Vq(ind) + Vch(:);
        end

    % block is red
    else
        % setup for child domains
        lench = len/2;
        dom(2:2:8) = dom(2:2:8) - lench;
    
        % add children to the stack
        j = 0;
        for j1 = 0:1
            for j2 = 0:1
                for j3 = 0:1
                    for j4 = 0:1
                        cc = [j1 j1 j2 j2 j3 j3 j4 j4] * lench;
                        domch = dom + cc;

                        if ~ivp || domch(7) <= 1e-12
                            is = is+1;
                            stack(is,:) = [domch tree(it)+j];
                        end

                        j = j+1;
                    end
                end
            end
        end
    end
end

% extrapolate some gaps between blocks
% Vq = inpaint_nans(Vq);
end