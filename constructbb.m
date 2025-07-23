function [tree, data, reds, redvol, level, lent, lend, errabs, errrel, errg, mvs] = ...
    constructbb(solve, asolve, k, p, tol, Nx, Nt, C, maxlevel, smooth, verbose, paropts, erropts)
%% CONSTRUCTBB
%   Construct the approximate Green's function in a hierarchical,
%   memory-efficient data structure (Algorithm 2), using a black-box
%   numerical solver to generate input-output data.
%
%   Input:
%   * solve  : @(f) -> double, forward numerical solver
%   * asolve : @(f) -> double, backward (adjoint) numerical solver
%   * k      : target rank
%   * p      : oversampling
%   * tol    : error tolerance
%   * Nx     : number of spatial gridpoints (must be of form 2^n+1)
%   * Nt     : number of temporal gridpoints (must be of form 2^n+1)
%   * (optional) C        : positive scalar value, or vector for different 
%                           C per level, for rank test = 1
%   * (optional) maxlevel : hierarchy level cap = 5
%   * (optional) smooth   : whether to artificially smooth random test
%                           functions = false
%   * (optional) verbose  : verbosity (values can be 0, 1 or 2) = 1
%   * (optional) worker   : # worker for CONSTRUCTPAR
%   * (optional, named) level : starting level = 0
%   * (optional, named) reds  : starting coord = [0 0 0 0]
%   * (optional, named) G     : @(x,t,y,s) -> double, actual Green's 
%                               function for error computation
% 
%   Output:
%   * tree   : integer array containing tree structure
%   * data   : cell array containing low-rank representations of blocks
%   * reds   : list of reds with tree indices
%   * redvol : final volume of red blocks
%   * level  : final height of TREE
%   * lent   : final length of TREE
%   * lend   : final length of DATA
%   * errabs : absolute approximation error in the L2 norm
%   * errrel : relative approximation error in the L2 norm, -1 if doerrrel = false
%   * errg   : error from only green blocks, not red blocks
%   * mvs    : row vector of number of matvecs used per level

arguments
    solve    function_handle
    asolve   function_handle
    k        (1,1) uint64
    p        (1,1) uint64
    tol      (1,1) double
    Nx       (1,1) double % must be 2^n+1
    Nt       (1,1) double % must be 2^n+1
    C        (1,1) double = 1
    maxlevel (1,1) double = 5
    smooth   logical = 0
    verbose  int8 = 1
    paropts.level   (1,1) double = 0
    paropts.reds    (:,4) double = [1 1 1 1]
    erropts.G       function_handle
end

assert(floor(log2(Nx-1)) == log2(Nx-1));
assert(floor(log2(Nt-1)) == log2(Nt-1));

mvs = 0;                         % total matvecs used
err = 0;                         % track squared L2 error
level = paropts.level;
tree = zeros(1000,1);            % tree structure of partition hierarchy
data = cell(100,3);              % store SVDs of low-rank blocks
lent = size(paropts.reds,1);     % current length of tree
lend = 0;                        % current length of data
reds = [paropts.reds (1:lent)']; % update reds with tree indices
redvol = 2^(-4*level) * lent;    % calculate redvol from # of reds
hx = 1/(Nx-1);                   % spatial grid length
ht = 1/(Nt-1);                   % temporal grid length

% precompute Green's function if G is given
doerr = ~isempty(fieldnames(erropts));
if doerr
    x = linspace(0,1,Nx);
    t = linspace(0,1,Nt);
    [xx,tt,yy,ss] = ndgrid(x,t,x,t);
    Gd = erropts.G(xx,tt,yy,ss);
end

% partition loop
while level <= maxlevel && redvol > tol^2
    errlevel = 0;                       % track error this level
    newreds = zeros(16*length(reds),5); % preallocate reds for next level
    chscale = 2^(-level-1);             % side length of child nodes
    lenr = size(reds,1);                % number of red blocks this level
    n1 = (Nx-1)*chscale+1;              % number of spatial gridpoints this level
    n2 = (Nt-1)*chscale+1;              % number of temporal gridpoints this level

    % 2D trapezoidal rule at this level
    W = trap2(n1,n2,hx,ht);

    % loop over reds at current level
    for i = 1:lenr
        cid = reds(i,1:4);          % get index coordinates of red block
        tree(reds(i,5)) = lent + 1; % update tree with index of its children

        if verbose > 1
            fprintf('%d: Red block at (%d, %d, %d, %d)\n', i, cid);
        end

        % coordinate intervals of children [x0 x1 t0 t1 y0 y1 s0 s1]
        children = zeros(16,8);
        base = [cid(1) cid(1) cid(2) cid(2) cid(3) cid(3) cid(4) cid(4)] ...
            + [0 1 0 0 0 1 0 0] * (n1-1) ...
            + [0 0 0 1 0 0 0 1] * (n2-1);
        j = 1;
        for j1 = 0:1
            for j2 = 0:1
                for j3 = 0:1
                    for j4 = 0:1
                        cc = [j1 j1 j2 j2 j3 j3 j4 j4] ...
                            .* ([1 1 0 0 1 1 0 0] * (n1-1) ...
                             +  [0 0 1 1 0 0 1 1] * (n2-1));
                        children(j,:) = base + cc;
                        j = j+1;
                    end
                end
            end
        end

        % check if children are green or red
        for j = 1:16
            ch = children(j,:); % get coordinate intervals of child

            % skip if we know it's zero from domain
            if ch(7) >= ch(4)
                continue;
            end

            % rSVD and use it for rank detection
            xtid = ch(1:4);
            ysid = ch(5:8);
            [U,S,V,mvssvd] = rsvdbb(solve,asolve,k,p,Nx,Nt,xtid,ysid,smooth);
            high = S(1,1) > 1e-12 && S(k,k) > C*tol*S(1,1);
            mvs = mvs + mvssvd;

            % if high rank, color red
            if high
                newreds(16*(i-1)+j,:) = [ch(1) ch(3) ch(5) ch(7) lent+j];

            % if low rank, save the rSVD output
            else
                % compute error from this level
                if doerr
                    Gdch = Gd(xtid(1):xtid(2),xtid(3):xtid(4), ...
                              ysid(1):ysid(2),ysid(3):ysid(4));
                    errch = W' * (reshape(Gdch,n1*n2,n1*n2)-U*S*V').^2 * W;
                    errlevel = errlevel + errch;
                end

                % only save if its not numerically zero, otherwise
                % tree(lent+j) = 0 and block reads as zero
                if S(1,1) > 1e-4
                    lend = lend + 1;        % update length of data array
                    tree(lent+j) = -lend;   % update tree with index of data
                    data{lend,1} = U;
                    data{lend,2} = diag(S); % only store the diagonal
                    data{lend,3} = V;
                end
            end
        end

        % update tree length
        lent = lent + 16;
    end

    % update errors
    if verbose
        fprintf("-------------------------\nLevel %d \n\n" + ...
            "# red blocks:  %d\n" + ...
            "Total red vol: %f\n" + ...
            "Total matvecs: %d\n", ...
            level, lenr, redvol, mvs);
    end
    if doerr
        err = err + errlevel;
        fprintf("Sq. error this level: %.4e\n" + ...
            "Sq. error cumulative: %.4e\n", ...
            errlevel, err);
    end

    % set next level reds
    level = level + 1;
    reds = newreds(any(newreds,2),:);
    redvol = size(reds,1) * 2^(-4*level);
end

if doerr
    errg = sqrt(err);
    
    % compute errors from remaining red blocks
    scale = 2^(-level);
    lenr = size(reds,1);
    n1 = (Nx-1)*scale+1;
    n2 = (Nt-1)*scale+1;
    W = trap2(n1,n2,hx,ht);
    
    for i = 1:lenr
        cid = reds(i,1:4);
        Gdch = Gd(cid(1):cid(1)+n1-1,cid(2):cid(2)+n2-1, ...
                  cid(3):cid(3)+n1-1,cid(4):cid(4)+n2-1);
        err = err + W' * reshape(Gdch,n1*n2,n1*n2).^2 * W;
    end
    
    % compute final absolute and relative error
    W = trap2(Nx,Nt,hx,ht);
    Gnorm = sqrt(W' * reshape(Gd,Nx*Nt,Nx*Nt).^2 * W);
    errabs = sqrt(err);
    errrel = errabs / Gnorm;
else
    errabs = -1;
    errrel = -1;
    errg = -1;
end
end

