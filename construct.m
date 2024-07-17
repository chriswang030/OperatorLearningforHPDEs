function [tree, data, reds, redvol, level, lent, lend, errabs, errrel, errg, mvs] = ...
    construct(G, k, p, q, tol, C, N, maxlevel, doerrrel, level, reds, verbose, worker)
%% CONSTRUCT
%   Construct the approximate Green's function in a hierarchical,
%   memory-efficient data structure (Algorithm 2).
%
%   Input:
%   * G      : @(x,t,y,s) -> double, actual Green's function
%   * k      : target rank
%   * p      : oversampling
%   * q      : iterations
%   * tol    : error tolerance
%   * (optional) C        : positive scalar value, or vector for different 
%                           C per level, for rank test = 4
%   * (optional) N        : discretization dimension = 15
%   * (optional) maxlevel : hierarchy level cap = 10
%   * (optional) doerrrel : compute relative error = true
%   * (optional) level    : starting level = 0
%   * (optional) reds     : starting coord = [0 0 0 0]
%   * (optional) verbose  : verbosity (values can be 0, 1 or 2) = 1
%   * (optional) worker   : # worker for CONSTRUCTPAR
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
%   * mvs    : number of matvecs used

arguments
    G        function_handle
    k        (1,1) double
    p        (1,1) double
    q        (1,1) double
    tol      (1,1) double
    C        (:,1) double = 4
    N        (1,1) double = 15      
    maxlevel (1,1) double = 10
    doerrrel logical = true
    level    (1,1) double = 0
    reds     (:,4) double = [0 0 0 0]
    verbose  int8 = 1
    worker   double = 0
end

mvs = 0;                        % total matvecs used
err = 0;                        % track squared L2 error
tree = zeros(1000,1);           % tree structure of partition hierarchy
data = cell(100,3);             % store SVDs of low-rank blocks
lent = size(reds,1);            % current length of tree
lend = 0;                       % current length of data
reds = [reds (1:lent)'];        % update reds with tree indices
redvol = 2^(-4*level) * lent;   % calculate redvol from # of reds

% set up base quadrature rule
[bb,wb] = legpts(N, [0 1]);
[xb,tb,yb,sb] = ndgrid(bb,bb,bb,bb);
Wb = reshape(wb'*wb, N^2, 1);

% partition loop
while level <= maxlevel && redvol > tol^2
    errlevel = 0;                       % track error this level
    newreds = zeros(16*length(reds),5); % preallocate reds for next level
    chlen = 2^(-level-1);               % side length of child nodes
    lenr = size(reds,1);                % number of red blocks this level
    ic = min(level+1, length(C));       % index for C if vector
    
    % loop over reds at current level
    for i = 1:lenr
        coord = reds(i,1:4);            % get coords of a red block
        tree(reds(i,5)) = lent + 1;     % update tree with index of its children

        if verbose > 1
            fprintf('%d: Red block at (%.4f, %.4f, %.4f, %.4f)\n', i, coord);
        end

        % compute coordinates of children
        children = zeros(16, 4);
        j = 1;
        for j1 = 0:1
            for j2 = 0:1
                for j3 = 0:1
                    for j4 = 0:1
                        cc = [j1*chlen j2*chlen j3*chlen j4*chlen];
                        children(j,:) = coord + cc;
                        j = j+1;
                    end
                end
            end
        end

        % check if children are green or red
        for j = 1:16
            ch = children(j,:);     % get coords of child

            if verbose > 1
                fprintf('\t%d.%-2d Child at (%.4f, %.4f, %.4f, %.4f): ', i, j, ch);
            end

            % skip if we know it's zero from domain
            if ch(4) >= ch(2)+chlen
                if verbose > 1 fprintf('zero\n'); end
                continue;
            end

            % scale quadrature to domain
            xx = chlen * xb + ch(1);
            tt = chlen * tb + ch(2);
            yy = chlen * yb + ch(3);
            ss = chlen * sb + ch(4);

            % form matrix from G
            Gch = reshape(G(xx,tt,yy,ss), N^2, N^2);
            W = chlen^2 * Wb;

            % rank detection, save some matvecs if no power iteration
            if q > 0
                [high, mvsrd] = rankd(Gch, k, p, q, tol, C(ic), W, W);
                mvs = mvs + mvsrd;
            else
                [U, S, V, mvssvd] = rsvd(Gch, k, p, q, W, W);
                high = (S(k,k) > C(ic) * tol * S(1,1));
                mvs = mvs + mvssvd;
            end

            % if high rank, color red
            if high
                if verbose > 1 fprintf('red\n'); end
                newreds(16*(i-1)+j,:) = [ch lent+j];

            % if low rank, approximate with rSVD
            else
                if verbose > 1 fprintf('green\n'); end

                % do rSVD if not done already
                if q > 0
                    [U, S, V, mvssvd] = rsvd(Gch, k, p, q, W, W);
                    mvs = mvs + mvssvd;
                end

                errch = W' * (Gch-U*S*V').^2 * W;
                errlevel = errlevel + errch;

                % only save if its not numerically zero, otherwise
                % tree(lent+j) = 0 and block reads as zero
                if S(1,1) > 1e-16
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
    err = err + errlevel;
    if verbose
        fprintf("-------------------------\nLevel %d (worker: %d)\n\n" + ...
            "# red blocks:  %d\n" + ...
            "Total red vol: %f\n" + ...
            "Sq. error this level: %.4e\n" + ...
            "Sq. error cumulative: %.4e\n", ...
            level, worker, lenr, redvol, errlevel, err);
    end

    % set next level reds
    level = level + 1;
    reds = newreds(any(newreds,2),:);
    redvol = size(reds,1) * 2^(-4*level);
end

errg = sqrt(err);

% compute errors from remaining red blocks
chlen = 2^(-level);
lenr = size(reds,1);
for i = 1:lenr
    ch = reds(i,1:4);   % get coords of red block

    % scale quadrature to domain
    xx = chlen * xb + ch(1);
    tt = chlen * tb + ch(2);
    yy = chlen * yb + ch(3);
    ss = chlen * sb + ch(4);
    
    % form matrix from G
    Gch = reshape(G(xx,tt,yy,ss), N^2, N^2);
    W = chlen^2 * Wb;

    % add to error
    err = err + W'*Gch.^2*W;
end

% compute final absolute error
errabs = sqrt(err);
errrel = -1;

% compute final relative error
if doerrrel
    [xx, wx] = legpts(100, [0 1]);
    [xx,tt,yy,ss] = ndgrid(xx,xx,xx,xx);
    Gch = reshape(G(xx,tt,yy,ss), 100^2, 100^2);
    W = reshape(wx'*wx, 100^2, 1);
    Gnorm = sqrt(W' * Gch.^2 * W);
    errrel = errabs / Gnorm;
end

end