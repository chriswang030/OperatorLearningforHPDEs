function [tree, data, redvol, height, errabs, errrel, mvs] = ...
    constructpar(G, k, p, q, tol, C, N, maxlevel, prelevel, verbose)
%% CONSTRUCTPAR
%   Use parallel workers to construct the approximate Green's function in a 
%   hierarchical, memory-efficient data structure (Algorithm 2).
%
%   Input:
%   * G      : @(x,t,y,s) -> double, actual Green's function
%   * k      : target rank
%   * p      : oversampling
%   * q      : power iterations
%   * tol    : error tolerance
%   * (optional) C        : positive scalar value, or vector for different 
%                           C per level, for rank test = 4
%   * (optional) N        : discretization dimension = 15
%   * (optional) maxlevel : hierarchy level cap = 10
%   * (optional) prelevel : levels to perform preparallelization = 0
%   * (optional) verbose  : verbosity (values can be 0, 1 or 2) = 0
% 
%   Output:
%   * tree   : integer array containing tree structure
%   * data   : cell array containing low-rank representations of blocks
%   * redvol : final volume of red blocks
%   * height : final height of tree
%   * errabs : absolute approximation error in the L2 norm
%   * errrel : relative approximation error in the L2 norm
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
    prelevel (1,1) double = 0
    verbose  (1,1) int8 = 0
end

% run pre-parallel
[tree, data, reds, ~, level, lent, lend, ~, ~, errg, mvs] = ...
    construct(G,k,p,q,tol,C,N,prelevel,verbose);

% update C if needed
pC = C(end);
if length(C) > prelevel+1
    pC = C(prelevel+2:end);
end

% initialize structures for parfor outputs
lenr = size(reds,1);
ptree = cell(lenr,1);
pdata = cell(lenr,1);
predvol = zeros(lenr,1);
pheight = zeros(lenr,1);
plent = zeros(lenr,1);
plend = zeros(lenr,1);
perrabs = zeros(lenr,1);
pmvs = zeros(lenr,1);

% parallel loop
parfor i = 1:lenr
%for i = 1:lenr
    coord = reds(i,1:4);            % get coords of a red block
    [pt, pd, ~, pr, ph, plt, pld, pe, ~, ~, pm] = ...
        construct(G,k,p,q,tol/lenr,pC*lenr,N,maxlevel,false,level,coord,verbose,i);
    ptree{i} = pt;
    pdata{i} = pd;
    predvol(i) = pr;
    pheight(i) = ph;
    plent(i) = plt;
    plend(i) = pld;
    perrabs(i) = pe;
    pmvs(i) = pm;
end

% reconstruct tree and data
if verbose
    fprintf("\nCompleted parallel runs, reconstructing data structures... ");
end

for i = 1:lenr
    tree(reds(i,5)) = lent + 1;
    pt = ptree{i};                              % copy tree from one worker
    pt = [pt ; zeros(plent(i)-length(pt),1)];   % extend by zero if needed
    pt = pt(2:end);                             % remove first element for splicing
    pt(pt > 0) = pt(pt > 0) + lent - 1;         % shift tree indices for positive entries
    pt(pt < 0) = pt(pt < 0) - lend;             % shift data indices for negative entries
    tree(lent+1:lent+plent(i)-1) = pt;          % copy over pt
    data(lend+1:lend+plend(i),:) = pdata{i};    % copy over pd
    lent = lent + plent(i) - 1;                 % update tree length
    lend = lend + plend(i);                     % update data length
end

if verbose
    fprintf("done\n");
end

% reconstruct other quantities
redvol = sum(predvol);
height = max(pheight);
errabs = sqrt(errg^2 + sum(perrabs.^2));
mvs = mvs + sum(pmvs);

% compute final relative error
[xx, wx] = legpts(100, [0 1]);
[xx,tt,yy,ss] = ndgrid(xx,xx,xx,xx);
Gch = reshape(G(xx,tt,yy,ss), 100^2, 100^2);
W = reshape(wx'*wx, 100^2, 1);
Gnorm = sqrt(W' * Gch.^2 * W);
errrel = errabs / Gnorm;

end