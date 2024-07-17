% machine name
mach = 'local';

% algorithm parameters
k = 2;          % target rank for rSVD
p = 1;          % oversampling parameter for rSVD
q = 0;          % power iterations for rSVD
tol = 0.04;     % error tolerance for partition error and rank detection
N = 6;          % discretization dimension
C = 1;          % constant factor used in rank detection
maxlevel = 7;   % cap on number of partition levels
prelevel = 0;   % number of levels pre-parallelization for CONSTRUCTPAR
verbose = 1;    % setting for verbosity (can take values 0, 1, or 2)

% plot settings
d = 100;        % dimension of grid for plot
y = 0.8;        % y-value in slice of Green's function plotted
s = 0.1;        % s-value in slice of Green's function plotted

% Green's function
c = 4;          % wave speed
G = @(x,t,y,s) wavegbc(x,t,y,s,c);

% run main algorithm
[tree, data, redvol, height, errabs, errrel, mvs] = ...
    constructpar(G,k,p,q,tol,C,N,maxlevel,prelevel,verbose);

% print summary
fprintf('\n\n------------ SUMMARY ------------\n');
fprintf('# matvecs: %d\n', mvs);
fprintf('Final reds volume: %.4e\n', redvol);
fprintf('Final abs. error:  %.4e\n', errabs);
fprintf('Final rel. error:  %.4e\n', errrel);

% save everything
time = datetime('now');
timestr = string(time, 'MM-dd-yyyy_HH-mm-ss');
save(sprintf('data/%s_%s.mat', mach, timestr), '-nocompression', '-v7.3');
 
% plot and save
surfplot(G,tree,data,N,d,y,s);
savefig(gcf, sprintf('data/%s_%s.fig', mach, timestr), 'compact');