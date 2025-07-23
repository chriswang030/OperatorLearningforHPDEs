% metadata
mach = 'local';
time = datetime('now');
timestr = string(time, 'MM-dd-yyyy_HH-mm-ss');

% save
dosave = 0;

% algorithm parameters
k = 2;          % target rank for rSVD
p = 1;          % oversampling parameter for rSVD
tol = 0.01;     % error tolerance for partition error and rank detection
Nx = 65;        % number of spatial gridpoints
Nt = 129;       % number of temporal gridpoints (consider CFL condition)
C = 1;          % constant factor used in rank detection
maxlevel = 2;   % cap on number of partition levels
smooth = 0;     % do not smooth random test functions
verbose = 1;

% plot settings
y = 0.8;        % y-value in slice of Green's function plotted
s = 0.1;        % s-value in slice of Green's function plotted

% Green's function
c = 2;          % wave speed
G = @(x,t,y,s) wavegbc(x,t,y,s,c);

% set solver and adjoint solver for speed-2 wave eqution
a = 4*ones(Nx,1);
solve  = @(f) solve_uw1(f,a);
asolve = @(f) flip(solve_uw1(flip(f,2),a),2);

% run main algorithm
[tree, data, reds, redvol, height, ~, ~, errabs, errrel, errg, mvs] = ...
    constructbb(solve,asolve,k,p,tol,Nx,Nt,C,maxlevel,smooth,verbose,G=G);

% print summary
fprintf('\n\n------------ SUMMARY ------------\n');
fprintf('# matvecs: %d\n', mvs);
fprintf('Final reds volume: %f\n', redvol);
fprintf('Final abs. error:  %f\n', errabs);
fprintf('Final rel. error:  %f\n', errrel);

% plot
surfplotbb(G,tree,data,Nx,Nt,y,s);

% save everything
if dosave
    save(sprintf('data/%s_%s.mat',mach,timestr),'-v7.3');
    savefig(gcf,sprintf('data/%s_%s.fig',mach,timestr),'compact');
end