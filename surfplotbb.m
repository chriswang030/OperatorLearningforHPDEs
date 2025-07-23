function surfplotbb(G,tree,data,Nx,Nt,y,s)
%% SURFPLOTBB
%   Plot pointwise approximation against true Green's function and error
%   plot, at specified (y,s)-slice, using data from CONSTRUCTBB.
%
%   Input:
%   * G    : @(x,t,y,s) -> double, actual Green's function
%   * tree : integer array with tree structure
%   * data : cell array with stored SVDs of low-rank blocks
%   * N    : side length of the matrix representing each block, so that
%            cells of DATA store SVDs of N^2 x N^2 matrices
%   * d    : resolution of pointwise plot
%   * y    : y-slice for plot
%   * s    : s-slice for plot

arguments
    G    function_handle
    tree (:,1) double
    data (:,3) cell
    Nx   (1,1) double
    Nt   (1,1) double
    y    (1,1) double
    s    (1,1) double
end

yid = y*(Nx-1)+1;
sid = s*(Nt-1)+1;

% initialize plot functions
x = linspace(0,1,Nx);
t = linspace(0,1,Nt);
[xx,tt] = meshgrid(x,t);
Gtrue = G(xx,tt,y,s);    % true Green's function
Gapprox = zeros(Nx,Nt);  % approximate Green's function
rects = zeros(100,4);
rectid = 0;

% initialize stack
stack = [1 1 1 1 1 0]; % store index coordinates (1:4), tree index (5), level(6) 
stid = 1;              % stack index

% stack loop
while stid > 0
    % retrieve next item from stack
    cid = stack(stid,1:4);    % index coordinates
    tid = stack(stid,5);      % tree index
    level = stack(stid,6);
    scale = 2^(-level);
    n1 = (Nx-1)*scale+1;
    n2 = (Nt-1)*scale+1;

    % update stack size
    stid = stid - 1;
    
    % block is green or zero
    if tid > length(tree) || tree(tid) <= 0
        % if block is green, retrieve approximate Green's function
        if tid <= length(tree) && tree(tid) < 0
            id = -tree(tid);   % green IDs were stored as negative integers
            U = data{id,1};
            S = diag(data{id,2});
            V = data{id,3};

            % save to full plot
            F = reshape(U*S*V',n1,n2,n1,n2);
            Gapprox(cid(1):cid(1)+n1-1,cid(2):cid(2)+n2-1) ...
                = F(:,:,floor(yid)-cid(3)+1,floor(sid)-cid(4)+1);
        end

        % save domain for plot
        rectid = rectid + 1;
        rects(rectid,:) = [x(cid(1)) t(cid(2)) scale scale];

    % block is red
    else
        chx = (n1-1)/2;
        cht = (n2-1)/2;

        % add children to the stack
        j = 0;
        for j1 = 0:1
            for j2 = 0:1
                for j3 = 0:1
                    for j4 = 0:1
                        cc = [j1 j2 j3 j4] ...
                            .* ([1 0 1 0] * chx ...
                             +  [0 1 0 1] * cht);
                        cidch = cid + cc;

                        % add to stack if (y,s) is contained in child block
                        if cidch(3) <= yid && yid <= cidch(3)+chx ...
                            && cidch(4) <= sid && sid <= cidch(4)+cht

                            stid = stid + 1;
                            stack(stid,:) = [cidch tree(tid)+j level+1];
                        end

                        j = j+1;
                    end
                end
            end
        end
    end
end

% transpose to fit meshgrid
Gapprox = Gapprox';

% plot
tiledlayout(1,3);

% plot approximate Green's function with partition blocks
h(1) = nexttile;
contourf(xx,tt,Gapprox,EdgeColor='none');
for i = 1:rectid
    rectangle(Position=rects(i,:),LineWidth=1);
end
title(sprintf('$\\tilde G(x,t;%.1f,%.1f)$',y,s),interpreter='latex');
xlabel('$x$',interpreter='latex');
yl = ylabel('$t$',interpreter='latex',rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% plot true Green's function
h(2) = nexttile;
surf(xx,tt,Gtrue);
contourf(xx,tt,Gtrue,EdgeColor='none');
title(sprintf('$G(x,t;%.1f,%.1f$)',y,s),interpreter='latex');
xlabel('$x$',interpreter='latex');
yl = ylabel('$t$',interpreter='latex',rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% plot error
h(3) = nexttile;
E = abs(Gapprox - Gtrue);
surf(xx, tt, E,  EdgeColor='none');
view(2);
title('$|G-\tilde G|$',interpreter='latex');
xlabel('$x$',interpreter='latex');
yl = ylabel('$t$',interpreter='latex',rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% set z-axes
zl = zlim(h(2));
zlim(h(1),zl);

% colorbars
[s1,l1] = bounds(Gtrue,'all');
[s2,l2] = bounds(Gapprox,'all');
s = min(s1,s2);
l = max(l1,l2);
clim(h(1),[s l]);
clim(h(2),[s l]);
cb1 = colorbar(h(1));
cb1.Layout.Tile = 'west';

l3 = max(E,[],'all');
if l3 < 1e-8 l3 = 1; end
clim(h(3),[0 l3]);
cb2 = colorbar(h(3));
cb2.Layout.Tile = 'east';

% fontsize
fontsize(gcf,24,'points');
end