function surfplot(G, tree, data, N, d, y, s)
%% SURFPLOT
%   Plot pointwise approximation against true Green's function and error
%   plot, at specified (y,s)-slice.
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
    N    (1,1) double
    d    (1,1) double
    y    (1,1) double
    s    (1,1) double
end

% initialize plot functions
F = zeros(d,d);
E = zeros(d,d);
xxu = linspace(0,1,d);
ttu = linspace(0,1,d);
[xxg,ttg] = meshgrid(xxu,ttu);
rects = zeros(100,4);
ir = 0;

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
        Gapprox = 0;
        err = 0;

        % if block is green, retrieve approximate Green's function
        if it <= length(tree) && tree(it) < 0
            id = -tree(it);    % green IDs were stored as negative integers
            U = data{id,1};
            S = data{id,2};
            V = data{id,3};

            % get query points for interpn
            [xxq,ttq,yyq,ssq] = ...
                ndgrid(xxu(dom(1) <= xxu & xxu < dom(2)), ...
                       ttu(dom(3) <= ttu & ttu < dom(4)), ...
                       y, s);
           
            % get legpts from how data is actually stored
            xx = legpts(N, dom(1:2));
            tt = legpts(N, dom(3:4));
            yy = legpts(N, dom(5:6));
            ss = legpts(N, dom(7:8));
            [xx,tt,yy,ss] = ndgrid(xx,tt,yy,ss);

            % interpolate into uniform grid
            Gapprox = reshape(U*diag(S)*V', N, N, N, N);
            Gapprox = interpn(xx,tt,yy,ss,Gapprox,xxq,ttq,yyq,ssq);

            % discretize actual Green's function and compute error
            Gch = G(xxq,ttq,y,s);
            err = abs(Gch-Gapprox);
        end

        % save plots of block
        Gapprox = Gapprox'; % transpose for linear indexing
        err = err';         % transpose for linear indexing
        F(dom(1) <= xxg & xxg < dom(2) & dom(3) <= ttg & ttg < dom(4)) ...
            = Gapprox(:);
        E(dom(1) <= xxg & xxg < dom(2) & dom(3) <= ttg & ttg < dom(4)) ...
            = err(:);

        % save domain for plot
        ir = ir + 1;
        rects(ir,:) = [dom(1) dom(3) len len];

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
                        cc = [j1*lench j1*lench j2*lench j2*lench ...
                              j3*lench j3*lench j4*lench j4*lench];
                        domch = dom + cc;

                        if domch(5) <= y && y <= domch(6) ...
                            && domch(7) <= s && s <= domch(8)
                            is = is + 1;
                            stack(is,:) = [domch tree(it)+j];
                        end

                        j = j+1;
                    end
                end
            end
        end
    end
end

% stitch blocks together for NaNs from interpn 
F = inpaint_nans(F);
E = inpaint_nans(E);

% plot setup
[xxu,ttu] = meshgrid(xxu,ttu);
tiledlayout(1,3);

% plot approximate Green's function with partition blocks
h(1) = nexttile;
contourf(xxu, ttu, F, EdgeColor='none');
for i = 1:ir
    rectangle(Position=rects(i,:), LineWidth=1);
end
title('$\tilde G$', interpreter='latex');
xlabel('$x$', interpreter='latex');
yl = ylabel('$t$', interpreter='latex', rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% plot true Green's function
h(2) = nexttile;
contourf(xxu, ttu, G(xxu,ttu,y,s), EdgeColor='none');
title('$G$', interpreter='latex');
xlabel('$x$', interpreter='latex');
yl = ylabel('$t$', interpreter='latex', rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% plot error
h(3) = nexttile;
surf(xxu, ttu, E,  EdgeColor='none');
view(2);
title('$\|G-\tilde G\|_{L^2}$', interpreter='latex');
xlabel('$x$', interpreter='latex');
yl = ylabel('$t$', interpreter='latex', rotation=0);
yl.Position(1) = yl.Position(1) - 0.05;

% set z-axes
zl = zlim(h(2));
zlim(h(1),zl);

% colorbars
[s1,l1] = bounds(G(xxu,ttu,y,s),'all');
[s2,l2] = bounds(F,'all');
s = min(s1,s2);
l = max(l1,l2);
clim(h(1),[s l]);
clim(h(2),[s l]);
cb1 = colorbar(h(1));
cb1.Layout.Tile = 'west';

l3 = max(E,[],'all');
clim(h(3),[0 l3*0.05]);
cb2 = colorbar(h(3));
cb2.Layout.Tile = 'east';

% fontsize
fontsize(24,'points');
end