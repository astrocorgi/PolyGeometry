%function [nodes,bounds,header,facet] = mesh_bound_init(a,b,e,f,w);
clear;clf;
a = 10000; % X length
b = 1000;  % Y length
e = 100;   % element size [m]
w = .5;    % thickness of water, as percent of Y length
%--------------------------------------------------------------------------
%                BC 32 : snow (erosion)
%           *       *       *       *
%          __*____*____*  *   *    *    *
%         |             ____*_  *   *
% BC 1:  O|                      ___ *    *
% roller  |              side 3      \  *
%        O| side 1                      \ *_________________ water level
%         |         side 2                \  <--   BC 32:
%        O|_______________________|________| <---- horizontal hydrostat.
%          O     O     O     O     ^^^^^^ 
%                                  ||||||                    
%               BC 16: roller       BC 32: vertical hydrostat.
%--------------------------------------------------------------------------
ny = floor(b/e);
nx = floor(a/e);
ns = (ny-1) + (nx-2); %taking out (0,b) and (0,0)&(a,0)

S1      = NaN*zeros(ny,3);
S2      = NaN*zeros(nx+1,3);
S3      = NaN*zeros(ns,3);

delndx      = [];
segmentchk  = 5;
nnx         = 2;

% side 1 -- contains (0,b) 
for n = 1:ny
    S1(n,1) = 0.0;
    S1(n,2) = 0.0 + e*n; %starts above 0, ends at b
    S1(n,3) = 1;
end;

% side 2 -- contains (0,0)and (a,0) <-- nx+1
for n = 1:nx+1
    S2(n,1) = 0.0 + e*(n-1);
    S2(n,2) = 0.0;
    S2(n,3) = 16;
end;

% side 3 -- elliptical 
for n = 1:ny-1 % loop S1 excludes (0,b)
    S3(n,1) = sqrt(a^2*(1-(S1(n,2)/b)^2));   % x = sqrt(a^2 * (1 - (y/b)^2)
    S3(n,2) =              S1(n,2)       ;   % y = y
    S3(n,3) = 32;
end;

for n = ny:ns+1 % loop S2 excludes (0,0) and (a,0)
    S3(n,1) =              S2(nnx,1)       ; % x = x
    S3(n,2) = sqrt(b^2*(1-(S2(nnx,1)/a)^2)); % y = sqrt(b^2 * (1 - (x/a)^2)
    nnx = nnx + 1;
    S3(n,3) = 32;
end;

S3 = sortrows(S3,1);

% smoothing
for smoothing = 1:segmentchk;

    for n = 2:ns
        distchk = sqrt((S3(n,1) - S3(n-1,1))^2 + (S3(n,2) - S3(n-1,2))^2);
        cond = distchk < e; 
        if cond; delndx(end+1) = n; end;
    end;

    for n = 1:length(delndx);
        nn = delndx(n);
    
        newx = 0.5*(S3(nn+1,1) + S3(nn-1,1));
        newy = 0.5*(S3(nn+1,2) + S3(nn-1,2));
        S3(nn,1) = newx;
        S3(nn,2) = newy;
    end;

end;

s1 = size(S1); s2 = size(S2); s3 = size(S3);
nnodes = s1(1) + s2(1) + s3(1);
nodecol = [0:1:nnodes-1]'; f2 = [0; nodecol(1:end-1)]; f3 = [nnodes-1; nodecol(2:end)];
%--------------------------------------------------------------------------
% max thickness of ice that can float in water of depth b*w:
hi = b*w*1030/917;
% node on side 3 where ice is <= hi:
nd3 = min(find(S3(:,2) <= hi));
% x-position of that max-thickness node:
xgl = S3(nd3,1);
% all the side 2 nodes >= that x-position:
nd2 = find(S2(:,1) >= xgl);
% BC flag assigned to side 2 nodes which should be floating:
S2(nd2,3) = 32;
% re-order side 2 nodes in reverse
tmp = S2([end:-1:1],:); S2 = tmp; 
% clock-wise ordering of all nodes:
S = [S1 ; S3 ; S2];

% put floating part below water?
S(:,2) = S(:,2) - hi;

head1 = [nnodes 2 0 0]; head2 = [nnodes 1];
nodes = [nodecol S(:,1) S(:,2)];
facet = [nodecol f2 f3 S(:,3)];


