clear;clf;
% this script creates an ordered perimeter of an ellipse
% GEOMETRIC INPUTS:
L = 10000;              % length of the X-domain
grounded = 8000/10000;  % % of the basal side which is grounded
Hi0 = 1500;             % this is the LHS thickness
dhill = 0;              % this is the LHS thickness of a bed-topo-wedge
tongH = 500;            % this is the RHS vert. amt. of ice below water
Yres = 50;              % vertical resolution
Xres = 50;              % horizontal resolution (try to keep it 1:1)
%__________________________________________________________________________
% 
% Hi0-tongH ________________________
%         o|                         -------___
%          |                       S2           \_
%          | S1                                    \
%         o|                                        \~~~~~~~~~~~~~~~~~~ 0.0
%          |                 S3                      \<-
% -tongH  o|_________________________________________|<--
%          o        o        o^   ^   ^   ^   ^  ^  ^ L
%                             |   |   |   |   |  |  |
%__________________________________________________________________________

% S1 (LHS)
ynodes = ceil(Hi0/Yres) + 1;
y1 = linspace(-tongH,Hi0,ynodes)';
x1 = zeros(ynodes,1);
f1 = ones(ynodes,1); 
f1(1) = 32; % this is the lower left corner of the BOTTOM side

% S3 (BOTTOM, count from RHS back to LHS)
% note : BOTTOM BCFLAC == 16 ; TOP BCFLAG == 32;
xnodes = ceil((L-Xres)/Xres) + 1;
x3 = linspace(L,Xres,xnodes)'; % this side won't include bottom corner
y3 = -tongH * ones(xnodes,1);
g_nds = find(x3 <= L*grounded);
slope3 = -dhill/(grounded*L);   % only use this stuff for bed topography
y3(g_nds) = (dhill - tongH) + slope3 * x3(g_nds);
f3 = 16*ones(xnodes,1);         % these nodes are 'floating'
f3(g_nds) = 4;                 % these nodes are 'grounded'

% S2 (QUARTER ELLIPSE)
[x2,y2,f2,elpsnodes] = makeQuarterEllipse(L,Hi0,tongH,x1,y1,x3,y3);

plot(x1,y1,'ko');hold on;
plot(x2,y2,'rx');
plot(x3,y3,'gs');

nnodes = xnodes + ynodes + elpsnodes;
nodcol = [0:nnodes-1]';
xs = [x1;x2;x3]; 
ys = [y1;y2;y3];
% getting perimeter into .poly file format
nodes  = [nodcol xs ys];
filler1 = [0; nodcol(1:end-1)];
filler2 = [nodcol(end); nodcol(2:end)];
BCflags = [f1;f2;f3];
facets = [nodcol filler1 filler2 BCflags];

% 'regional attributes' is required by TetGen
% and has the following format:
% [
%   #number materials
%   x-position of nth material
%   y-position of nth material
%   0   <--- what that means I don't know
%  -1   <--- some parameter I don't use
% ]
regional = [1 L/2 -0.5*(tongH) 0    -1];

% these commands print out the information from above to the MatLab
% screen and write that printout as a file (.poly) in the below 
% directory -- so make sure the MatLab window is big enough that the
% output doesn't get all janky looking.

format long g
cd \\Utig2\disk_student\logan\DynEarthSol\dv\DynEarthSol3D
fid = fopen('test.poly','w+');
diary 'test.poly'
[nnodes 2 0 0]
nodes
[nnodes 1]
facets
[0]
[1]
regional
diary off
fclose(fid);
