function [S1,S2,S4] = ellipseChunk(L,grounded,Hi0,dhill,tongH,Yres,Xres);

% S1 (LHS)
ynodes = ceil(Hi0/Yres) + 1;
%ynodes = 10;
y1 = linspace(-tongH + dhill,Hi0,ynodes)';
x1 = zeros(ynodes,1);
f1 = ones(ynodes,1); 
f1(1) = 16; % this is the lower left corner of the BOTTOM side

% S4 (BOTTOM, count from RHS back to LHS)
% note : BOTTOM BCFLAC == 16 ; TOP BCFLAG == 32;
xnodes = ceil((L-Xres)/Xres) + 1;
%xnodes = 50;
x4 = linspace(L,Xres,xnodes)'; % this side won't include bottom corner
y4 = -tongH * ones(xnodes,1);
g_nds = find(x4 <= L*grounded);
slope4 = -dhill/(grounded*L);   % only use this stuff for bed topography
y4(g_nds) = (dhill - tongH) + slope4 * x4(g_nds);
f4 = 32*ones(xnodes,1);         % these nodes are 'floating'
f4(g_nds) = 16;                 % these nodes are 'grounded'

% S2 (QUARTER ELLIPSE)
[x2,y2,f2,elpsnodes] = makeQuarterEllipse(L,Hi0,tongH,x1,y1,x4,y4);

chunkNodes = xnodes + ynodes + elpsnodes;
nodcol = [0:chunkNodes-1]';
S1 = [x1 y1 f1]; 
S2 = [x2 y2 f2];
S4 = [x4 y4 f4];
