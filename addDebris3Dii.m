%% Here's how the node list is going to be set up:
%  . == basal nodes
%  o == surface nodes
% (xxx) == node # in node list 
%
%       o(2*b+s+1)____
%      /              ____
% (2*b+1)  .(b+1)         ---_____             
%    o____                         --------____   (2*b + 2*s)        
%  ./         ------_________                ----____o
%  1\____                    ___________________    / \. (2*b)
%        \__________                            ---o  /(2*b + s)
%                   -----------____________________\./ (b)
%
%%
clear;clf;
% first: what are you going to call this poly file?
str = 'lense23D.poly';

nd = 100;                           % number nodes in X-dir
L = 1000;                           % length of domain in X-dir
x = linspace(0,L,nd);               % X-dir


% simple testing:
% bas = zeros(length(x),1);
% srf = 100*ones(length(x),1);
% srf = srf(2:end-1);
% xsrf = x(2:end-1);

bas = 300*exp(-x/400);              % basal topo
a =  -28.3223;
b = 2.37e-03;
c = 328.3223;
srf = a*exp(b*x)+c;                 % surface topo
srf = srf(2:end-1);                 % trimming the surface topo
xsrf = x(2:end-1);                  % x corresponding to # surface pts.
%% dy determination

debrisStart = 300;
debrisThickness = 10;
xmax = max(x);


xdebris = x(x > debrisStart);
ydebris = [srf(xsrf > debrisStart) bas(end)];

%find the slope of the debris thickness increase
m = debrisThickness/(xmax - debrisStart);
dy = m*xdebris; %y = mx

dtopo = ydebris + dy; 
debris_end = [debrisThickness+xmax, bas(end)]; %[x, y]

xdebris = [xdebris debris_end(1)];
ydebris = [dtopo debris_end(2)];

width = 100;                        % width between Y0 / Y1 faces

b = length(bas); s = length(srf);   % number basal / surface nodes
ib = [1:length(x)];                 % this is a single index of basal nodes
is = [1:length(xsrf)];              % this is a single index of surface nodes


%% ---------------------------- NODES ----------------------------------%%
nodes = NaN*ones(b+b+s+s+2*length(dtopo), 4);       % empty node matrix to be filled in

% BASAL nodes: ------------------------------------------------------------
% front side:
for i = 1 : b
    ii = i;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);             % x
    nodes(i,3) = width;             % y, front face
    nodes(i,4) = bas(ii);           % z, basal
end;
% back side: (note the loop counter)
for i = (b+1) : (2*b)
    ii = i-b;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);             % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = bas(ii);           % z, basal
end;

% SURFACE nodes: ----------------------------------------------------------
% front side:
for i = (2*b) + 1 : (2*b) + s
    ii = i-2*b;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xsrf(ii);          % x
    nodes(i,3) = width;                 % y, back face
    nodes(i,4) = srf(ii);           % z, surface
end;
% back side (bounds of loop are written this way for clarity):
for i = (2*b + s) + 1 : (2*b + s) + s
    ii = i-(2*b)-s;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xsrf(ii);          % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = srf(ii);           % z, surface
end;

% DEBRIS nodes: ----------------------------------------------------------
% front side:
for i = (2*b) + 1 : (2*b) + s
    ii = i-2*b;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xsrf(ii);          % x
    nodes(i,3) = width;                 % y, back face
    nodes(i,4) = srf(ii);           % z, surface
end;
% back side (bounds of loop are written this way for clarity):
for i = (2*b + 2*s) + 1 : (2*b + 2*s) + 
    ii = i-(2*b)-s;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xsrf(ii);          % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = srf(ii);           % z, surface
end;


% DEBRIS nodes: -----------------------------------------------------------
% we don't need to add any nodes, just select a left / right position for
% the surface and basal sides

dx = mean(diff(x));
bl = 400; nbl = find((nodes(:,2) >= bl-dx) & (nodes(:,2) <= bl+dx)); 
          nbl = [nbl(1) nbl(3)]; % this is front / back side BOTTOM
br = 450; nbr = find((nodes(:,2) >= br-dx) & (nodes(:,2) <= br+dx)); 
          nbr = [nbr(1) nbr(3)]; % this is front / back side BOTTOM
sl = 425; nsl = find((nodes(:,2) >= sl-dx) & (nodes(:,2) <= sl+dx)); 
          nsl = [nsl(5) nsl(7)]; % this is front / back side SURFACE
sr = 475; nsr = find((nodes(:,2) >= sr-dx) & (nodes(:,2) <= sr+dx)); 
          nsr = [nsr(5) nsr(7)]; % this is front / back side SURFACE
% PLOT for sanity check:
plot3(nodes(:,2)  ,nodes(:,3)  ,nodes(:,4),'k.'); hold on;
plot3(nodes(nbl,2),nodes(nbl,3),nodes(nbl,4),'ro')
plot3(nodes(nsl,2),nodes(nsl,3),nodes(nsl,4),'ro')
plot3(nodes(nbr,2),nodes(nbr,3),nodes(nbr,4),'bs')
plot3(nodes(nsr,2),nodes(nsr,3),nodes(nsr,4),'bs')

%% ---------------------------- FACETS ---------------------------------%%

% Z0 (basal) --------------------------------------------------------------
% z0 face first: you will have (b-1) quadrilateral faces
z0header = NaN*ones(b-1,3); % 3 numbers make up a header
z0facets = NaN*ones(b-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 16;                  % Z0 flag;
                

for f = 1 : b-1
                             % Here's what the header means:
    z0header(f,1) = 1;       % "1 facet follows this header"
    z0header(f,2) = 0;       % "with 0 attributes"?
    z0header(f,3) = flag;    % "and this (Z0 = 16) flag"

    z0facets(f,1) = 4;      % # of points to follow (4, it's a quadrilateral)
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z0facets(f,2) = f+1;       %   b+f----b+f+1
    z0facets(f,3) = f;         %   |   Z0   |
    z0facets(f,4) = f+b;       %   |        |
    z0facets(f,5) = f+1+b;     %   f------f+1
end;

% Z1 (surface) ------------------------------------------------------------
% z1 face next: you will have (s-1) + 1 (left side) + 1 (right side) quadrilateral faces
z1header = NaN*ones(1 + s-1 + 1 ,3); % 3 numbers make up a header
z1facets = NaN*ones(1 + s-1 + 1 ,4); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 32;                  % Z1 flag;

% start w/ left side
z1header(1,1) = 1;
z1header(1,2) = 0;
z1header(1,3) = flag;

z1facets(1,1) = 4;              %  (b+1)     (2b+s+1)
z1facets(1,2) = 1;              %  ._________o
z1facets(1,3) = 2*b + 1;        %  |         |
z1facets(1,4) = 2*b + s + 1;    %  |(1)      |
z1facets(1,5) = b + 1;          %  ._________o (2b+1)

counter = 2; % start after the first facet listed above

% loop over surface
for f = (2*b + 1) : (2*b + s) - 1
    
                             % Here's what the header means:
    z1header(counter,1) = 1;       % "1 facet follows this header"
    z1header(counter,2) = 0;       % "with 0 attributes"?
    z1header(counter,3) = flag;    % "and this (Z0 = 16) flag"


    z1facets(counter,1) = 4;      % # of points to follow 
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z1facets(counter,2) = f;        %   f+s----f+s+1
    z1facets(counter,3) = f+1;      %   |        |
    z1facets(counter,4) = f+s+1;    %   |        |
    z1facets(counter,5) = f+s;      %   f------f+1
    
    counter = counter + 1;
end;

% end w/ right side
z1header(end,1) = 1;
z1header(end,2) = 0;
z1header(end,3) = flag;

z1facets(end,1) = 4;              %  (2b+2s)   (2b)
z1facets(end,2) = b;              %  o_________.
z1facets(end,3) = 2*b;            %  |         |
z1facets(end,4) = 2*b + 2*s;      %  |(2b+s)   |
z1facets(end,5) = 2*b + s;        %  o_________. (b)

% Y0 (front side) ---------------------------------------------------------
% on the front side you'll have (s-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
y0header = NaN*ones(s-1,3); % 3 numbers make up a header
y0facets = NaN*ones(s-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 4;                   % Y0 flag;


% this is a really important array: 
%   the 1st column are front basal node numbers excluding the ends
%   the 2nd column are front surface node numbers:
indices = [[2:b-1]' [2*b+1:2*b+s]']; % call it I in the pic below

for f = 1 : s-1
    y0header(f,1) = 1;
    y0header(f,2) = 0;
    y0header(f,3) = flag;

    y0facets(f,1) = 4;                   %  I(f,2)   I(f+1,2)
    y0facets(f,2) = indices(f,1);        %  o_________o
    y0facets(f,3) = indices(f+1,1);      %  |         |
    y0facets(f,4) = indices(f+1,2);      %  |I(f,1)   |
    y0facets(f,5) = indices(f,2);        %  ._________. I(f+1,1)
end

% Y1 (back side) ----------------------------------------------------------
% on the back side you'll have (s-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
y1header = NaN*ones(s-1,3); % 3 numbers make up a header
y1facets = NaN*ones(s-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 8;                   % Y0 flag;


% this is a really important array: 
%   the 1st column are front basal node numbers excluding the ends
%   the 2nd column are front surface node numbers:
indices = [[b+2:b+b-1]' [2*b+s+1:2*b+2*s]']; % call it I in the pic below

for f = 1 : s-1
    y1header(f,1) = 1;
    y1header(f,2) = 0;
    y1header(f,3) = flag;               % Note that outward normal is INTO
                                        % page for this order:
    y1facets(f,1) = 4;                  %  I(f,2)   I(f+1,2)
    y1facets(f,2) = indices(f+1,1);     %  o_________o
    y1facets(f,3) = indices(f,1);       %  |         |
    y1facets(f,4) = indices(f,2);       %  |I(f,1)   |
    y1facets(f,5) = indices(f+1,2);     %  ._________. I(f+1,1)
end

% Y0 / Y1 TRIANGULAR FACETS ON LEFT AND RIGHT SIDE ------------------------
% there are 4 of these

yTRIhead = NaN*ones(4,3);
yTRIface = NaN*ones(4,4); % now we have 3 pts defining a facet, so 4 cols

% Y0 (front) triangles
flag = 4;

% left:
yTRIhead(1,:) = [1 0 flag];
yTRIface(1,:) = [3      1       2       (2*b+1)];    
% right:
yTRIhead(2,:) = [1 0 flag];
yTRIface(2,:) = [3      (b-1)   b       (2*b+s)];      

% Y1 (back) triangles 
flag = 8;

% left:
yTRIhead(3,:) = [1 0 flag];
yTRIface(3,:) = [3      (b+2)   (b+1)   (2*b+s+1)];    
% right:
yTRIhead(4,:) = [1 0 flag];
yTRIface(4,:) = [3      (2*b)   (2*b-1) (2*b+2*s)];

% DEBRIS FACETS left / right side -----------------------------------------
debhead = NaN*ones(2,3);
debface = NaN*ones(2,5);
% there are only 2 facets for an inner debris band
% left:
debhead(1,:) = [1 0 0]; % there is no boundary flag for an internal facet
debface(1,:) = [4 nbl(1) nbl(2) nsl(2) nsl(1)];

% right:
debhead(2,:) = [1 0 0]; % there is no boundary flag for an internal facet
debface(2,:) = [4 nbr(1) nbr(2) nsr(2) nsr(1)];

%% BECAUSE C++ / TETGEN LIBRARY NUMBER EVERYTHING STARTING FROM 0, WE 
%  SUBTRACT 1 EVERYWHERE.

% concatenate all the header / facet lists for QUADRILATERAL ELEMENTS:
header4 = [z0header;
          z1header;
          y0header;
          y1header;
          debhead];

facets4 = [z0facets;
          z1facets;
          y0facets;
          y1facets
          debface];

% now subtract 1 from the node list:
nodes(:,1) = nodes(:,1) - 1;
% now subtract 1 from the facets4 and yTRIface lists:
facets4 (:,2:5) = facets4 (:,2:5) - 1;
yTRIface(:,2:4) = yTRIface(:,2:4) - 1;

%% NOW WE WRITE TO FILE
format short g
fid = fopen(str,'w');

% node header:
fmt    = '%d %d %d %d\n';
header = [length(nodes) 3   0   0];
fprintf(fid, fmt, header);

% write nodes:
fmt    = '%d %f %f %f\n';
for n = 1:length(nodes)
    fprintf(fid, fmt, nodes(n,:));
end;

% facet header:
nofacets = [length(header4) + length(yTRIhead)];
header = [nofacets 1];
fmt    = '%d %d\n';
fprintf(fid, fmt, header);

% write quadrilateral facets 
% (alternate header and facet row by row)
quads = length(header4);
fmthead  = '%d %d %d\n';
fmtfacet = '%d %d %d %d %d\n';
for n = 1:quads
    fprintf(fid, fmthead , header4(n,:));
    fprintf(fid, fmtfacet, facets4(n,:));
end;

% same for triangular elements:
tris = length(yTRIhead);
fmtfacet = '%d %d %d %d\n';
for n = 1:tris
    fprintf(fid, fmthead , yTRIhead(n,:));
    fprintf(fid, fmtfacet, yTRIface(n,:));
end;

% our mesh has 0 holes:
fprintf(fid, '%d\n', 0);

% number of regions?
REGIONS = 3; % ICE --- ROCK --- ICE
fprintf(fid, '%d\n', REGIONS);

% region description:
%         [region #     x0      y0          z0                      nmat    maximum elem. size]
regions = [0            bl-dx   width/2     nodes(nbl(1)-3,4)+dx     0           1e9;
           1            bl+dx   width/2     nodes(nbr(1)-2,4)+dx     1           1e9;
           2            br+dx   width/2     nodes(nbr(1)+5,4)+dx     0           1e9];

fmt = '%d %f %f %f %d %f\n';
for n = 1:REGIONS
    fprintf(fid,fmt,regions(n,:));
end;
fclose(fid);
      











































