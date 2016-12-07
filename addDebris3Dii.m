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
str = 'lense3D_2.poly';

x_nodes = 100;                           % number nodes in X-dir
L = 1000;                                % length of domain in X-dir
x = linspace(0,L,x_nodes);               % X-dir


% simple testing:
% bas = zeros(length(x),1);
% srf = 100*ones(length(x),1);
% srf = srf(2:end-1);
% xsrf = x(2:end-1);

basal = 300*exp(-x/400);            % basal position function
a =  -28.3223;                      % ???
n_basal = 2.37e-03;                 % ???
c = 328.3223;                       % ???
surface = a*exp(n_basal*x)+c;       % surface position function
surface = surface(2:end-1);         % trimming the surface topo
x_surface = x(2:end-1);             % x corresponding to # surface pts.
%% dz determination

debris_start = 300;
debris_thickness = 10;
xmax = max(x);


xdebris = x(x > debris_start);
zdebris = [surface(x_surface > debris_start) basal(end)];

%find the slope of the debris thickness increase
m = debris_thickness/(xmax - debris_start);
dz = m*xdebris; %y = mx, equal to thickness of debris on top of surface

zdebris = zdebris + dz; 
debris_end = [debris_thickness+xmax, basal(end)]; %[x, z]

xdebris = [xdebris debris_end(1)+1.5*debris_thickness]; %Make it a bit extra thick at the base by adding 1.5* debris thickness to x value
zdebris = [zdebris debris_end(2)];

glacier_width = 100;                        % width between Y0 / Y1 faces

n_basal = length(basal);        %num basal nodes
n_surface = length(surface);    %num surface nodes
n_debris = length(xdebris);
%ib = 1:length(x);               % this is a single index of basal nodes
%is = 1:length(x_surface);       % this is a single index of surface nodes


%% Variables for important intersections

int1_b = n_basal
int2_b = 2*n_basal
int3_s = 2*n_basal + 1
int4_s = 2*n_basal + n_surface
int5_s = 2*n_basal + n_surface + 1
int6_s = 2*n_basal + 2*n_surface
int7_d = 2*n_basal + 2*n_surface + 1
int8_d = 2*n_basal + 2*n_surface + n_debris
int9_d = 2*n_basal + 2*n_surface + n_debris + 1
int10_d = 2*n_basal + 2*n_surface + 2*n_debris


%% ---------------------------- NODES ----------------------------------%%
nodes = NaN*ones(n_basal+n_basal+n_surface+n_surface+2*length(zdebris), 4);       % empty node matrix to be filled in

% BASAL nodes: ------------------------------------------------------------
% front side (step one, 1 to b):
for i = 1 : n_basal
    ii = i;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);             % x
    nodes(i,3) = glacier_width;     % y, front face
    nodes(i,4) = basal(ii);         % z, basal
end;
% back side: (note the loop counter) (step two, b+1 to 2b)
for i = (n_basal+1) : (2*n_basal)
    ii = i-n_basal;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);             % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = basal(ii);           % z, basal
end;

% SURFACE nodes: ----------------------------------------------------------
% front side: (step 4, 2b+1 to 2b+s)
for i = (2*n_basal) + 1 : (2*n_basal) + n_surface
    ii = i-2*n_basal;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);          % x
    nodes(i,3) = glacier_width;                 % y, back face
    nodes(i,4) = surface(ii);           % z, surface
end;
% back side (bounds of loop are written this way for clarity): (step 5,
% 2b+s+1 to 2b+2s)
for i = (2*n_basal + n_surface) + 1 : (2*n_basal + 2*n_surface)
    ii = i-(2*n_basal)-n_surface;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);          % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = surface(ii);           % z, surface
end;

% DEBRIS nodes: ----------------------------------------------------------
% front side:
for i = (2*n_basal + 2*n_surface) + 1 : (2*n_basal + 2*n_surface) + n_debris
    ii = i-2*n_basal-2*n_surface;               %incrementing through this section only
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xdebris(ii);     % x location of this node
    nodes(i,3) = glacier_width;     % y set to glacier width for all nodes
    nodes(i,4) = zdebris(ii);       % z, surface (I think this should be debris surface)
end;
% back side (bounds of loop are written this way for clarity):
for i = (2*n_basal + 2*n_surface + n_debris) + 1 : (2*n_basal + 2*n_surface) + 2*n_debris
    ii = i-(2*n_basal)-2*n_surface-n_debris; % incrementing through this section only
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xdebris(ii);     % x
    nodes(i,3) = 0;                 % y, back face
    nodes(i,4) = zdebris(ii);       % z, surface
end;


% DEBRIS nodes: -----------------------------------------------------------
% we don't need to add any nodes, just select a left / right position for
% the surface and basal sides

dx = mean(diff(x));
bl = 400; nbl = find((nodes(:,2) >= bl-dx) & (nodes(:,2) <= bl+dx)); %basal left?
          nbl = [nbl(1) nbl(3)]; % this is front / back side BOTTOM
br = 450; nbr = find((nodes(:,2) >= br-dx) & (nodes(:,2) <= br+dx)); %basal right?
          nbr = [nbr(1) nbr(3)]; % this is front / back side BOTTOM
sl = 425; nsl = find((nodes(:,2) >= sl-dx) & (nodes(:,2) <= sl+dx)); %surface left?
          nsl = [nsl(5) nsl(7)]; % this is front / back side SURFACE
sr = 475; nsr = find((nodes(:,2) >= sr-dx) & (nodes(:,2) <= sr+dx)); %surface right?
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
z0header = NaN*ones(n_basal-1,3); % 3 numbers make up a header
z0facets = NaN*ones(n_basal-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 16;                  % Z0 flag;
                

for f = 1 : n_basal-1
                             % Here's what the header means:
    z0header(f,1) = 1;       % "1 facet follows this header"
    z0header(f,2) = 0;       % "with 0 attributes"?
    z0header(f,3) = flag;    % "and this (Z0 = 16) flag"

    z0facets(f,1) = 4;      % # of points to follow (4, it's a quadrilateral)
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z0facets(f,2) = f+1;       %   b+f----b+f+1
    z0facets(f,3) = f;         %   |   Z0   |
    z0facets(f,4) = f+n_basal;       %   |        |
    z0facets(f,5) = f+1+n_basal;     %   f------f+1
end;

% Z1 (surface) ------------------------------------------------------------
% z1 face next: you will have (s-1) + 1 (left side) + 1 (right side) quadrilateral faces
z1header = NaN*ones(1 + n_surface-1 + 1 ,3); % 3 numbers make up a header
z1facets = NaN*ones(1 + n_surface-1 + 1 ,4); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 32;                  % Z1 flag;

% start w/ left side
z1header(1,1) = 1;
z1header(1,2) = 0;
z1header(1,3) = flag;

z1facets(1,1) = 4;              %  (b+1)     (2b+s+1)
z1facets(1,2) = 1;              %  ._________o
z1facets(1,3) = 2*n_basal + 1;        %  |         |
z1facets(1,4) = 2*n_basal + n_surface + 1;    %  |(1)      |
z1facets(1,5) = n_basal + 1;          %  ._________o (2b+1)

counter = 2; % start after the first facet listed above

% loop over surface
for f = (2*n_basal + 1) : (2*n_basal + n_surface) - 1
    
                             % Here's what the header means:
    z1header(counter,1) = 1;       % "1 facet follows this header"
    z1header(counter,2) = 0;       % "with 0 attributes"?
    z1header(counter,3) = flag;    % "and this (Z0 = 16) flag"


    z1facets(counter,1) = 4;      % # of points to follow 
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z1facets(counter,2) = f;        %   f+s----f+s+1
    z1facets(counter,3) = f+1;      %   |        |
    z1facets(counter,4) = f+n_surface+1;    %   |        |
    z1facets(counter,5) = f+n_surface;      %   f------f+1
    
    counter = counter + 1;
end;

% end w/ right side
z1header(end,1) = 1;
z1header(end,2) = 0;
z1header(end,3) = flag;

z1facets(end,1) = 4;              %  (2b+2s)   (2b)
z1facets(end,2) = n_basal;              %  o_________.
z1facets(end,3) = 2*n_basal;            %  |         |
z1facets(end,4) = 2*n_basal + 2*n_surface;      %  |(2b+s)   |
z1facets(end,5) = 2*n_basal + n_surface;        %  o_________. (b)

% Y0 (front side) ---------------------------------------------------------
% on the front side you'll have (s-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
y0header = NaN*ones(n_surface-1,3); % 3 numbers make up a header
y0facets = NaN*ones(n_surface-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 4;                   % Y0 flag;


% this is a really important array: 
%   the 1st column are front basal node numbers excluding the ends
%   the 2nd column are front surface node numbers:
indices = [[2:n_basal-1]' [2*n_basal+1:2*n_basal+n_surface]']; % call it I in the pic below

for f = 1 : n_surface-1
    y0header(f,1) = 1;
    y0header(f,2) = 0;
    y0header(f,3) = flag;

    y0facets(f,1) = 4;                   %  I(f,2)   I(f+1,2)
    y0facets(f,2) = indices(f,1);        %  o_________o
    y0facets(f,3) = indices(f+1,1);      %  |         |
    y0facets(f,4) = indices(f+1,2);      %  |I(f,1)   |
    y0facets(f,5) = indices(f,2);        %  ._________. I(f+1,1)
end

% Y1 (back side) ;) -------------------------------------------------------
% on the back side you'll have (s-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
y1header = NaN*ones(n_surface-1,3); % 3 numbers make up a header
y1facets = NaN*ones(n_surface-1,5); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 8;                   % Y1 flag;


% this is a really important array: 
%   the 1st column are front basal node numbers excluding the ends
%   the 2nd column are front surface node numbers: (should maybe say "back"
%   not "front"?
indices = [[n_basal+2:n_basal+n_basal-1]' [2*n_basal+n_surface+1:2*n_basal+2*n_surface]']; % call it I in the pic below

for f = 1 : n_surface-1
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
yTRIface(1,:) = [3      1       2       (2*n_basal+1)];    
% right:
yTRIhead(2,:) = [1 0 flag];
yTRIface(2,:) = [3      (n_basal-1)   n_basal       (2*n_basal+n_surface)];      

% Y1 (back) triangles 
flag = 8;

% left:
yTRIhead(3,:) = [1 0 flag];
yTRIface(3,:) = [3      (n_basal+2)   (n_basal+1)   (2*n_basal+n_surface+1)];    
% right:
yTRIhead(4,:) = [1 0 flag];
yTRIface(4,:) = [3      (2*n_basal)   (2*n_basal-1) (2*n_basal+2*n_surface)];

%% DEBRIS FACETS left / right side -----------------------------------------
debhead = NaN*ones(2,3);
debface = NaN*ones(2,5);
% there are only 2 facets for an inner debris band
% left:
debhead(1,:) = [1 0 0]; % there is no boundary flag for an internal facet
debface(1,:) = [4 nbl(1) nbl(2) nsl(2) nsl(1)];

% right:
debhead(2,:) = [1 0 0]; % there is no boundary flag for an internal facet
debface(2,:) = [4 nbr(1) nbr(2) nsr(2) nsr(1)];

% Y0D (front debris side) -------------------------------------------------
% on the front side you'll have (d-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
flag = 4;                   % Y0 flag;


% this is a really important array: 
%   the 1st column are front surface nodes excluding those before the
%   debris and at the ends
%   the 2nd column are the debris surface nodes excluding the start and the
%   end?
surf_index1 = 2*n_basal + 1 + (n_surface - n_debris) + 1; %+1 for indexing, +1 to skip the first
surf_index2 = 2*n_basal + n_surface-1; %+1 for indexing? -1 to skip the last
debris_index1 = (2*n_basal + 2*n_surface) + 1 + 1; %+1 for indexing, +1 to skip the first 
debris_index2 = 2*n_basal + 2*n_surface + n_debris - 1;
indices = [(surf_index1:surf_index2)' (debris_index1:debris_index2)']; % call it I in the pic below

%preallocate
y0header_debris = NaN*ones(length(indices),3);
y0facets_debris = NaN*ones(length(indices),5);

% Debris Y0
for f = 1 : length(indices)-1
    y0header_debris(f,1) = 1;
    y0header_debris(f,2) = 0;
    y0header_debris(f,3) = flag;

    y0facets_debris(f,1) = 4;                   %  I(f,2)   I(f+1,2)
    y0facets_debris(f,2) = indices(f,1);        %  o_________o
    y0facets_debris(f,3) = indices(f+1,1);      %  |         |
    y0facets_debris(f,4) = indices(f+1,2);      %  |I(f,1)   |
    y0facets_debris(f,5) = indices(f,2);        %  ._________. I(f+1,1)
end


% Y1D (back debris side) -------------------------------------------------
% on the front side you'll have (d-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
flag = 8;                   % Y1 flag;


% this is a really important array: 
%   the 1st column are front surface nodes excluding those before the
%   debris and at the ends
%   the 2nd column are the debris surface nodes excluding the start and the
%   end?
surf_index1 = 2*n_basal + n_surface + 1 + (n_surface - n_debris) + 1; %+1 for indexing, +1 to skip the first
surf_index2 = 2*n_basal + 2*n_surface -1; %+1 for indexing
debris_index1 = (2*n_basal + 2*n_surface) + n_debris + 1 + 1; %+1 for indexing, +1 to skip the first 
debris_index2 = 2*n_basal + 2*n_surface + 2*n_debris - 1 ; % -1 to skip the last
indices = [(surf_index1:surf_index2)' (debris_index1:debris_index2)']; % call it I in the pic below

%preallocate
y1header_debris = zeros(length(indices),3);
y1facets_debris = zeros(length(indices),5);

% Debris Y0
for f = 1 : length(indices)-1
    y1header_debris(f,1) = 1;
    y1header_debris(f,2) = 0;
    y1header_debris(f,3) = flag;

    y1facets_debris(f,1) = 4;                   %  I(f,2)   I(f+1,2)
    y1facets_debris(f,2) = indices(f,1);        %  o_________o
    y1facets_debris(f,3) = indices(f+1,1);      %  |         |
    y1facets_debris(f,4) = indices(f+1,2);      %  |I(f,1)   |
    y1facets_debris(f,5) = indices(f,2);        %  ._________. I(f+1,1)
end

%% Figure out debris layer triangular facets

% Y0D / Y1D TRIANGULAR FACETS ON LEFT AND RIGHT SIDE ----------------------
% there are 4 of these

yTRIhead_debris = NaN*ones(4,3);
yTRIface_debris = NaN*ones(4,4); % now we have 3 pts defining a facet, so 4 cols

% Y0 (front) triangles
flag = 4;

% left:
yTRIhead_debris(1,:) = [1 0 flag];
yTRIface_debris(1,:) = [3      2*n_basal+1+(n_surface-n_debris)     2*n_basal+2+(n_surface - n_debris)   (2*n_basal+2*n_surface+1)]; %not sure what the first numbers are here
% right:
yTRIhead_debris(2,:) = [1 0 flag];
yTRIface_debris(2,:) = [3      (2*n_basal+n_surface)   (2*n_basal+2*n_surface+n_debris-1)  (2*n_basal+2*n_surface+n_debris)]; %what does the 3 mean

% Y1 (back) triangles 
flag = 8;

% left:
yTRIhead_debris(3,:) = [1 0 flag];
yTRIface_debris(3,:) = [3      (2*n_basal+2*n_surface-n_debris+2)   (2*n_basal+2*n_surface-n_debris+1)   (2*n_basal+2*n_surface+n_debris+1)];    
% right:
yTRIhead_debris(4,:) = [1 0 flag];
yTRIface_debris(4,:) = [3      (2*n_basal+2*n_surface)   (2*n_basal+2*n_surface-1) (2*n_basal+2*n_surface+2*n_debris)];

%% BECAUSE C++ / TETGEN LIBRARY NUMBER EVERYTHING STARTING FROM 0, WE 
%  SUBTRACT 1 EVERYWHERE.

% concatenate all the header / facet lists for QUADRILATERAL ELEMENTS:
header4 = [z0header;
          z1header;
          y0header;
          y1header;
          debhead;
          y0header_debris;
          y1header_debris];

facets4 = [z0facets;
          z1facets;
          y0facets;
          y1facets
          debface;
          y0facets_debris;
          y1facets_debris];

tri_facets = [yTRIface;
              yTRIface_debris];
          
tri_headers = [yTRIhead;
               yTRIhead_debris];

      
% now subtract 1 from the node list:
nodes(:,1) = nodes(:,1) - 1;
% now subtract 1 from the facets4 and yTRIface and yTRIface_debris lists:
facets4 (:,2:5) = facets4 (:,2:5) - 1;
tri_facets(:,2:4) = tri_facets(:,2:4) - 1;


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
n_facets = length(header4) + length(tri_headers); %why are we counting facets by headers
header = [n_facets 1];
fmt    = '%d %d\n';
fprintf(fid, fmt, header);

% write quadrilateral facets 
% (alternate header and facet row by row)
n_quads = length(header4);
fmthead  = '%d %d %d\n';
fmtfacet = '%d %d %d %d %d\n';
for n = 1:n_quads
    fprintf(fid, fmthead , header4(n,:));
    fprintf(fid, fmtfacet, facets4(n,:));
end;

% same for triangular elements:
n_tris = length(tri_headers);
fmtfacet = '%d %d %d %d\n';
for n = 1:n_tris
    fprintf(fid, fmthead , tri_headers(n,:));
    fprintf(fid, fmtfacet, tri_facets(n,:));
end;

% our mesh has 0 holes:
fprintf(fid, '%d\n', 0);

% number of regions?
REGIONS = 3; % ICE --- ROCK --- ICE
fprintf(fid, '%d\n', REGIONS);

% region description:
%         [region #     x0      y0          z0                      nmat    maximum elem. size]
regions = [0            bl-dx   glacier_width/2     nodes(nbl(1)-3,4)+dx     0           1e9;
           1            bl+dx   glacier_width/2     nodes(nbr(1)-2,4)+dx     1           1e9;
           2            br+dx   glacier_width/2     nodes(nbr(1)+5,4)+dx     0           1e9];

fmt = '%d %f %f %f %d %f\n';
for n = 1:REGIONS
    fprintf(fid,fmt,regions(n,:));
end;
fclose(fid);
      












































