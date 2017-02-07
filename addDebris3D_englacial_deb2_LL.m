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
str = 'lense3D_endeb3.poly';

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
debris_thickness = 20;
englacial_debris_start = 350;
englacial_debris_thickness = 10;
debris_thickness = 10;
xmax = max(x);


xdebris = x(x > debris_start);
zdebris = [surface(x_surface > debris_start) basal(end)];

%find the slope of the debris thickness increase
m = debris_thickness/(xmax - debris_start);
dz = m*(xdebris-debris_start+150); %y = mx, equal to thickness of debris on top of surface

zdebris = zdebris + dz; 
debris_end = [debris_thickness+xmax, basal(end)]; %[x, z]

xdebris = [xdebris debris_end(1)+1.1*debris_thickness]; %Make it a bit extra thick at the base by adding 1.5* debris thickness to x value
zdebris = [zdebris debris_end(2)];

glacier_width = 100;                        % width between Y0 / Y1 faces

n_basal = length(basal);        %num basal nodes
n_surface = length(surface);    %num surface nodes
n_debris = length(xdebris);

% useful corner nodes:

b1 = 1;
b2 = n_basal;
b3 = n_basal+1;
b4 = 2*n_basal;

s1 = 2*n_basal + 1;
s2 = 2*n_basal + n_surface;
s3 = 2*n_basal + n_surface + 1;
s4 = 2*n_basal + 2*n_surface;

d1 = 2*n_basal + 2*n_surface + 1;
d2 = 2*n_basal + 2*n_surface + n_debris;
d3 = 2*n_basal + 2*n_surface + n_debris + 1;
d4 = 2*n_basal + 2*n_surface + 2*n_debris;

ds1 = 2*n_basal + n_surface - n_debris + 3;
ds2 = 2*n_basal + 2*n_surface - n_debris + 3;

corners = [b1 b2 b3 b4 s1 s2 s3 s4 d1 d2 d3 d4 ds1 ds2];

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
    nodes(i,2) = x_surface(ii);          % x
    nodes(i,3) = glacier_width;                 % y, back face
    nodes(i,4) = surface(ii);           % z, surface
end;
% back side (bounds of loop are written this way for clarity): (step 5,
% 2b+s+1 to 2b+2s)
for i = (2*n_basal + n_surface) + 1 : (2*n_basal + 2*n_surface)
    ii = i-(2*n_basal)-n_surface;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x_surface(ii);          % x
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

% PLOT for sanity check:
plot3(nodes(:,2)  ,nodes(:,3)  ,nodes(:,4),'k.', 'MarkerSize', 10); hold on;

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
    z0facets(f,2) = f+1;             %   b+f----b+f+1
    z0facets(f,3) = f;               %   |   Z0   |
    z0facets(f,4) = f+n_basal;       %   |        |
    z0facets(f,5) = f+1+n_basal;     %   f------f+1
end;

% Z1 ice (surface) ------------------------------------------------------------
% z1 face next: you will have (s-1) + 1 (left side) + 1 (right side) quadrilateral faces
z1header = NaN*ones(1 + n_surface-1 + 1 ,3); % 3 numbers make up a header
z1facets = NaN*ones(1 + n_surface-1 + 1 ,4); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 32;                  % Z1 flag;

% start w/ left side
z1header(1,1) = 1;
z1header(1,2) = 0;
z1header(1,3) = flag;

z1facets(1,1) = 4;                          %  (b+1)     (2b+s+1)
z1facets(1,2) = 1;                            %  ._________o
z1facets(1,3) = 2*n_basal + 1;                %  |         |
z1facets(1,4) = 2*n_basal + n_surface + 1;    %  |(1)      |
z1facets(1,5) = n_basal + 1;                  %  ._________o (2b+1)

counter = 2; % start after the first facet listed above

% loop over surface
for f = (2*n_basal + 1) : (2*n_basal + n_surface) - 1
    
                             % Here's what the header means:
    z1header(counter,1) = 1;       % "1 facet follows this header"
    z1header(counter,2) = 0;       % "with 0 attributes"?
    if (f >= ds1)
        flag = 0;
    end
    z1header(counter,3) = flag;    % "and this (Z0 = 16) flag"


    z1facets(counter,1) = 4;      % # of points to follow 
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z1facets(counter,2) = f;                %   f+s----f+s+1
    z1facets(counter,3) = f+1;              %   |        |
    z1facets(counter,4) = f+n_surface+1;    %   |        |
    z1facets(counter,5) = f+n_surface;      %   f------f+1
    
    counter = counter + 1;
end;

% end w/ right side
z1header(end,1) = 1;
z1header(end,2) = 0;
z1header(end,3) = 0;

z1facets(end,1) = 4;                            %  (2b+2s)   (2b)
z1facets(end,2) = n_basal;                      %  o_________.
z1facets(end,3) = 2*n_basal;                    %  |         |
z1facets(end,4) = 2*n_basal + 2*n_surface;      %  |(2b+s)   |
z1facets(end,5) = 2*n_basal + n_surface;        %  o_________. (b)

% Z1 debris (surface) ------------------------------------------------------------
% z1 face next: you will have (s-1) + 1 (left side) + 1 (right side) quadrilateral faces
z1debris_header = NaN*ones(1 + n_debris  ,3); % 3 numbers make up a header
z1debris_facets = NaN*ones(1 + n_debris  ,4); % 1 number saying the number of nodes we'll have
                            % and 4 nodes comprising the actual facet
flag = 32;                  % Z1 flag;

% start w/ left side
z1debris_header(1,1) = 1;
z1debris_header(1,2) = 0;
z1debris_header(1,3) = flag;

z1debris_facets(1,1) = 4;                          
z1debris_facets(1,2) = ds1-1;                         
z1debris_facets(1,3) = d1;                
z1debris_facets(1,4) = d3;   
z1debris_facets(1,5) = ds2-1;                 

counter = 2; % start after the first facet listed above

% loop over surface
for f = d1 : d2 - 1
    
                             % Here's what the header means:
    z1debris_header(counter,1) = 1;       % "1 facet follows this header"
    z1debris_header(counter,2) = 0;       % "with 0 attributes"?
    z1debris_header(counter,3) = flag;    % "and this (Z0 = 32) flag"


    z1debris_facets(counter,1) = 4;      % # of points to follow 
                            % "listed in this order" (follows right-hand
                            % rule, normal vector OUTWARD):
    z1debris_facets(counter,2) = f;                %   f+s----f+s+1
    z1debris_facets(counter,3) = f+1;              %   |        |
    z1debris_facets(counter,4) = f+n_debris+1;     %   |        |
    z1debris_facets(counter,5) = f+n_debris;       %   f------f+1
    
    counter = counter + 1;
end;

% last facet tucks under, connecting last debris node to last basal ice
% mode
flag = 16;
z1debris_header(end,1) = 1;
z1debris_header(end,2) = 0;
z1debris_header(end,3) = flag;

z1debris_facets(end,1) = 4;                          
z1debris_facets(end,2) = d4;                         
z1debris_facets(end,3) = d2;                
z1debris_facets(end,4) = b2;   
z1debris_facets(end,5) = b4;

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

% %% DEBRIS FACETS interior        ------------------------------------------
% 
deb_englacial_header = NaN*ones(2,3); 

deb_englacial_facets = NaN*ones(2,5);

%find the intersecting points for the interior debris layer
%I mixed up Y and Z again. should probably fix and propagate through
surf_left_x = x_surface(x_surface > englacial_debris_start);
surf_left_y = surface(x_surface > englacial_debris_start);
surf_left_x = surf_left_x(1);
surf_left_y = surf_left_y(1);
n_surf_lefta = getNode(nodes,surf_left_x,0,surf_left_y);
n_surf_leftb = getNode(nodes,surf_left_x,100,surf_left_y);


surf_right_x = x_surface(x_surface > englacial_debris_start + englacial_debris_thickness);
surf_right_y = surface(x_surface > englacial_debris_start + englacial_debris_thickness);
surf_right_x = surf_right_x(1);
surf_right_y = surf_right_y(1);
n_surf_righta = getNode(nodes,surf_right_x,0,surf_right_y);
n_surf_rightb = getNode(nodes,surf_right_x,100,surf_right_y);


basal_left_x = x(x>debris_start);
basal_left_y = basal(x>debris_start);
basal_left_x = basal_left_x(1);
basal_left_y = basal_left_y(1);
n_basal_lefta = getNode(nodes,basal_left_x,0,basal_left_y);
n_basal_leftb = getNode(nodes,basal_left_x,100,basal_left_y);


basal_right_x = x(x>debris_start+englacial_debris_thickness);
basal_right_y = basal(x>debris_start+englacial_debris_thickness);
basal_right_x = basal_right_x(1);
basal_right_y = basal_right_y(1);
n_basal_righta = getNode(nodes,basal_right_x,0,basal_right_y);
n_basal_rightb = getNode(nodes,basal_right_x,100,basal_right_y);

% there are only 2 facets for an inner debris band going all the way
% through the glacier (top-to-bottom)

% left:
nd = 45;

deb_englacial_header(1,:) = [1 0 0]; % there is no boundary flag for an internal facet

%deb_englacial_facets(1,:) = [4 n_basal_lefta n_basal_leftb n_surf_leftb n_surf_lefta];
deb_englacial_facets(1,:) = [4 y0facets(nd,2) y0facets(nd,5) y1facets(nd,4) y1facets(nd,3)];


% right:

deb_englacial_header(2,:) = [1 0 0]; % there is no boundary flag for an internal facet

%deb_englacial_facets(2,:) = [4 n_basal_righta n_basal_rightb  n_surf_rightb n_surf_righta];
deb_englacial_facets(2,:) = [4 y0facets(nd,3) y0facets(nd,4) y1facets(nd,5) y1facets(nd,2)];


%Add to plot of debris nodes
x_intersects = [basal_left_x,basal_right_x,surf_left_x,surf_right_x];
y_intersects = [basal_left_y,basal_right_y,surf_left_y,surf_right_y];
plot3(x_intersects, [0 0 0 0]  ,y_intersects,'ro'); hold on;
plot3(x_intersects, [100 100 100 100], y_intersects,'ro'); hold on;

%% DEBRIS FACETS left / right side ----------------------------------------

% Y0D (front debris side) -------------------------------------------------
% on the front side you'll have (d-1) quadrilateral elements and (2)
% triangular elements (on the left and right sides -- deal w/ those last)
flag = 4;                   % Y0 flag;

% this is a really important array: 
%   the 1st column are front surface nodes excluding those before the
%   debris and at the ends
%   the 2nd column are the debris surface nodes excluding the start and the
%   end?
surf_indices = [(ds1:s2) b2]; %+1 for indexing, +1 to skip the first
debris_indices = [d1:d2-1]; %+1 for indexing, +1 to skip the first 

indices = [surf_indices' debris_indices']; % call it I in the pic below

%preallocate
y0header_debris = NaN*ones(length(indices)-1,3);
y0facets_debris = NaN*ones(length(indices)-1,5);

% Debris Y0
for f = 1 : length(indices)-1
    y0header_debris(f,1) = 1;
    y0header_debris(f,2) = 0;
    y0header_debris(f,3) = flag;

    y0facets_debris(f,1) = 4;                   %  I(f,2)   I(f+1,2)
    y0facets_debris(f,2) = indices(f,1);        %  o_________o
    y0facets_debris(f,3) = indices(f,2);      %  |         |
    y0facets_debris(f,4) = indices(f+1,2);      %  |I(f,1)   |
    y0facets_debris(f,5) = indices(f+1,1);        %  ._________. I(f+1,1)
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

surf_indices = [(ds2:s4) b4];
debris_indices = [d3:d4-1];

indices = [surf_indices' debris_indices']; % call it I in the pic below

%preallocate
y1header_debris = zeros(length(indices)-1,3);
y1facets_debris = zeros(length(indices)-1,5);

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
yTRIface_debris(1,:) = [3 d1 ds1 ds1-1]; %not sure what the first numbers are here
% right:
yTRIhead_debris(2,:) = [1 0 flag];
yTRIface_debris(2,:) = [3 b2 d2-1 d2]; %what does the 3 mean

% Y1 (back) triangles 
flag = 8;

% left:
yTRIhead_debris(3,:) = [1 0 flag];
yTRIface_debris(3,:) = [3 d3 ds2-1 ds2];    
% right:
yTRIhead_debris(4,:) = [1 0 flag];
yTRIface_debris(4,:) = [3 d4 d4-1 b4];

%% BECAUSE C++ / TETGEN LIBRARY NUMBER EVERYTHING STARTING FROM 0, WE 
%  SUBTRACT 1 EVERYWHERE.

% concatenate all the header / facet lists for QUADRILATERAL ELEMENTS:
header4 = [z0header;
          z1header;
          z1debris_header;
          y0header;
          y1header;
          y0header_debris;
          y1header_debris;
          deb_englacial_header];
      
facets4 = [z0facets;
          z1facets;
          z1debris_facets;
          y0facets;
          y1facets;
          y0facets_debris;
          y1facets_debris;
          deb_englacial_facets]; %

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
REGIONS = 2; % ICE --- ROCK --- ICE
fprintf(fid, '%d\n', REGIONS);

% region description:
%         [region #     x0        y0                z0      nmat    maximum elem. size]
regions = [0            mean(x)   glacier_width/2   100     0           1e9;
           1            mean(x)   glacier_width/2   240     1           1e9];

fmt = '%d %f %f %f %d %f\n';
for n = 1:REGIONS
    fprintf(fid,fmt,regions(n,:));
end;
fclose(fid);
