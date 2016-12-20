% Ellipse 2D
%Making an idealized 2D model for DynEarthSol

clear;clf;
% first: what are you going to call this poly file?
str = 'ellipse2D.poly';

x_nodes = 100;                           % number nodes in X-dir
L = 1000;                                % length of domain in X-dir
x = linspace(0,L,x_nodes);               % X-dir


basal = 300*exp(-x/400);            % basal position function
a =  -28.3223;                      % ???
n_basal = 2.37e-03;                 % ???
c = 328.3223;                       % ???
surface = a*exp(n_basal*x)+c;       % surface position function
surface = surface(2:end-1);         % trimming the surface topo
x_surface = x(2:end-1);             % x corresponding to # surface pts.

%% debris geometries

debris_start = 300;
debris_thickness = 20;
englacial_debris_start = 350;
englacial_debris_thickness = 20;
debris_thickness = 10;
xmax = max(x);

xdebris = x(x > debris_start);
zdebris = [surface(x_surface > debris_start) basal(end)];

%find the slope of the debris thickness increase
m = debris_thickness/(xmax - debris_start);
dz = m*(xdebris-debris_start+150); %y = mx, equal to thickness of debris on top of surface

zdebris = zdebris + dz; 
debris_end = [debris_thickness+xmax, basal(end)]; %[x, z]

xdebris = [xdebris debris_end(1)+1.5*debris_thickness]; %Make it a bit extra thick at the base by adding 1.5* debris thickness to x value
zdebris = [zdebris debris_end(2)];

n_basal = length(basal);        %num basal nodes
n_surface = length(surface);    %num surface nodes
n_debris = length(xdebris);

%% Nodes

nodes = NaN*ones(n_basal+n_surface+length(zdebris), 3); % empty node matrix to be filled in

% BASAL nodes: ------------------------------------------------------------
for i = 1 : n_basal
    ii = i;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x(ii);             % x
    nodes(i,3) = basal(ii);         % z, basal
    
end;

% SURFACE nodes: ----------------------------------------------------------
for i = (n_basal) + 1 : n_basal + n_surface
    ii = i-n_basal;
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = x_surface(ii);     % x
    nodes(i,3) = surface(ii);       % z, surface
end;

% DEBRIS nodes: ----------------------------------------------------------
for i = (n_basal + n_surface) + 1 : (n_basal + n_surface) + n_debris
    ii = i-n_basal-n_surface;       %incrementing through this section only
    nodes(i,1) = i;                 % this is the node number
    nodes(i,2) = xdebris(ii);       % x location of this node
    nodes(i,3) = zdebris(ii);       % z, surface (I think this should be debris surface)
end;

%% Facets

xmin_debris = x_surface(x_surface>debris_start);
debris_intersect = getNode(nodes,xmin_debris(1,2),xmin_debris(1,3));










