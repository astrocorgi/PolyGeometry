% This script creates a poly file for the interpolated version of
% arolla.poly. Load the "polyA_interp.mat" file into the workspace and then
% run. The output of this script is run through the addDebris function in
% interp_debris_runscript

nodes = polyA_interp(:,1);
x = polyA_interp(:,2);
y = polyA_interp(:,3);

polyB_interp = zeros(length(x),4);

%Loop everything around
for k = 0:length(x)-2
    polyB_interp(k+1,1) = k;
    polyB_interp(k+1,2) = k;
    polyB_interp(k+1,3) = k+1;
    if k < 498 %derived from plot
        polyB_interp(k+1,4) = 16;
    else
        polyB_interp(k+1,4) = 32;
    end
end

%connecting the end to the start
end_node = polyA_interp(end,1);
polyB_interp(end,1) = polyB_interp(end-1,1)+1;
polyB_interp(end,2) = end_node;
polyB_interp(end,3) = 0;
polyB_interp(end,4) = 32;

flagPlot(polyA_interp,polyB_interp);

save polyA_interp