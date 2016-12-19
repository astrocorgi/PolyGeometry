function node_num = getNode(nodes,x,y,z) %gets
%This function gets the node number for a given x,y input.
%It looks in the nodes matrix for the node number. Should work for both 2D
%and 3D files. The node must already exis and the x and y coordinates must
%be exactly correct
if nargin==3
    num = nodes(nodes(:,2)==x & nodes(:,nodes(:,3)==y),:);
    node_num = num(1);
else if nargin == 4
        num = nodes(nodes(:,2)==x & nodes(:,3)==y & nodes(:,4)==z,:);
        node_num = num(1);
    end
    
end

