function [ x,y ] = getCoord( polyA,node_num )
%This function gets the x,y coordinates for a given node input.
%It looks in the polyA matrix for the coordinates. Only works in 2D for now.
%The node must already exist. Vector input is acceptable.
    x = zeros(length(node_num),1);
    y = zeros(length(node_num),1);

    for k = 1:length(node_num)
        index = polyA(:,1)==node_num(k);
        x(k) = polyA(index,2);
        y(k) = polyA(index,3);
    end
end

