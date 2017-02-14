function [ x,y ] = getCoord( polyA,node_num )
%This function gets the x,y coordinates for a given node input.
%It looks in the polyA matrix for the coordinates. Should work for both 2D
%and . The node must already exist and the x and y coordinates must
%be exactly correct

    index = polyA(:,1)==node_num;
    x = polyA(index,2) ;
    y = polyA(index,3);

end

