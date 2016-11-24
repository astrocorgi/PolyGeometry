function [ aout, bout ] = addDebris3D( polyA,polyB,debrisStart,debrisThickness )
%ADDDEBRIS3D places a debris layer on top of a given 3D poly file. Inputs are
%the first half of a poly file (i.e. the XYZ coordinates) and the second
%half of a poly file (facet definitions), where you want the debris to
%start, and how thick you want the debris to be at the toe
% Facet notation:  
% 1 == X0 side
% 2 == X1
% 4 == Y0
% 8 == Y1
% 16 == Z0
% 32 == Z1
% 
% for example, if you have a box like this:
% 
%         /* Define 8 corner points of the box, with this order:
%          *         4 ------- 7
%          *        /         /|
%          *       /         / 6
%          *      0 ------- 3 /
%          *      |         |/
%          *      1 ------- 2
%          *
%          * Cut-out diagram with boundary flag:
%          *             4 ------- 7
%          *             | BOUNDZ1 |
%          *   4 ------- 0 ------- 3 ------- 7 ------- 4
%          *   | BOUNDX0 | BOUNDY0 | BOUNDX1 | BOUNDY1 |
%          *   5 ------- 1 ------- 2 ------- 6 ------- 5
%          *             | BOUNDZ0 |
%          *             5 ------- 6
%          */
% 
% Then you could define the BOUNDX0 facet in the following way 
% in your poly file:
% 
% 1   0   1
% 4   0   1   5   4
% 
% Confirm for yourself that when that bottom stencil gets folded up that the order of the 4 points forming the facet gives an outward normal vector according to the right-hand rule... =)
% 
% Btw, that little pic is from mesh.cxx, in case you ever need to remind yourself. 
% 


close all
facet_flag = polyB(:,4);

%for the input poly file, find the upper facets

uppers = facet_flag == 

end

