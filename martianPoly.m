close all
clear 
clc
mola = load('cms_mola_transect1.txt');
subsurf = load('cms_subsurface_s1_smoothed_transect1.txt');

mola_x = distdim(mola(:,1),'deg','km','mars');
subsurf_x = distdim(subsurf(:,1),'deg','km','mars');

figure
plot(mola_x,mola(:,2));
hold on
plot(subsurf_x,subsurf(:,2));

figure
subsurf = [subsurf_x,subsurf(:,2)];
subsurf_flip = flipud(subsurf);

combined = [mola_x, mola(:,2);subsurf_flip(:,1),subsurf_flip(:,2)];

plot(combined(:,1),combined(:,2));

combined = deleteRow(combined,24.468051803350210);


plot(combined(:,1),combined(:,2));

combined = deleteRow(combined,24.006390448605742);


function new_matrix = deleteRow(input_matrix,x_value)
    %this function removes a row from the combined x,y data matrix
    %according to an input x-value. Used to trim the measured geometry down
    %to something worth inputting to the model. Could potentially bug out
    %if two y values share an x value, but that seems really unlikely given
    %the number of sig figs and by how we set up the arrays.
    
    index = input_matrix(:,1)==x_value;
    input_matrix(index,:) = []; 
    new_matrix = input_matrix;
end