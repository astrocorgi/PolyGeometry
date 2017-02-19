
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