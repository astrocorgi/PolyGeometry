function [ PolyA_Out, PolyB_Out ] = addDebris(polyA,polyB,debrisStart,debrisThickness)
%addDebris places a debris layer on top of a given poly file. Inputs are
%the first half of a poly file (i.e. the X-Y coordinates) and the second
%half of a poly file (facet definitions), where you want the debris to
%start, and how thick you want the debris to be at the toe
%16 = bottom facet, 32 = top facet
%   Detailed explanation goes here
    
    close all
    facet_flag = polyB(:,4);
    %for the input poly file, find all the upper facets
    uppers = facet_flag == 32;
    disp(uppers);
    upper_polyA = polyA(uppers,:);
    disp(upper_polyA);
    
    %allocate memory for new y-coord array
    y = zeros(length(uppers),1);
    %use linear equation for debris thickness change
    xmax = max(upper_polyA(:,2));
    delta_x = xmax - debrisStart;
    delta_y = debrisThickness;
    m = delta_y/delta_x; %slope
    
    %create new points
    index = upper_polyA(:,2) > debrisStart; %if the x-coordinate of the point is greater than where we want to start the debris layer...
    upper_xfilt_polyA = upper_polyA(index,:); %then store it in a new array
    x = upper_xfilt_polyA(:,2); %get the x coordinates of the upper facets 
    y = m*x + upper_xfilt_polyA(:,3); %y = mx + b to create the new y coordinates for the debris layer
    
    %sort the x,y in ascending x
    A = [x,y];
    [X,Y]=sort(A(:,1));
    B=A(Y,:);
    x = B(:,1);
    y = B(:,2);

    num_point = (polyA(end,1)+1:polyA(end,1)+length(y))'; %create the point numbers
    
    polyA_append = [num_point,x,y]; %Combine to add to poly file
    disp(polyA_append);
    plot(x,y,'*');
    hold on
    plot(polyA(:,2),polyA(:,3));
    xlabel('X');
    ylabel('Y');
    legend('Debris layer','Initial poly geometry');
    title('Plot of the new poly x-y coordinates');
    
    %change all "32" facets greater than the debris start
    num = 1;
    for i=1:length(polyB(:,1))
        if polyB(i,4)==32
            %check if the respective point in polyA is greater than the
            %debris start. Match by point number
            point_number = polyB(i,2);
            row_index = polyA(:,1) == point_number;
            x_val = polyA(row_index,2);
            if x_val > debrisStart
                %change the 32 flag to an internal flag (0)
                polyB(i,4) = 0;
                x_pointnum(num,:) = [x_val,point_number]; %collect the x and point number values into one array to use later
                num = num+1;
            end
        end
    end
    
    
    
    %%Create the new polyB segment
    %Find the points intersecting the debris layer
    xmin_index = x_pointnum(:,1)==min(x_pointnum(:,1));
    intersect1 = x_pointnum(xmin_index,2);
    disp('debris startpoint');
    disp(intersect1);
    
    xmax_index = x_pointnum(:,1)==max(x_pointnum(:,1));
    intersect2 = x_pointnum(xmax_index,2);
    disp('debris endpoint');
    disp(intersect2);
    %disp(polyB);
    
    polyB_append = zeros(length(polyA_append(:,1))+1,4); %allocate new array
    point_num = polyB(end,1)+1;
    
    %Create the facet array for the debris layer
    for i = 1:length(polyA_append(:,1))+1
        if (i==1) %the first facet starts at the debris/ice intersection
            polyB_append(i,:) = [point_num,polyA_append(i,1),intersect1,32];
        elseif (i==length(polyA_append(:,1))+1)
            polyB_append(i,:) = [point_num,polyA_append(i-1,1),intersect2,32];
        else %otherwise the node connects to 
            polyB_append(i,:) = [point_num,polyA_append(i,1),polyA_append(i,1)-1,32];
        end
        point_num = point_num+1;
    end
    PolyA_Out = [polyA;polyA_append];
    PolyB_Out = [polyB;polyB_append];
    %disp('polyA new');
    %disp(PolyA_Out);
    %disp('polyB new');
    %disp(PolyB_Out);
end

