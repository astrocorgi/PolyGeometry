function [ PolyA_Out, PolyB_Out ] = addDebrisMartian(polyA,polyB,debrisStart,debrisThickness,output_file)
%addDebris places a debris layer on top of a given poly file. Inputs are
%the first half of a poly file (i.e. the X-Y coordinates) and the second
%half of a poly file (facet definitions), where you want the debris to
%start, and how thick you want the debris to be at the toe
%16 = bottom facet, 32 = top facet
%   Detailed explanation goes here
    
    close all
    facet_flag = polyB(:,4);
    uppers = facet_flag==32;
    num = sum(uppers);
    
    polyB(facet_flag==32,4) = 0; %set all upper facets to 0's
    internal_nodes = polyB(facet_flag==32,2);
    [internal_x, internal_y]= getCoord(polyA,internal_nodes);
    x_pointnum = [internal_x,internal_nodes]; %collect the x and point number values for the top nodes into one array to use later
        
    %use linear equation for debris thickness increase
    xmax = max(internal_x);
    delta_x = xmax - debrisStart;
    delta_y = debrisThickness;
    m = delta_y/delta_x; %slope
    
    %create new points
    y_new = m*(internal_x-debrisStart) + internal_y+20; %y = mx + b to create the new y coordinates for the debris layer
    
    %add one extra point for debris-layer hitting ground in front of ice
    internal_x(length(internal_x)+1) = max(polyA(:,2))+debrisThickness*2; %x = glacier toe+debris thickness plus a little bit extra (this is arbitrary, I don't think anyone has characterized debris accumulation at DCG toes)
    y_new(length(y_new)+1) = min(polyA(:,3)); %this assumes the toe is the lowest point. this is almost always true, but still might bug out later
    
    %sort the x,y in ascending x
    A = [internal_x,y_new];
    [X,Y] = sort(A(:,1));
    B = A(Y,:);
    internal_x = B(:,1);
    y_new = B(:,2);

    num_point = (polyA(end,1)+1:polyA(end,1)+length(y_new))'; %create the point numbers
    
    polyA_append = [num_point,internal_x,y_new]; %Combine to add to poly file
    
    %add some padding to the front so the debris layer is extra thick there
    %NOTE: AWKWARD HARD CODING HERE
    polyA_append(end-10:end,2) = polyA_append(end-10:end,2) + 200;
    polyA_append(end-10:end-1,3) = polyA_append(end-10:end-1,3) - 15; %The front coords were too tall
    
    disp(polyA_append);
    plot(polyA_append(:,2),y_new,'*');
    hold on
    plot(polyA(:,2),polyA(:,3));
    xlabel('X');
    ylabel('Y');
    legend('Debris layer','Initial poly geometry');
    title('Plot of the new poly x-y coordinates');
    
    %change all "32" facets to 0 facets
%     num = 1;
%     for i=1:length(polyB(:,1))
%         if polyB(i,4)==32
%             %check if the respective point in polyA is greater than the
%             %debris start. Match by point number
%             point_number = polyB(i,2);
%             row_index = polyA(:,1) == point_number;
%             x_val = polyA(row_index,2);
%             %if x_val > debrisStart
%                 %change the 32 flag to an internal flag (0)
%                 polyB(i,4) = 0;
%                 x_pointnum(num,:) = [x_val,point_number]; %collect the x and point number values into one array to use later
%                 num = num+1;
%             %end
%         end
%     end
        
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
        if (i==1) %the first facet starts at the debris/ice intersection and is a 1
            polyB_append(i,:) = [point_num,polyA_append(i,1),intersect1,1];
        elseif (i==length(polyA_append(:,1))+1) %connecting the ice to the debris toe, is a bottom facet
            polyB_append(i,:) = [point_num,polyA_append(i-1,1),intersect2+1,16];
        else %otherwise the node connects to 
            polyB_append(i,:) = [point_num,polyA_append(i,1),polyA_append(i,1)-1,32];
        end
        point_num = point_num+1;
    end
    
    
    PolyA_Out = [polyA;polyA_append];
    PolyB_Out = [polyB;polyB_append];
    
    flagPlot(PolyA_Out,PolyB_Out);
    %disp('polyA new');
    %disp(PolyA_Out);
    %disp('polyB new');
    %disp(PolyB_Out);
    
    
    %Find points for material determination
    internal_index = PolyB_Out(:,4) == 0; %all the internal boundaries
    internals = PolyB_Out(internal_index,:);
    internal_midpoint = floor(length(internals(:,1))/2);
    internal_node = internals(internal_midpoint,2);
    [x_material,y_internal] = getCoord(PolyA_Out,internal_node);
    y_mat1 = y_internal - 0.01;
    y_mat2 = y_internal + 0.01;
    
    %Write poly file
    fileID = fopen(output_file,'w');
    fprintf(fileID,'%d %d %d %d \n', length(PolyA_Out(:,1)),2,0,0);%coordinate header
    format_spec = '%d %.15g %.15g \n';
    for k = 1:length(PolyA_Out(:,1))
        fprintf(fileID,format_spec,PolyA_Out(k,1),PolyA_Out(k,2),PolyA_Out(k,3));
    end
    fprintf(fileID,'\n'); %add line break
    fprintf(fileID,'%d %d \n',length(PolyB_Out(:,1)),1); %facet header
    format_spec = '%d %d %d %d \n';
    for k=1:length(PolyB_Out(:,1))
        fprintf(fileID,format_spec,PolyB_Out(k,1),PolyB_Out(k,2),PolyB_Out(k,3),PolyB_Out(k,4));
    end
    fprintf(fileID,'%d \n',0); %material header
    fprintf(fileID,'%d \n',2); %number of materials = 2
    
    fprintf(fileID,'%d %.15g %.15g %d %d \n',0,x_material,y_mat1,0,1e4);
    fprintf(fileID,'%d %.15g %.15g %d %d \n',1,x_material,y_mat2,1,1e4);
    
    
end

