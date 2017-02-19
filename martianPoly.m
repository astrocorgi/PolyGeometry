close all
clear all
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

%Convert x axis from kilometres to metres
combined(:,1) = combined(:,1)*1000;

plot(combined(:,1),combined(:,2)); %plot result
%combined = deleteRow(combined,24468.051803350210); %select point for removal
combined(54,:) = [];

plot(combined(:,1),combined(:,2)); %plot result
%combined = deleteRow(combined,24006.390448605742); %select point for removal
combined(53,:) = [];
plot(combined(:,1),combined(:,2)); %plot result

figure
index = combined(:,1) < 15000;
combined_crop = combined(index,:);
combined_crop(33,:) = []; 
combined_crop(32,:) = []; 
combined_crop(31,:) = []; 

combined_crop(:,2) = combined_crop(:,2)-0.03*combined_crop(:,1)+500;

toe_insert = [14250,375;14750,150];
combined_crop = [combined_crop(1:30,:); toe_insert(:,:); combined_crop(31:end,:)]; 


plot(combined_crop(:,1),combined_crop(:,2));


%% assigning node numbers


nodes = [0:length(combined_crop(:,1))-1]';
x = combined_crop(:,1);
y = combined_crop(:,2);

polyA = [nodes,x,y];
polyB = zeros(length(x),4);

%Loop everything around
for k = 0:length(x)-2
    polyB(k+1,1) = k;
    polyB(k+1,2) = k;
    polyB(k+1,3) = k+1;
    if k < 32 %derived from plot
        polyB(k+1,4) = 32;
    else
        polyB(k+1,4) = 16;
    end
end

%connecting the end to the start
end_node = polyA(end,1);
polyB(end,1) = polyB(end-1,1)+1;
polyB(end,2) = end_node;
polyB(end,3) = 0;
polyB(end,4) = 1;

flagPlot(polyA,polyB);

save polyfile_martian1

%% Now add the debris

addDebris(polyA,polyB,0,50,'euripus_ideal1.poly')


