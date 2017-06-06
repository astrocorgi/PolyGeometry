% this script loads the data and runs the functions needed to make a poly
% geometry file for dynearthsol2D input. It takes the picks from SHARAD,
% loops em together, cleans it up a little, then makes it into a poly file.
% Finally, the debris is added on top.


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

toe_insert = [14250,360;14750,150];
combined_crop = [combined_crop(1:30,:); toe_insert(:,:); combined_crop(31:end,:)]; 


%% assigning node numbers

x = combined_crop(:,1);
y = combined_crop(:,2); %fff

xx = interp(x,5);
yy = interp(y,5);

plot(xx,yy,'b*')


polyA = [xx,yy]; %incomplete polyA


%Delete and adjust specific points that look bad
polyA(168,:) = [];
polyA(168,:) = [];
polyA(168,:) = [];
polyA(168,:) = [];
polyA(167,:) = [];
polyA(166,:) = [];
polyA(165,:) = [];
polyA(164,:) = [];
polyA(163,:) = [];
polyA(162,:) = [];
polyA(142,2) = polyA(142,2)+2;
polyA(143,2) = polyA(143,2)+3;
polyA(144,2) = polyA(144,2)+4;
polyA(145,2) = polyA(145,2)+5;
polyA(146,2) = polyA(146,2)+4;
polyA(147,2) = polyA(147,2)+3;
polyA(147:159,2) = polyA(147:159,2)-4;

plot(polyA(:,1),polyA(:,2),'ro');

polyB = zeros(length(polyA(:,1)),4);

nodes = [0:length(polyA(:,1))-1]';
polyA = [nodes, polyA(:,1), polyA(:,2)];

%Loop everything around
for k = 0:length(polyA(:,1))-2
    polyB(k+1,1) = k;
    polyB(k+1,2) = k;
    polyB(k+1,3) = k+1;
    if k < 160 %derived from plot
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

%subtracting a parabola from the bottom layer to make a compressive regime
%points 160 to 719 are on the bottom
l = 719-160;
dist = max(x);
interval = dist/l;
x_parabola = 0:interval:dist;
y_parabola = -0.0000050*(x_parabola-14916/2).^2+278.1;
polyA(161:720,3) = polyA(161:720,3) - y_parabola';

%fixing individual basal points that look weird
polyA(162,3) = polyA(162,3) - 10;
polyA(163,3) = polyA(163,3) - 8;
polyA(164,3) = polyA(164,3) - 5;
polyA(165,3) = polyA(165,3) - 5;
polyA(166,3) = polyA(166,3) - 6;
polyA(167:186,3) = polyA(167:186,3) - 2;
polyA(187,3) = polyA(187,3) - 1;

flagPlot(polyA,polyB);

save polyfile_martian_comp

%% Now add the debris

addDebrisMartian(polyA,polyB,0,5,'euripus_ideal_comp.poly')
