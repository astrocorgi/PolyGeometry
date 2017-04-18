%% Generate cookie dough model for the NPLD

close all
layer_thickness = 5;
ice_thickness = 30;
total_thickness = 3000;

%Import x,y MOLA outline
load npld_xy
x = DistanceKm * 1000; %distance in metres
y = MOLA512ppdElevationabove73N;


plot(x,y)
title('Original profile');
xlabel('Distance, m');
ylabel('Elevation, m');

fit_coeff = [-2.11002559631190e-54,1.47864923144361e-47,-4.42231384212260e-41,7.35502315970651e-35,-7.42105431419209e-29,4.64159898149896e-23,-1.75848368158209e-17,3.73402853309855e-12,-3.71695936188702e-07,0.0147665674494390,-4676.62516115479];
y_fit = fit_coeff(1)*x.^10 + fit_coeff(2)*x.^9 + fit_coeff(3)*x.^8 + fit_coeff(4)*x.^7 + fit_coeff(5)*x.^6 + fit_coeff(6)*x.^5 + fit_coeff(7)*x.^4 + fit_coeff(8)*x.^3 + fit_coeff(9)*x.^2 + fit_coeff(10)*x + fit_coeff(11);
figure
plot(x,y_fit);

% Trim the ends off by x-values
x_crop = x(x>9.422e4);
y_crop = y_fit(x>9.422e4);
x_crop = x_crop(x_crop<1.276e6);
y_crop = y_crop(x_crop<1.276e6);

%Trim the end off by y-value (making it level)
x_crop = x_crop(y_crop > -4547);
y_crop = y_crop(y_crop > -4547);
plot(x_crop,y_crop);

y_crop = y_crop - min(y_crop); %changing elevation to 0

%% Create the polyA file

%assign node numbers
polyA = [[0:length(x_crop)-1]',x_crop,y_crop];

%% Create the polyB file

polyB = zeros(length(x_crop),4);

%loop everything around to make facets
for k = 0:length(polyA(:,1))-2
    polyB(k+1,1) = k;
    polyB(k+1,2) = k;
    polyB(k+1,3) = k+1;
    polyB(k+1,4) = 32; %all these facets are top facets
end

%Connecting the end to the start
end_node = polyA(end,1);
polyB(end,1) = polyB(end-1,1)+1;
polyB(end,2) = end_node;
polyB(end,3) = 0;
polyB(end,4) = 0; %facet flag

figure
flagPlot(polyA,polyB);

%% Add the internal facets

%find nodes one ice thickness up from bottom
%first_node = y_crop nearest to

%n_layers = max(y_crop)/(ice_thickness+layer_thickness);
%for k = 1:n_layers;
    
%end














