function flagPlot(polyA_out,polyB_out )
%FLAGPLOT takes an input polyB file and plots the flag condition for each
%node pairing.
%   Detailed explanation goes here

n = length(polyB_out(:,1));
figure(3);
hold on

for k=1:n
    flag = polyB_out(k,4);
    node_a = polyB_out(k,2);
    node_b=polyB_out(k,3);
    [x_a,y_a] = getCoord(polyA_out,node_a);
    [x_b,y_b] = getCoord(polyA_out,node_b);
    if flag==0
        plot([x_a x_b],[y_a y_b],'r-');
    elseif flag==1
        plot([x_a x_b],[y_a y_b],'m-');
    elseif flag==16
        plot([x_a x_b],[y_a y_b],'b-');
    elseif flag==32
        plot([x_a x_b],[y_a y_b],'g-');
    end %if
end


end

