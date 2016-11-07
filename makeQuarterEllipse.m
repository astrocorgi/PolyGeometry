function [x2,y2,f2,elpsnodes] = makeQuarterEllipse(L,Hi0,tongH,x1,y1,x3,y3)

a = L;
b = Hi0 + tongH;
ny = length(y1);
nx = length(x3);
elpsnodes = (ny-2) + (nx-1);   %taking out (0,b) and (0,0)&(a,0)

% initialize the output variables
x2 = NaN*ones(elpsnodes,1); y2 = x2; f2 = x2;

yy1 = y1 + tongH;
yy3 = y3 + tongH;

delndx      = [];       % things this script needs to smooth
segmentchk  = 10;
nnx         = 1;

% side 2 -- elliptical 
for n = 2:ny-1 % from side 1
    x2(n-1) = sqrt(a^2 * (1-(yy1(n)/b)^2));   % x = sqrt(a^2 * (1 - (y/b)^2)
    y2(n-1) =                yy1(n)       ;   % y = y
    f2(n-1) = 32;
end;

for n = ny-1:elpsnodes % from side 3
    x2(n) =              x3(nnx+1)       ; % x = x
    y2(n) = sqrt(b^2*(1-(x3(nnx+1)/a)^2)); % y = sqrt(b^2 * (1 - (x/a)^2)
    nnx = nnx + 1;
    f2(n) = 32;
end;

y2 = y2 - tongH;

allTogether = [x2 y2 f2];
allTogether = sortrows(allTogether,1);

x2 = allTogether(:,1);
y2 = allTogether(:,2);
f2 = allTogether(:,3);

e = mean([abs(diff(y1)) ; abs(diff(x3))]);

% temporarily add low.right corner node so ellipse fills all the
% way down:
x2(end+1) = x3(1);
y2(end+1) = y3(1);
elpsnodes = elpsnodes + 1;

% smoothing
for smoothing = 1:segmentchk;
    for n = 2:elpsnodes
        distchk = sqrt((x2(n) - x2(n-1))^2 + (y2(n) - y2(n-1))^2);
        cond = distchk < .5*e; 
        if cond; delndx(end+1) = n; end;
    end;

    for n = 1:length(delndx);
        nn = delndx(n);
    
        newx = 0.5*(x2(nn+1) + x2(nn-1));
        newy = 0.5*(y2(nn+1) + y2(nn-1));
        x2(nn) = newx;
        y2(nn) = newy;
    end;
end;

x2 = x2(1:end-1);
y2 = y2(1:end-1);
elpsnodes = elpsnodes - 1;
