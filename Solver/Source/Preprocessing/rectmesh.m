function[x,y,node] = rectmesh(xdiv,ydiv,L,H,xc,yc)

%%%%%%%%%%%%%%% Generate a square mesh %%%%%%%%%%%%%%%%%%%%%%%%
numnodx = xdiv + 1;
numnody = ydiv + 1;

numnod = numnodx*numnody;

x = zeros(numnod,1);
y = zeros(numnod,1);

numsquare = xdiv*ydiv;
node = zeros(3,2*numsquare);
nrowtriangles = 2*xdiv;


for j=1:numnody
    x(((j-1)*numnodx + 1): j*numnodx) = linspace(-L/2+xc,L/2+xc,numnodx);
    y(((j-1)*numnodx + 1): j*numnodx) = (j-1)*(H/ydiv)-H/2+yc;
end

for ysquare=1:ydiv
    nodefirstsquare = [(ysquare-1)*numnodx+1      (ysquare-1)*numnodx+2
                       (ysquare-1)*numnodx+2      (ysquare-1)*numnodx+(numnodx+2)
                       (ysquare-1)*numnodx+(numnodx+1)      (ysquare-1)*numnodx+(numnodx+1)];
    for xsquare=1:xdiv
        xcol = 2*xsquare-1;
        ycol = (ysquare-1)*nrowtriangles + xcol;
        node(1:3,ycol:(ycol+1)) = nodefirstsquare + (xsquare-1);
    end
end

% TRI = node';
% triplot(TRI,x,y);
