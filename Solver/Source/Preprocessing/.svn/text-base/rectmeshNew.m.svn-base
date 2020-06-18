function[x,y,node] = rectmeshNew(xdiv,ydiv,L,H,xc,yc)
%xdiv = 23; ydiv = 23; L = 1; H = 1; xc = 0.5; yc = 0.5;
%%%%%%%%%%%%%%% Generate a square mesh %%%%%%%%%%%%%%%%%%%%%%%%
numnodx = xdiv + 1;
numnody = ydiv + 1;

numnod = numnodx*numnody;

x = zeros(numnod,1);
y = zeros(numnod,1);

numsquare = (6*numnodx-2)*ydiv;
node = zeros(3,2*numsquare);
nrowtriangles = 2*(6*numnodx-2);


for j=1:numnody
%     x(((j-1)*numnodx + 1): j*numnodx) = linspace(-L/2+xc,L/2+xc,numnodx);
    x(((j-1)*(6*numnodx-1) + 1): j*(6*numnodx-1)) = union(linspace(-L/2+xc,(xc-0.06),(5.5)*numnodx),linspace((xc-0.06),L/2+xc,(1/2)*numnodx));
%     y(((j-1)*numnodx + 1): j*numnodx) = (j-1)*(H/ydiv)-H/2+yc;
    y(((j-1)*(6*numnodx-1) + 1): j*(6*numnodx-1)) = (j-1)*(H/ydiv)-H/2+yc;
end

for ysquare=1:ydiv
    nodefirstsquare = [(ysquare-1)*(6*numnodx-1)+1      (ysquare-1)*(6*numnodx-1)+2
                       (ysquare-1)*(6*numnodx-1)+2      (ysquare-1)*(6*numnodx-1)+(6*numnodx+1)
                       (ysquare-1)*(6*numnodx-1)+(6*numnodx)      (ysquare-1)*(6*numnodx-1)+(6*numnodx)];
    for xsquare=1:(6*numnodx-2)
        xcol = 2*xsquare-1;
        ycol = (ysquare-1)*nrowtriangles + xcol;
        node(1:3,ycol:(ycol+1)) = nodefirstsquare + (xsquare-1);
    end
end

% TRI = node';
% triplot(TRI,x,y);
