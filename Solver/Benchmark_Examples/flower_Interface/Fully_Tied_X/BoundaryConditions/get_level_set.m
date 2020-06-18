function [ls ] = get_level_set(x,y,PARAMS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global pp

position = [x y];
maxdist(1:size(position,1)) = 1000;
for k = 1:size(position,1)
    for j=1:size(pp,1)
        dist = sqrt((position(k,1)-pp(j,1))*(position(k,1)-pp(j,1)) + (position(k,2)-pp(j,2))*(position(k,2)-pp(j,2)));
        if dist<maxdist(k)
            maxdist(k) = dist;
        end
    end
end

yc=0.02*sqrt(5);
xc=0.02*sqrt(5);

for k=1:size(position,1)
    %compute theta at point
    thetap = atan2((position(k,2)-yc),(position(k,1)-xc));
    %compute rad at point
    radp = sqrt((position(k,1)-xc)*(position(k,1)-xc) + (position(k,2)-yc)*(position(k,2)-yc));
    %radius of actual point
    rad_act = 0.5+0.2*sin(5*thetap);
    if (radp > rad_act)
        ls(k,1) = maxdist(k);
    else
        ls(k,1) = -maxdist(k);
    end
    
end

end

