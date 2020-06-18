%%%%%%%%% For cut elements get the volume lying in the negative and positive subdomains
function [pos_area,neg_area] = sub_areas(node,e,x,y,ls,xint,yint)

xe = x(node(1:3,e));
ye = y(node(1:3,e));
lse = ls(node(1:3,e));

posId = find(lse>0.0);

if(length(posId)~=1)
    negId = find(lse<0.0);
    x_neg = [xint'; xe(negId)]; y_neg = [yint'; ye(negId)];
    neg_area = polyarea(x_neg, y_neg);
    pos_area = polyarea(xe, ye)-neg_area;
else
    x_pos = [xint'; xe(posId)]; y_pos = [yint'; ye(posId)];
    pos_area = polyarea(x_pos, y_pos);
    neg_area = polyarea(xe, ye)-pos_area;
end
