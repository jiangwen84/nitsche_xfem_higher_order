function [int_flag] = intersection_test2(node,e,x,y,PARAMS)
%%%%%%% Flags indicating intersection with elements %%%%%
int_flag = 0;    %%%% Flag for a cut element

p0 = [x(node(1,e)) y(node(1,e))];
p1 = [x(node(2,e)) y(node(2,e))];
p2 = [x(node(3,e)) y(node(3,e))];

p = [p0;p1;p2;...
    0.5*p0+0.5*p1;0.5*p1+0.5*p2;0.5*p0+0.5*p2;...
    0.25*p0+0.75*p1;0.75*p0+0.25*p1;...
    0.25*p1+0.75*p2;0.75*p1+0.25*p2;...
    0.25*p2+0.75*p0;0.75*p2+0.25*p0];

xp = p(:,1);
yp = p(:,2);

lsp = get_level_set(xp,yp,PARAMS);

for i=2:size(xp,1)
    schange=lsp(i-1)*lsp(i);   %%%%% Sign change in the level set
    if (schange<=0)
        int_flag = 1;
    end
end
