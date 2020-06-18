function [ls ] = get_level_set(x,y,PARAMS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%ls = x - (PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin);

% x0 = 2.3; y0 = -0.5;
% x1 = 6.7; y1 = 0.5;
% l  = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
% phi  = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
% ls = phi/l;

ls1 = (x-0.501).^2+  (y-0.109).^2- 0.201^2;
ls2 = (x-1.501).^2+ (y+0.109).^2- 0.201^2;
ls3 = (x-2.499).^2+ (y-0.109).^2- 0.201^2;
ls4 = (x-3.501).^2+ (y+0.109).^2- 0.201^2;

ls5 = (x-4.499).^2+ (y-0.109).^2- 0.201^2;
ls6 = (x-5.501).^2+ (y+0.109).^2- 0.201^2;
ls7 = (x-6.499).^2+ (y-0.109).^2- 0.201^2;
ls8 = (x-7.501).^2+ (y+0.109).^2- 0.201^2;

ls = min([ls1 ls2 ls3 ls4 ls5 ls6 ls7 ls8],[],2);
 
end