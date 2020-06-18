function [t_jump] = get_jump_flux(x,y,elm,sib,PARAMS)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[~,tElm,~] = exactsolution(x,y,elmDomain,PARAMS);
%[~,tSib,~] = exactsolution(x,y,sibDomain,PARAMS);

%if(elmDomain==1)
%  t_jump = tElm - tSib;
%else
%   t_jump = tSib - tElm;
%end


if(elm.domain==2)
    normals = elm.normals;
else
    normals = sib.normals;
end
r = x^2 + y^2;
t_jump = ( 4*r - 0.1/r - 2 )*( x*normals(1) + y*normals(2) );




end

