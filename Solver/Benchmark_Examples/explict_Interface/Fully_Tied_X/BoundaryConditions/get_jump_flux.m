function [t_jump] = get_jump_flux(x,y,elm,sib,PARAMS)

%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% elmDomain = elm.domain;
% sibDomain = sib.domain;
% 
% [~,tElm,~] = exactsolution(x,y,elmDomain,PARAMS);
% [~,tSib,~] = exactsolution(x,y,sibDomain,PARAMS);
% 
% if(elmDomain==1)
%   t_jump = tElm - tSib;
% else
%    t_jump = tSib - tElm;
% end

if(elm.domain==2)
    normals = elm.normals;
else
    normals = sib.normals;
end
    
%gradUex = [-cos(x)*cos(y), sin(x)*sin(y)];

t_jump = -cos(x)*cos(y)*normals(1) + sin(x)*sin(y)*normals(2);



end

