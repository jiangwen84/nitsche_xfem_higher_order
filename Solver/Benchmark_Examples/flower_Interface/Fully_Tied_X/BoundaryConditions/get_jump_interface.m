function [u_jump] = get_jump_interface(x,y,elmDomain,sibDomain,PARAMS)

[uElm,~,~] = exactsolution(x,y,elmDomain,PARAMS);
[uSib,~,~] = exactsolution(x,y,sibDomain,PARAMS);
u_jump = uElm-uSib;

%r = x^2 + y^2;

%u_jump = 0.1*r^2 - 0.01*log(2*sqrt(r) ) - r;
end
