function [u_jump] = get_jump_interface(x,y,elmDomain,sibDomain,PARAMS)

[uElm,~,~] = exactsolution(x,y,elmDomain,PARAMS);
[uSib,~,~] = exactsolution(x,y,sibDomain,PARAMS);
u_jump = uElm-uSib;
