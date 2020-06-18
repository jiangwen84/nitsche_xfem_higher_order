function [E,nu,yieldParam] = mat_param(PARAMS)
%%%%%%%%% Material parameters for the positive and negative subdomains %%%%
E = [1 10]; nu = zeros(2,1);

%%%%% yieldParam format for each row -> 
%            [Domain a, Domain b, Characteristic property for interface ab]
%     Characteristic property may be yield traction or
%     Coulomb's frictional coefficient
yieldParam = [1 1 1e16; 1 2 1e16; 2 2 1e16];