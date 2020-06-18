function [uex,FluxEx,gradU_voigt] = exactsolution(x,y,domain,PARAMS)

[E,~,~] = mat_param(PARAMS);
Xmax = PARAMS.Xmax; Xmin = PARAMS.Xmin; xint = PARAMS.xint;
Ymax = PARAMS.Ymax; Ymin = PARAMS.Ymin; yint = PARAMS.yint;

uex = []; FluxEx = []; gradU_voigt = [];

if(domain==1)
    r2 = x^2 + y^2;
    uex = r2;
    gradUex = [2*x 2*y];
    FluxEx = E(domain)*gradUex*[x/sqrt(r2);y/sqrt(r2)];
    
    gradU_voigt = gradUex';
elseif(domain==2)
    r2 = x^2 + y^2;
    uex = 0.1*r2^2 - 0.01*log(2*sqrt(r2));
    gradUex = [0.2*r2*2*x - 0.01*x/r2 0.2*r2*2*y - 0.01*y/r2];
    FluxEx = E(domain)*gradUex*[x/sqrt(r2);y/sqrt(r2)];
    gradU_voigt = gradUex';
end

end
