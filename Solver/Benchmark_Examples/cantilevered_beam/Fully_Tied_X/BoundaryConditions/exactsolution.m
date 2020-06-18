function [uex,FluxEx,gradU_voigt] = exactsolution(x,y,domain,PARAMS)

[E,~,~] = mat_param(PARAMS);
Xmax = PARAMS.Xmax; Xmin = PARAMS.Xmin; xint = PARAMS.xint;
Ymax = PARAMS.Ymax; Ymin = PARAMS.Ymin;

uex = []; FluxEx = []; gradU_voigt = [];

L = Xmax - Xmin;
D = Ymax - Ymin;
I = D^3/12;

nu = 0.0;
P = 1;

uex = [P*y/(6*E(1)*I)*((6*L-3*x)*x+(2+nu)*(y^2-D^2/4)); -P/(6*E(1)*I)*(3*nu*y^2*(L-x)+(4+5*nu)*D^2*x/4+(3*L-x)*x^2)];
gradU_voigt =[P*y/(6*E(1)*I)*(6*L-6*x);-P/(6*E(1)*I)*(6*nu*y*(L-x));P/(6*E(1)*I)*((6*L-3*x)*x+(2+nu)*(y^2-D^2/4))+P*y/(6*E(1)*I)*((2+nu)*2*y)+ -P/(6*E(1)*I)*(-3*nu*y^2+(4+5*nu)*D^2/4+6*L*x-3*x^2)];

end
