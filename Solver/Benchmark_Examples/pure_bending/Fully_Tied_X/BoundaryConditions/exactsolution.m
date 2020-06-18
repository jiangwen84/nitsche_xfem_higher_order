function [uex,FluxEx,gradU_voigt] = exactsolution(x,y,domain,PARAMS)

[E,~,~] = mat_param(PARAMS);
Xmax = PARAMS.Xmax; Xmin = PARAMS.Xmin; xint = PARAMS.xint;
Ymax = PARAMS.Ymax; Ymin = PARAMS.Ymin;
%%%%%%%%% Analytical solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c2_1 = 1/(Xmax-Xmin); c1_1 = -Xmin/(Xmax-Xmin);
% c2_2 = E(1)/E(2)*c2_1;
% c1_2 = c1_1+(1-E(1)/E(2))*c2_1*(xint*(Xmax-Xmin)+Xmin);
% 
% uex = []; FluxEx = []; gradU_voigt = [];
% if(domain==1)
%     if(PARAMS.ndof==2)        
%         uex = [c1_1 + c2_1*x; 0.0*y];
%         gradU_voigt = [c2_1;0;0];
%         gradUex = [c2_1, 0; 0, 0];
%         if(x==Xmin)
%             FluxEx = E(domain)*gradUex*[-1;0];
%         elseif(x==Xmax)
%             FluxEx = E(domain)*gradUex*[1;0];
%         end        
%     else        
%         uex = c1_1 + c2_1*x;
%         gradUex = [c2_1;0.0];
%         gradU_voigt = gradUex;
%         if(x==Xmin)
%             FluxEx = E(domain)*gradUex'*[-1;0];
%         elseif(x==Xmax)
%              FluxEx = E(domain)*gradUex'*[1;0];
%         end
%     end
% elseif(domain==2)
%     if(PARAMS.ndof==2)
%         uex = [c1_2 + c2_2*x; 0.0*y];
%         gradU_voigt = [c2_2;0;0];
%         gradUex = [c2_2, 0; 0, 0];
%         if(x==Xmin)
%             FluxEx = E(domain)*gradUex*[-1;0];
%         elseif(x==Xmax)
%             FluxEx = E(domain)*gradUex*[1;0];
%         end
%     else      
%         uex = c1_2 + c2_2*x;
%         gradUex = [c2_2;0.0];
%         gradU_voigt = gradUex;
%         if(x==Xmin)
%             FluxEx = E(domain)*gradUex'*[-1;0];
%         elseif(x==Xmax)
%             FluxEx = E(domain)*gradUex'*[1;0];
%         end
%     end

uex = []; FluxEx = []; gradU_voigt = [];

L = Xmax - Xmin;
H = Ymax - Ymin;
p = 1/H;
nu = 0;
if(domain==1)
         
        uex = [2*p*x*y/E(1); -p*(x^2-nu*y^2)/E(1)];
        gradU_voigt = [2*p*y/E(1); 2*p*nu*y/E(1); 0];
        gradUex = [2*p*y/E(1), 2*p*x/E(1); -2*p*x/E(1), 2*nu*p*y/E(1) ];
        if(y==Ymin)
            FluxEx = E(domain)*gradUex*[0;-1];
        elseif(y==Ymax)
            FluxEx = E(domain)*gradUex*[0;1];
        end        
elseif(domain==2)
        uex = [2*p*x*y/E(2); -p*(x^2-nu*y^2)/E(2)];
        gradU_voigt = [2*p*y/E(2); 2*p*nu*y/E(2); 0];
        gradUex = [2*p*y/E(2), 2*p*x/E(2); -2*p*x/E(2), 2*nu*p*y/E(2) ];
        if(y==Ymin)
            FluxEx = E(domain)*gradUex*[0;-1];
        elseif(y==Ymax)
            FluxEx = E(domain)*gradUex*[0;1];
        end
end
end
