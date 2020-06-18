function [uex,FluxEx,gradU_voigt] = exactsolution(x,y,domain,PARAMS)

[E,~,~] = mat_param(PARAMS);
Xmax = PARAMS.Xmax; Xmin = PARAMS.Xmin; xint = PARAMS.xint;
Ymax = PARAMS.Ymax; Ymin = PARAMS.Ymin; yint = PARAMS.yint;

uex = []; FluxEx = []; gradU_voigt = [];


%%%%%%%% Analytical solutions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

if(domain==1)
    r = sqrt(x^2 + y^2);
  %  if(PARAMS.ndof==2)        
        uex = 0;
       % gradU_voigt = [y^2/r^3;x^2/r^3;1/r];
        gradUex = [0 0];
        %if(x==Xmin)
        %    FluxEx = E(domain)*gradUex*[-1;0];
        %elseif(x==Xmax)
        %    FluxEx = E(domain)*gradUex*[1;0];
        %end
        
        FluxEx = 0;%E(domain)*gradUex*[x/r;y/r];
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
else
    %r = sqrt(x^2 + y^2);
    %ur = 1 + log(2*r);
  %  if(PARAMS.ndof==2)
        uex = 0;
%        gradU_voigt = [x^2/r^3 + ur*y^2/r^3;y^2/r^3 + ur*x^2/r^3; (1+ur)/r];
%        gradUex = [x^2/r^3 + ur*y^2/r^3, (1-ur)*x*y/r^3; (1-ur)*x*y/r^3, y^2/r^3 + ur*x^2/r^3];
         gradUex = [0 0];%[x/r^2 y/r^2];
        FluxEx = 0;%E(domain)*gradUex*[x/r;y/r];
end

%        if(x==Xmin)
%            FluxEx = E(domain)*gradUex*[-1;0];
%        elseif(x==Xmax)
%            FluxEx = E(domain)*gradUex*[1;0];
%        end
%     else      
%         uex = c1_2 + c2_2*x;
%         gradUex = [c2_2;0.0];
%         gradU_voigt = gradUex;
%         if(x==Xmin)
%             FluxEx = E(domain)*gradUex'*[-1;0];
%         elseif(x==Xmax)
%             FluxEx = E(domain)*gradUex'*[1;0];
%         end
%    end



end
