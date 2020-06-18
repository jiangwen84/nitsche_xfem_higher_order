function [psi, eta] = inverse(xpos,ypos,xe,ye)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if(length(xe) == 3)
    %No iteration is needed
    Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];
    Npar = Amat\[1 xpos ypos]';
    psi = Npar(1);
    eta = Npar(2);
    
elseif(length(xe)==6)
    
    psi = 1/3;
    eta = 1/3;
    enorm = 1;
    iter = 0;
    tol = 1e-16;
    iterMax = 1000;
    Fext = [xpos;ypos];
    
    while(enorm>tol && iter<iterMax)    
        iter = iter + 1;
        Npar = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
        Fint = [Npar*xe';Npar*ye'];
        
        Res = Fext - Fint;
             
        GN = [4*psi-1    0     4*psi+4*eta-3  4*eta    -4*eta      4-8*psi-4*eta;
              0    4*eta-1  4*eta+4*psi-3  4*psi 4-4*psi-8*eta    -4*psi];
        
        J = [xe; ye]*GN'; 
    
        delta = J\Res;
    
        psi = psi + delta(1,1);
        eta = eta + delta(2,1);
    
        if(iter==1)
            Res0 = Res;
            delta0 = delta;
        end
        enorm = (delta'*Res)/(delta0'*Res0);
    end
                  
end
end
