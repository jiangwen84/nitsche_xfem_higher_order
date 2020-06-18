function [fe] = get_bodyforce(elem,NODE,ndof,nlink,PARAMS,weight,gp)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

       fe = zeros(ndof*nlink,1);
      
       x0 = [NODE(elem.nodes).X];
       y0 = [NODE(elem.nodes).Y];
       x0 = x0'; y0 = y0';
       
       for i = 1:6
           psi = gp(i,1);
           eta = gp(i,2); 
           
           if(length(x0) == 3)
               Npar = [psi eta 1-psi-eta];
           elseif(length(x0)==6)
               Npar = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
           end
           
           Posx = Npar*x0; Posy = Npar*y0;
       
           [fBar] = bodyforce(Posx,Posy,elem.domain,PARAMS);
           
           fe = fe + weight(i)*Npar'*fBar;
       end
end