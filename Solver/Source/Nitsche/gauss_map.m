function  [gauss_x,gauss_y] = gauss_map(elem,NODE)
%UNTITLED2 Summary of this function goes here
xe = [NODE(elem.nodes).X]';
ye = [NODE(elem.nodes).Y]';
gauss_six = [0.445948490915965   0.445948490915965   
             0.445948490915965   0.108103018168070   
             0.108103018168070   0.445948490915965   
             0.091576213509771   0.091576213509771  
             0.091576213509771   0.816847572980459   
             0.816847572980459   0.091576213509771];
       
 
       gauss_x = zeros(6,1);
       gauss_y = gauss_x;
       
       for i = 1:6
           psi = gauss_six(i,1);
           eta = gauss_six(i,2);
           if(length(xe) == 3)
               N = [psi eta 1-psi-eta];
           elseif(length(xe)==6)
               N = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
           end
           gauss_x(i) = N*xe;
           gauss_y(i) = N*ye;
       end

end

