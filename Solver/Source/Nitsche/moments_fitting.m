function [ weight ]= moments_fitting(xe,ye,gauss_x,gauss_y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

ls = length(xe);
            
         A = zeros(6);
         inter = zeros(6,1);
         
         A(1,:) = ones(1,6);
         A(2,:) = gauss_x;
         A(3,:) = gauss_y;
         A(4,:) = gauss_x.*gauss_y;
         A(5,:) = gauss_x.*gauss_x;
         A(6,:) = gauss_y.*gauss_y;
       
      if(ls==4)
                
          gp = [-sqrt(3)/3 -sqrt(3)/3
                        sqrt(3)/3 -sqrt(3)/3
                        sqrt(3)/3  sqrt(3)/3
                       -sqrt(3)/3  sqrt(3)/3];
        
         w_4 = [1 1 1 1]; 
       
        for i = 1:4
            psi = gp(i,1);
            eta = gp(i,2);
            N = 0.25 * [(1-psi)*(1-eta)  (1+psi)*(1-eta)  (1+psi)*(1+eta)  (1-psi)*(1+eta)];                   % shape functions matrix  
            GN    = 0.25 * [eta-1  1-eta   1+eta   -eta-1;
                        psi-1  -psi-1  1+psi    1-psi];
            Jacobian = GN*[xe ye];       % compute Jacobian matrix 
            detJ  = det(Jacobian);     % Jacobian         
            inter(1) = inter(1) + w_4(i)*N*detJ*[1, 1, 1, 1]';
            inter(2) = inter(2) + w_4(i)*N*detJ*xe;
            inter(3) = inter(3) + w_4(i)*N*detJ*ye;
            inter(4) = inter(4) + w_4(i)*xe'*(N'*N)*ye*detJ;
            inter(5) = inter(5) + w_4(i)*xe'*(N'*N)*xe*detJ;
            inter(6) = inter(6) + w_4(i)*ye'*(N'*N)*ye*detJ;
        end
        
        weight = A\inter;
        
      elseif(ls==3)
          gp = [0.1666666666, 0.1666666666;
              0.6666666666, 0.1666666666;
              0.1666666666, 0.6666666666];
          
          w_3 = [1/6 1/6 1/6];
              
          for i = 1:3
              psi = gp(i,1);
              eta = gp(i,2);
              
              N = [psi eta 1-psi-eta];
              
              GN    =  [1 0 -1 ;
                  0 1 -1];
              
              Jacobian     = GN*[xe ye];        % Compute Jacobian matrix 
              detJ  = det(Jacobian);
              
            inter(1) = inter(1) + w_3(i)*N*detJ*[1, 1, 1]';
            inter(2) = inter(2) + w_3(i)*N*detJ*xe;
            inter(3) = inter(3) + w_3(i)*N*detJ*ye;
            inter(4) = inter(4) + w_3(i)*xe'*(N'*N)*ye*detJ;
            inter(5) = inter(5) + w_3(i)*xe'*(N'*N)*xe*detJ;
            inter(6) = inter(6) + w_3(i)*ye'*(N'*N)*ye*detJ;
          
          end
          weight = A\inter;
      elseif(ls==6)
          gp = [ 0.445948490915965   0.445948490915965   
                 0.445948490915965   0.108103018168070   
                 0.108103018168070   0.445948490915965   
                 0.091576213509771   0.091576213509771  
                 0.091576213509771   0.816847572980459   
                 0.816847572980459   0.091576213509771];
          
          w_6 = [0.111690794839006 0.111690794839006 0.111690794839006 0.054975871827661 0.054975871827661 0.054975871827661]; 
              
          for i = 1:6
              psi = gp(i,1);
              eta = gp(i,2);
              
              N  = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
              
              GN = [4*psi-1    0     4*psi+4*eta-3  4*eta    -4*eta      4-8*psi-4*eta;
                        0    4*eta-1  4*eta+4*psi-3  4*psi 4-4*psi-8*eta    -4*psi];
              
              Jacobian     = GN*[xe ye];        % Compute Jacobian matrix 
              detJ  = det(Jacobian);
              
            inter(1) = inter(1) + w_6(i)*N*detJ*ones(6,1);
            inter(2) = inter(2) + w_6(i)*N*detJ*xe;
            inter(3) = inter(3) + w_6(i)*N*detJ*ye;
            inter(4) = inter(4) + w_6(i)*xe'*(N'*N)*ye*detJ;
            inter(5) = inter(5) + w_6(i)*xe'*(N'*N)*xe*detJ;
            inter(6) = inter(6) + w_6(i)*ye'*(N'*N)*ye*detJ;
          
          end
          weight = A\inter;
          
      end
end

