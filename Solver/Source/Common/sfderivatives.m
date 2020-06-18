function [B] = sfderivatives(loc_conn,NODE,ndof,psi,eta)

%%%%%%%%% Nodal coordinates for the element %%%%%%%%%%
x = [NODE(loc_conn).X];
y = [NODE(loc_conn).Y];

C = [x' y'];

if(length(x)==3)
    
    GN = [1  0  -1;
          0  1  -1];
    Jacob     = GN*C;       % compute Jacobian matrix 
%   detJ  = det(Jacob);     % Compute Jacobian matrix 
      
    if (ndof == 2)   
        BB     = Jacob\GN;       % compute the derivative of the shape functions
        
        B1x     = BB(1,1);
        B2x     = BB(1,2);
        B3x     = BB(1,3);
        
        B1y     = BB(2,1);
        B2y     = BB(2,2);
        B3y     = BB(2,3); 
         
        B = [ B1x      0     B2x     0      B3x    0;
                0     B1y     0     B2y      0     B3y; 
              B1y     B1x    B2y    B2x     B3y    B3x]; 
    else
        B = Jacob\GN;
    end

elseif(length(x)==6)
     GN = [4*psi-1    0     4*psi+4*eta-3  4*eta    -4*eta      4-8*psi-4*eta;
              0    4*eta-1  4*eta+4*psi-3  4*psi 4-4*psi-8*eta    -4*psi];
          
      Jacob = GN*C;
      %det = det(Jacob);
      if(ndof==2)
          BB = Jacob\GN;
          
          B1x = BB(1,1); B2x = BB(1,2);
          B3x = BB(1,3); B4x = BB(1,4);
          B5x = BB(1,5); B6x = BB(1,6);    
        
          B1y = BB(2,1); B2y = BB(2,2);
          B3y = BB(2,3); B4y = BB(2,4);
          B5y = BB(2,5); B6y = BB(2,6);
         
        B = [ B1x      0     B2x     0      B3x    0    B4x 0   B5x   0   B6x 0
                0     B1y     0     B2y      0     B3y   0  B4y  0   B5y   0  B6y
              B1y     B1x    B2y    B2x     B3y    B3x  B4y B4x  B5y B5x  B6y B6x]; 
      else
          B = Jacob\GN;
      end
end
      
    
end
    