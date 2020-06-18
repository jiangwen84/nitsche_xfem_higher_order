function [keclass]  = stiffness(elem,NODE,ndof,nlink,weight,Dmat,gp)
keclass = zeros(ndof*nlink,ndof*nlink);
    
    for i = 1:6
        psi = gp(i,1);
        eta = gp(i,2);
        %%%%% Calculate B Matrix %%%%%%%%%%

        [B] = sfderivatives(elem.nodes,NODE,ndof,psi,eta);
        keclass = keclass + weight(i)*B'*Dmat*B;
     end

end
