function [stress,vonMises] = getStress(ELEM,NODE,PARAMS,u)

if(PARAMS.ndof==2)
    stress = zeros(length(ELEM),3); vonMises = zeros(length(ELEM),1);
else
    stress = zeros(length(ELEM),2); vonMises = zeros(length(ELEM),1);
end

gauss_six = [0.445948490915965   0.445948490915965   
             0.445948490915965   0.108103018168070   
             0.108103018168070   0.445948490915965   
             0.091576213509771   0.091576213509771  
             0.091576213509771   0.816847572980459   
             0.816847572980459   0.091576213509771];

for e=1:length(ELEM)
    
    psi = sum(gauss_six(:,1))/6;
    eta = sum(gauss_six(:,2))/6;
    %%%%% Calculate B Matrix %%%%%%%%%%
    [B] = sfderivatives(ELEM(e).nodes,NODE,PARAMS.ndof,psi,eta);
    
    Dmat = get_elas_tensor(ELEM(e),PARAMS);

    [uLoc] = get_local_solution(ELEM(e),PARAMS,u);
    stress(e,:)=Dmat*B*uLoc;
    if(PARAMS.ndof==2)
        vonMises(e) = sqrt(0.5*((stress(e,1)-stress(e,2))^2 + stress(e,1)^2 +...
            stress(e,2)^2) + 3*stress(e,3)^2);
    end
end

end