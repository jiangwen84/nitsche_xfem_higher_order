function [fe_gam] = neumann_force_terms(elem,bdy_normal,NODE,xint,yint,PARAMS)

ndof =  PARAMS.ndof;
nlink = PARAMS.nlink;

xe = [NODE(elem.nodes).X];
ye = [NODE(elem.nodes).Y];
%Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];

% Gauss quadrature points
gauss = [-3^(-0.5), 3^(-0.5)];



%%% Get the Jacobian
jcob_sub = 0.5*sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);

fe_gam = zeros(nlink*ndof,1);

% Dmat = get_elas_tensor(elem,PARAMS);

PARAMS.bdy_norm = bdy_normal;

for gp=1:2
    tau = gauss(gp);
    
    %%% Get the parent shape functions
    N1_sub = 0.5*(1-tau); N2_sub = 0.5*(1+tau);

    xpos = N1_sub*xint(1) + N2_sub*xint(2);
    ypos = N1_sub*yint(1) + N2_sub*yint(2);
    
    [psi, eta] = inverse(xpos,ypos,xe,ye);
    
    if(length(xe)==3)
        Npar = [psi eta 1-psi-eta];
    elseif(length(xe)==6)
        Npar = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
    end

    %Npar = Amat\[1 xpos ypos]';
    
    
    
    %%% Get the parent shape function on the embedded interface
    if(ndof==2)
        if(length(xe)==3)
            NparGamma = [Npar(1),0,Npar(2),0,Npar(3),0;
                         0,Npar(1),0,Npar(2),0,Npar(3)];
        elseif(length(xe)==6)
            NparGamma = [Npar(1),0,Npar(2),0,Npar(3),0,Npar(4),0,Npar(5),0,Npar(6),0;
                         0,Npar(1),0,Npar(2),0,Npar(3),0,Npar(4),0,Npar(5),0,Npar(6)];
        end
    else
        NparGamma = Npar;
    end
        
    if(PARAMS.ndof==1)
        [~,FluxEx,~] = exactsolution(xpos,ypos,elem.domain,PARAMS);
%         FluxEx = Dmat*bdy_normal*nabUex;
    else
        [~,FluxEx,~] = exactsolution(xpos,ypos,elem.domain,PARAMS);
    end
    
    
   %FluxEx = FluxEx*0; %Yingjie
   FluxEx = [0;0]; %WJ
   
   fe_gam = fe_gam + jcob_sub*NparGamma'*FluxEx;
end

end