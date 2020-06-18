function [ke_diag,fe_gam] = robin_terms_bdy(elem,NODE,xint,yint,PARAMS)

ndof =  PARAMS.ndof;
nlink = PARAMS.nlink;

xe = [NODE(elem.nodes).X];
ye = [NODE(elem.nodes).Y];
%Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];

% Gauss quadrature points
%gauss = [-3^(-0.5), 3^(-0.5)];

gauss = [-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053];
 w_4 = [0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];

fe_gam = zeros(nlink*ndof,1);
ke_diag = zeros(nlink*ndof,nlink*ndof);

%%% Get the Jacobian
jcob_sub = 0.5*sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);
kepen = zeros(nlink*ndof,nlink*ndof);

for gp=1:4
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
    NparGamma = Npar;    
    
    %%%%%%%%%%%%%%%%%%%% Penalty Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kepen = kepen + w_4(gp)*jcob_sub*(NparGamma'*PARAMS.htCoeff*NparGamma);
    fe_gam = fe_gam + w_4(gp)*jcob_sub*NparGamma'*PARAMS.htCoeff*PARAMS.ambTemp;
end

ke_diag = ke_diag + kepen;
end