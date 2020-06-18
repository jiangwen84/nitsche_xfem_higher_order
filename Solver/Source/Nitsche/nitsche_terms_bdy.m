function [ke_diag,fe_gam] = nitsche_terms_bdy(elem,bdy_normal,NODE,xint,yint,PARAMS)

ndof =  PARAMS.ndof;
nlink = PARAMS.nlink;

xe = [NODE(elem.nodes).X];
ye = [NODE(elem.nodes).Y];
%Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];

% Gauss quadrature points
%gauss = [-3^(-0.5), 3^(-0.5)];

fe_gam = zeros(nlink*ndof,1);
ke_diag = zeros(nlink*ndof,nlink*ndof);

alpha_n = 0; alpha_t = 0;

gauss = [-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053];
w_4 = [0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];

%%% Get the penalty or stabilization parameter(s)
if(ndof==2)
    %%%% Transformation matrix from (x,y) to (n,tau)
    %%%% alpha_(x,y) = P'*alpha_(n,tau)*P
    P = [bdy_normal(1,1)  bdy_normal(1,2);
        bdy_normal(1,2) -bdy_normal(1,1)];

    if(PARAMS.dirBdy_n)
        if(PARAMS.nit_bdy)
            if(bdy_normal(1,1)==0)
                alpha_n = elem.stab_dirich(2);
            elseif(bdy_normal(1,2)==0)
                alpha_n = elem.stab_dirich(1);
            end
        else
            alpha_n = PARAMS.pen_bdy;
        end
    end
    
    if(PARAMS.dirBdy_t)
        if(PARAMS.nit_bdy)
            if(bdy_normal(1,1)==0)
                alpha_t = elem.stab_dirich(2);
            elseif(bdy_normal(1,2)==0)
                alpha_t = elem.stab_dirich(1);
            end
        else
            alpha_t = PARAMS.pen_bdy;
        end
    end
      
    alpha_e = P'*[alpha_n 0; 0 alpha_t]*P;
    %alpha_e = [100 0; 0 100];
else
    P = 1;
    if(PARAMS.nit_bdy)
        if(bdy_normal(1,1)==0)
            alpha_e = elem.stab_dirich(2);
        elseif(bdy_normal(1,2)==0)
            alpha_e = elem.stab_dirich(1);
        end
    else
        alpha_e = PARAMS.pen_bdy;
    end
end

%%% Get the tractions
tracXY_elm1 = get_traction(elem,P,bdy_normal(1:2),PARAMS);

%%% Get the Jacobian
jcob_sub = 0.5*sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);

kepen = zeros(nlink*ndof,nlink*ndof);
kenit_elm = zeros(nlink*ndof,nlink*ndof);

for gp=1:4
    tau = gauss(gp);
    
    %%% Get the parent shape functions
    N1_sub = 0.5*(1-tau); N2_sub = 0.5*(1+tau);

    xpos = N1_sub*xint(1) + N2_sub*xint(2);
    ypos = N1_sub*yint(1) + N2_sub*yint(2);

    [psi, eta] = inverse(xpos,ypos,xe,ye);
    
    if(length(xe) == 3)
        Npar = [psi eta 1-psi-eta];
    elseif(length(xe)==6)
        Npar = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
    end
    
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
    
    
    [B] = sfderivatives(elem.nodes,NODE,ndof,psi,eta);
    tracXY_elm = tracXY_elm1*B;
    
    [uex] = exactsolution(xpos,ypos,elem.domain,PARAMS);
    
    %%%%%%%%%%%%%%%%%%%% Nitsche Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(PARAMS.nit_bdy)
        kenit_elm = kenit_elm - w_4(gp)*jcob_sub*(NparGamma'*tracXY_elm);
        fe_gam = fe_gam - w_4(gp)*jcob_sub*tracXY_elm'*uex;
    end
    
    %%%%%%%%%%%%%%%%%%%% Penalty Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kepen = kepen + w_4(gp)*jcob_sub*(NparGamma'*alpha_e*NparGamma);
    fe_gam = fe_gam + w_4(gp)*jcob_sub*NparGamma'*alpha_e*uex;
end

ke_diag = ke_diag + kepen + kenit_elm + kenit_elm';
end

function [tracXY] = get_traction(elem,P,normal,PARAMS)

Dmat = get_elas_tensor(elem,PARAMS);
if(PARAMS.ndof==1)
    normal_voigt = normal;
elseif(PARAMS.ndof==2)
    %%%%%% Discretized stress tensor = normal_voigt*Dmat*B*d %%%%%%%
    normal_voigt = [normal(1),0,normal(2);
        0,normal(2),normal(1)];
end
tracNT = P*normal_voigt*Dmat;
if(~PARAMS.dirBdy_n), tracNT(1,:) = zeros(1,PARAMS.nlink*PARAMS.ndof); end
if(~PARAMS.dirBdy_t && PARAMS.ndof==2), tracNT(2,:) = zeros(1,PARAMS.nlink*PARAMS.ndof); end
tracXY = P'*tracNT;

end