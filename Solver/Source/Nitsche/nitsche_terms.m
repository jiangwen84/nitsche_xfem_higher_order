function [ke_diag,ke_offdiag,fe_gam] = nitsche_terms(elem,siblings,NODE,VERT,PARAMS)

% mat_param;

ndof =  PARAMS.ndof;
nlink = PARAMS.nlink;
nit_n = PARAMS.nit_n; pen_n = PARAMS.pen_n;
nit_t = PARAMS.nit_t; pen_t = PARAMS.pen_t;

xe = [NODE(elem.nodes).X];
ye = [NODE(elem.nodes).Y];
%Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];

% Gauss quadrature points
%gauss = [-3^(-0.5), 3^(-0.5)];

gauss = [-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053];
w_4 = [0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];


fe_gam = zeros(nlink*ndof,1);
ke_diag = zeros(nlink*ndof,nlink*ndof);
ke_offdiag = zeros(nlink*ndof,nlink*ndof,2);

for iSib=1:length(siblings)
    xint = [VERT(elem.vertices(iSib,1:2)).X];
    yint = [VERT(elem.vertices(iSib,1:2)).Y];
    
    %%% Get the penalty or stabilization parameter(s)
    if(ndof==2)
        %%%% Transformation matrix from (x,y) to (n,tau)
        %%%% alpha_(x,y) = P'*alpha_(n,tau)*P
        P = [elem.normals(iSib,1)  elem.normals(iSib,2);
            elem.normals(iSib,2) -elem.normals(iSib,1)];
        if(nit_n), alpha_n = elem.stab(iSib); else alpha_n = pen_n; end
        if(nit_t), alpha_t = elem.stab(iSib); else alpha_t = pen_t; end
        alpha_e = P'*[alpha_n 0; 0 alpha_t]*P;
        %alpha_e = [100 0; 0 100];
    else
        P = 1;
        if(nit_n)
            alpha_e = elem.stab(iSib);
        else
            alpha_e = pen_n; 
        end
    end
    
    %%% Get the tractions
    tracXY_elm1 = get_B_normal(elem,P,elem.normals(iSib,1:2),PARAMS,elem.gamma(iSib));
    [gamId] = get_sibling_weight_id(elem,siblings,iSib,PARAMS);
    tracXY_sib1 = get_B_normal(siblings(iSib),P,-elem.normals(iSib,1:2),PARAMS,siblings(iSib).gamma(gamId));
     
    %%% Get the Jacobian
    jcob_sub = 0.5*sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);
    
    kepen = zeros(nlink*ndof,nlink*ndof);
    kenit_elm = zeros(nlink*ndof,nlink*ndof);
    kenit_sib = zeros(nlink*ndof,nlink*ndof);
    fe_gamiSib = zeros(nlink*ndof,1);
    
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
         
       % Npar = Amat\[1 xpos ypos]';
        
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
        
        
        %%%%%%%%%%%%%%%%%B Matrix%%%%%%%%%%%%%%%%%%%%%
        [B] = sfderivatives(elem.nodes,NODE,ndof,psi,eta);
        
        
        tracXY_elm = tracXY_elm1*B;
        tracXY_sib = tracXY_sib1*B;
        
        if(PARAMS.ndof==1 && PARAMS.jump==1)
            [u_jump] = get_jump_interface(xpos,ypos,elem.domain,siblings(iSib).domain,PARAMS);
        end
        
        if(PARAMS.tjump == 1)
            [t_jump] = get_jump_flux(xpos,ypos,elem,siblings(iSib),PARAMS);
        end
        
        %t_jump = -2;
        %%%%%%%%%%%%%%%%%%%% Nitsche Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kenit_elm = kenit_elm - w_4(gp)*jcob_sub*(NparGamma'*tracXY_elm);
        kenit_sib = kenit_sib - w_4(gp)*jcob_sub*(NparGamma'*tracXY_sib);
        if(PARAMS.ndof==1 && PARAMS.jump==1)
            fe_gamiSib = fe_gamiSib - w_4(gp)*jcob_sub*u_jump*tracXY_elm';
        end
        
         if(PARAMS.tjump == 1)
            fe_gamiSib = fe_gamiSib + w_4(gp)*jcob_sub*siblings(iSib).gamma(gamId)*NparGamma'*t_jump;
        end
        
        %%%%%%%%%%%%%%%%%%%% Penalty Terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kepen = kepen + w_4(gp)*jcob_sub*(NparGamma'*alpha_e*NparGamma);
        if(PARAMS.ndof==1 && PARAMS.jump==1)
            fe_gamiSib = fe_gamiSib + w_4(gp)*jcob_sub*NparGamma'*alpha_e*u_jump;
        end
    end

    %%% Symmetry terms contribution to tangent stiffness
    ke_offdiag(:,:,iSib) = -kepen - (kenit_sib + kenit_elm');
    ke_diag = ke_diag + kepen + (kenit_elm + kenit_elm');
    fe_gam = fe_gam+fe_gamiSib;
end

end

function [gamId] = get_sibling_weight_id(elem,siblings,iSib,PARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each element with a junction, interface is piecewise linear.
% For each piecewise linear segment, weights have to sum to unity.
% Here, we find the corresponding segment weight to elem.gamma(iSib) in siblings(iSib).gamma
gamId = find(siblings(iSib).siblings==PARAMS.e);
sumOfWeights = elem.gamma(iSib) + siblings(iSib).gamma(gamId);
%%%%%% Sum of weights on each segment should be unity for variational
%%%%%% consistency. Throw error otherwise.
if(abs(sumOfWeights-1)>1e-15)
    err = MException('ConsistencyChk:Variationally inconsistent',...
        'Weights do not sum to one in element %d for sibling %d',PARAMS.e, elem.siblings(iSib));
    throw(err);
end
%%%%%% Additional check: normals. In a similar fashion, for each piecewise
%%%%%% linear segment, normals have to be equal in magnitude and
%%%%%% opposing in sign between the element and the sibling.
diffNormal = -elem.normals(iSib,1:2)-siblings(iSib).normals(gamId,1:2);
if(norm(diffNormal)>1e-15)
    err = MException('NormalSignChk:WrongSign',...
        'Wrong sign for segment normal in element %d for sibling %d',PARAMS.e, elem.siblings(iSib));
    throw(err);
end
end

