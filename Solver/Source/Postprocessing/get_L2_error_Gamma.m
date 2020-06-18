function [L2UjumpErr,L2NabUjumpErr] = get_L2_error_Gamma(NODE,ELEM,VERT,PARAMS,disp)

L2UjumpErr = 0; L2UjumpEx = 0;
L2NabUjumpErr = 0; L2NabUjumpEx = 0;

for e=1:length(ELEM)
    if(~isempty(ELEM(e).vertices))
        [L2UjumpErrElm,L2UjumpExElm,L2NabUjumpErrElm,L2NabUjumpExElm] =...
            get_error_elem_Gamma(ELEM(e),ELEM(ELEM(e).siblings(1)),NODE,...
            VERT,PARAMS,disp);
    else
        L2UjumpErrElm = 0;
        L2UjumpExElm = 0;
        L2NabUjumpErrElm = 0;
        L2NabUjumpExElm = 0;
    end
    L2UjumpErr = L2UjumpErr + L2UjumpErrElm;
    L2UjumpEx = L2UjumpEx + L2UjumpExElm;
    L2NabUjumpErr = L2NabUjumpErr + L2NabUjumpErrElm;
    L2NabUjumpEx = L2NabUjumpEx + L2NabUjumpExElm;
end

if(L2UjumpEx<1.0e-15)
    fprintf('Computing absolute error of flux on Gamma\n');
    L2UjumpErr = sqrt(L2UjumpErr)
else
    fprintf('Computing relative error of jump solution on Gamma\n');
    L2UjumpErr = sqrt(L2UjumpErr/L2UjumpEx)
end
if(L2NabUjumpEx<1.0e-15)
    fprintf('Computing absolute error of flux on Gamma\n');
    fprintf('L2NabUjumpErrElm = %7.4e\n',L2NabUjumpErr);
    L2NabUjumpErr = sqrt(L2NabUjumpErr)
else
    fprintf('Computing relative error of flux on Gamma\n');
    L2NabUjumpErr = sqrt(L2NabUjumpErr/L2NabUjumpEx)
end
end

function [L2UjumpErrElm,L2UjumpExElm,L2NabUErrElm,L2NabUExElm] =...
    get_error_elem_Gamma(elem,sibling,NODE,VERT,PARAMS,u)

%% Separate primal variable from Lagrange multipliers (if necessary)
disp = u(1:length(NODE)*PARAMS.ndof);
lambda = u(length(NODE)*PARAMS.ndof+1:end);

xe = [NODE(elem.nodes).X];
ye = [NODE(elem.nodes).Y];

mat_param;
ndof = PARAMS.ndof;

[ueElm] = get_nodal_solution(elem,PARAMS,disp);
[ueSib] = get_nodal_solution(sibling,PARAMS,disp);
[Amat] = interpolationmatrix(elem,NODE,PARAMS);

xint = [VERT(elem.vertices(1,1:2)).X];
yint = [VERT(elem.vertices(1,1:2)).Y];

if(ndof==2)
    P = [elem.normals(1,1)  elem.normals(1,2);
        elem.normals(1,2) -elem.normals(1,1)];
else
    P = 1;
end

if(ndof==1)
    normal_voigt = elem.normals;
elseif(ndof==2)
    %%%%%% Discretized stress tensor = normal_voigt*Dmat*B*d %%%%%%%
    normal_voigt = [elem.normals(1,1),0,elem.normals(1,2);
        0,elem.normals(1,2),elem.normals(1,1)];
end

L2UjumpErrElm = 0; L2UjumpExElm = 0;
L2NabUExElm = 0;L2NabUErrElm = 0;

jcobGam = 0.5*sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);

gauss = [-0.861136311594053, -0.339981043584856, 0.339981043584856,  0.861136311594053];
w_4 = [0.347854845137454 0.652145154862546 0.652145154862546 0.347854845137454];

for gp=1:4
    tau = gauss(gp);
    
    %tracXY_elm1 = get_B_normal(elem,P,elem.normals(1,1:2),PARAMS,elem.gamma(1));
    %tracXY_sib1 = get_B_normal(sibling,P,-elem.normals(1,1:2),PARAMS,sibling.gamma(1));
    
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
    
    uhElm = NparGamma*ueElm; uhSib = NparGamma*ueSib;
    
    %%%%%%%%%%%%%%%%%B Matrix%%%%%%%%%%%%%%%%%%%%%
    DmatElm = get_elas_tensor(elem,PARAMS);
    DmatSib = get_elas_tensor(sibling,PARAMS);
    [B] = sfderivatives(elem.nodes,NODE,ndof,psi,eta);
    
    %tracXY_elm = tracXY_elm1*B;
    %tracXY_sib = tracXY_sib1*B;
    
    tracXY_elm = DmatElm*B*ueElm;
    tracXY_sib = DmatSib*B*ueSib;
        
    if(PARAMS.nit_n == 1 || PARAMS.nit_t == 1)
        %avNabUh = -(elem.gamma*tracXY_elm*ueElm + sibling.gamma*tracXY_sib*ueSib);
        avNabUh = -(elem.gamma*normal_voigt*tracXY_elm + sibling.gamma*normal_voigt*tracXY_sib);
    else
        avNabUh = 0;
    end
    
    uhjump = (uhElm - uhSib);    
    uexjump = get_jump_interface(xpos,ypos,elem.domain,sibling.domain,PARAMS);
    
    if(PARAMS.nit_n==1 || PARAMS.nit_t==1)
        avNabUh = avNabUh + elem.stab*(uhjump-uexjump);
    elseif(isempty(lambda))
        avNabUh = avNabUh + PARAMS.pen_n*(uhjump-uexjump);
    end
    
    [~,~,nabUex_voigtElm] = exactsolution(xpos,ypos,elem.domain,PARAMS);
    [~,~,nabUex_voigtSib] = exactsolution(xpos,ypos,sibling.domain,PARAMS);

    nabUexElm = -normal_voigt*DmatElm*nabUex_voigtElm;
    nabUexSib = (-normal_voigt)*DmatSib*nabUex_voigtSib;
    avNabUex = elem.gamma*nabUexElm + sibling.gamma*nabUexSib;
        
    L2UjumpErrElm = L2UjumpErrElm...
        + jcobGam*w_4(gp)*(uhjump-uexjump)'*(uhjump-uexjump);
    L2UjumpExElm = L2UjumpExElm + jcobGam*w_4(gp)*(uexjump'*uexjump);
    L2NabUErrElm = L2NabUErrElm...
        + jcobGam*w_4(gp)*(avNabUh - avNabUex)'*(avNabUh - avNabUex);
    L2NabUExElm = L2NabUExElm + jcobGam*w_4(gp)*(avNabUex'*avNabUex);
end
end