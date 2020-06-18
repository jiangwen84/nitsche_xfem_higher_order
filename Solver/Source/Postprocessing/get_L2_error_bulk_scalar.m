function[L2Err] = get_L2_error_bulk_scalar(NODE,ELEM,VERT,VERT_ELE,PARAMS,u)
%% Remove Lagrange multipliers if necessary
disp = u(1:length(NODE)*PARAMS.ndof); mat_param;

L2Ex = 0.0; L2Err = 0.0;NrgErr = 0.0; NrgEx = 0.0;

VOL = 0.0; %% debug

for e=1:length(ELEM)
    
    %% Get "interpolation" matrix and interpolated flux (constant in elem)    
    [Amat] = interpolationmatrix(ELEM(e),NODE,PARAMS);   
    %% Transform the physical part of elem into a cube
    [quadNod] = getDegenQuadCoord(ELEM(e),NODE,VERT, VERT_ELE, PARAMS);
    if(length(quadNod)~=4)
        err = MException('degeCubChk:numberofNodeChk',...
            'for ELEM %d, only %d for degenerated cube',...
            e,length(quadNod));
        throw(err);
    end
    
    %% Set element quantities to zero
    volElm = 0.0;  L2ExElm = 0.0; L2ErrElm = 0.0; NrgErrElm = 0.0; NrgExElm = 0.0;
    %% Loop over the Gauss points in the 3 spatial directions
    for iPsi=1:4
        for iEta=1:4
            [jcob,Pos,Weight] = get_4thOrderGaussQuad_info(quadNod,iPsi,iEta);
            
            [nabUh] = elementflux_new(ELEM(e),NODE,PARAMS,disp,Pos.X,Pos.Y,Amat);
            
            %% Get interpolated solution and exact solution and flux
            [uex,~,nabUex] = exactsolution(Pos.X,Pos.Y,ELEM(e).domain,PARAMS);
            [uh] = interpolatedsolution(Pos.X,Pos.Y,ELEM(e),Amat,PARAMS,disp);
            %% Element volume
            volElm = volElm + Weight.psi*Weight.eta*jcob;
            
            %% L2 error
            L2ErrElm = L2ErrElm...
                + Weight.psi*Weight.eta*jcob*(uh-uex)'*(uh-uex);
            L2ExElm = L2ExElm...
                + Weight.psi*Weight.eta*jcob*(uex)'*(uex);
            
            %             %% Energy semi-norm
            NrgErrElm = NrgErrElm...
                + Weight.psi*Weight.eta*jcob*(nabUh-nabUex)'*(nabUh-nabUex);
            NrgExElm = NrgExElm...
                + Weight.psi*Weight.eta*jcob*(nabUex)'*(nabUex);
        end
    end
    
    [~,volElmEx] = convhulln([[quadNod.X]',[quadNod.Y]']);
    if(abs(volElm-volElmEx)>1.0e-11)
        fprintf('ELEM(%d), abs(volElm-volEx) = %e\n',e,abs(volElm-volElmEx));
    end
    L2Err = L2Err + L2ErrElm; L2Ex = L2Ex + L2ExElm;
    NrgErr = NrgErr + NrgErrElm;
    NrgEx = NrgEx + NrgExElm;
    VOL = VOL + volElm;
end

L2Err = sqrt(L2Err/L2Ex)
NrgErr = sqrt(NrgErr/NrgEx)
%L2Err = sqrt(L2Err);
end


function[uh] = interpolatedsolution(xpos,ypos,elem,Amat,PARAMS,disp)

ndof = PARAMS.ndof; nlink = PARAMS.nlink;


xe = Amat(2,:);
ye = Amat(3,:);

[psi, eta] = inverse(xpos,ypos,xe,ye);

% get the tet shape functions
%Npar = Amat\[1 xpos ypos]';

 if(length(xe)==3)
     Npar = [psi eta 1-psi-eta];
 elseif(length(xe)==6)
     Npar = [2*psi^2-psi 2*eta^2-eta 2*psi^2+2*eta^2+4*psi*eta-3*psi-3*eta+1 4*psi*eta 4*eta*(1-psi-eta) 4*psi*(1-psi-eta)];
 end

% get the nodal displacements ue and the interpolated displacement uh
uh = zeros(ndof,1);
for j=1:ndof
    for i=1:nlink
        uh(j,1) = uh(j,1) + Npar(i)*disp(ndof*(elem.nodes(i)-1)+j);
    end
end
end


function[nabUh] = elementflux_new(elem,NODE,PARAMS,disp,xpos,ypos,Amat)

ndof = PARAMS.ndof; nlink = PARAMS.nlink;

xe = Amat(2,:);
ye = Amat(3,:);

[psi, eta] = inverse(xpos,ypos,xe,ye);

% get the B-matrix
if(nlink==4)
    [B] = sfderivatives(elem.nodes,[NODE.X],[NODE.Y],[NODE.Z],ndof);
elseif(nlink==3 || nlink==6)
    [B] = sfderivatives(elem.nodes,NODE,ndof,psi,eta);
end

[ue] = get_nodal_solution(elem,PARAMS,disp);

% get the gradient of the primal variable
nabUh = B*ue;

end