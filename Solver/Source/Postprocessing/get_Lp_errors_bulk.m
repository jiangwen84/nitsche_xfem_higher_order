function[L2Err,SigXXErr,H1Err] = get_Lp_errors_bulk(NODE,ELEM,VERT,VERT_ELE,PARAMS,u)
%% Remove Lagrange multipliers if necessary
disp = u(1:length(NODE)*PARAMS.ndof);
% [E,nu,~] = mat_param(PARAMS);

L2Ex = 0.0; L2Err = 0.0; SigXXErr = 0.0; SigXXEx = 0.0; H1Err = 0.0; H1Ex = 0.0; volEx = 0.0;
% NrgErr = 0.0; NrgEx = 0.0;
VOL = 0.0;

for e=1:length(ELEM)
    %% Set element quantities to zero
    volElm = 0.0;  L2ExElm = 0.0; L2ErrElm = 0.0; % NrgErrElm = 0.0; NrgExElm = 0.0;
    SigXXErrElm = 0.0; SigXXExElm = 0.0; H1ErrElm = 0.0; H1ExElm = 0.0;
    %% Get "interpolation" matrix and interpolated flux (constant in elem)
    %[nabUh] = elementflux(ELEM(e),NODE,PARAMS,disp);
    nabUh = [0 0 0]';
    [Amat] = interpolationmatrix(ELEM(e),NODE,PARAMS);
    Dmat = get_elas_tensor(ELEM(e),PARAMS);
    Sigh = Dmat*nabUh;
    %% Transform the physical part of elem into a cube
    [quadNod] = getDegenQuadCoord(ELEM(e),NODE,VERT, VERT_ELE, PARAMS);
    
    if(~isempty([quadNod(:,2).X]))
        totSub = 2;
    else
        totSub = 1;
    end
    
    for iSub=1:totSub
        if(length(quadNod(:,iSub))~=4)
            err = MException('degeCubChk:numberofNodeChk',...
                'for ELEM %d, only %d for degenerated cube',...
                e,length(quadNod(:,iSub)));
            throw(err);
        end
        
        %% Loop over the Gauss points in the 3 spatial directions
        for iPsi=1:4
            for iEta=1:4
                [jcob,Pos,Weight] = get_4thOrderGaussQuad_info(quadNod(:,iSub),iPsi,iEta);
                
                %% Get interpolated solution and exact solution and flux
                [uex,~,nabUex] = exactsolution(Pos.X,Pos.Y,ELEM(e).domain,PARAMS);
                [uh] = interpolatedsolution(Pos.X,Pos.Y,ELEM(e),Amat,PARAMS,disp);
                SigEx = Dmat*nabUex;
                %% Element volume
                volElm = volElm + Weight.psi*Weight.eta*jcob;
                
                %% L2 error
                L2ErrElm = L2ErrElm...
                    + Weight.psi*Weight.eta*jcob*(uh-uex)'*(uh-uex);
                L2ExElm = L2ExElm...
                    + Weight.psi*Weight.eta*jcob*(uex)'*(uex);
                
                %             %% Energy semi-norm
                %             NrgErrElm = NrgErrElm...
                %                 + Weight.psi*Weight.eta*jcob*(nabUh-nabUex)'*(nabUh-nabUex);
                %             NrgExElm = NrgExElm...
                %                 + Weight.psi*Weight.eta*jcob*(nabUex)'*(nabUex);
                
                %% Energy semi-norm
                SigXXErrElm = SigXXErrElm...
                    + Weight.psi*Weight.eta*jcob*(Sigh(1)-SigEx(1))'*(Sigh(1)-SigEx(1));
                SigXXExElm = SigXXExElm...
                    + Weight.psi*Weight.eta*jcob*(SigEx(1))'*(SigEx(1));
                
                %% H1 norm
                H1ErrElm = H1ErrElm...
                    + Weight.psi*Weight.eta*jcob*(...
                    (nabUh-nabUex)'*(nabUh-nabUex)+(uh-uex)'*(uh-uex));
                H1ExElm = H1ExElm...
                    + Weight.psi*Weight.eta*jcob*(...
                    (nabUex)'*(nabUex)+(uex)'*(uex));
            end
        end
    end
%     if(abs(volElm-ELEM(e).volume)>1.0e-11)
%         fprintf('ELEM(%d), abs(volElm-volEx) = %e\n',e,abs(volElm-ELEM(e).volume));
%     end
    [~,volElmEx] = convhulln([[quadNod.X]',[quadNod.Y]']);
    if(abs(volElm-volElmEx)>1.0e-11)
        fprintf('ELEM(%d), abs(volElm-volEx) = %e\n',e,abs(volElm-volElmEx));
    end
    
    L2Err = L2Err + L2ErrElm; L2Ex = L2Ex + L2ExElm; volEx = volEx + volElm;
    SigXXErr = SigXXErr + SigXXErrElm;
    SigXXEx = SigXXEx + SigXXExElm;
    H1Err = H1Err + H1ErrElm;
    H1Ex = H1Ex + H1ExElm;
    %     NrgErr = NrgErr + NrgErrElm;
    %     NrgEx = NrgEx + NrgExElm;
    
    %VOL = VOL + volElm;
end

L2Err = sqrt(L2Err/L2Ex)

SigXXErr = sqrt(SigXXErr/SigXXEx);
H1Err = sqrt(H1Err/H1Ex);
% NrgErr = sqrt(NrgErr/NrgEx);
end

function[nabUh] = elementflux(elem,NODE,PARAMS,disp)

ndof = PARAMS.ndof; nlink = PARAMS.nlink;

% get the B-matrix
if(nlink==4)
    [B] = sfderivatives(elem.nodes,[NODE.X],[NODE.Y],[NODE.Z],ndof);
elseif(nlink==3)
    [B] = sfderivatives(elem.nodes,NODE,ndof);
end

[ue] = get_nodal_solution(elem,PARAMS,disp);

% get the gradient of the primal variable
nabUh = B*ue;

end

function[uh] = interpolatedsolution(xpos,ypos,elem,Amat,PARAMS,disp)

% ndof = PARAMS.ndof; nlink = PARAMS.nlink;
%
% % get the tet shape functions
% Npar = Amat\[1 xpos ypos]';
%
% % get the nodal displacements ue and the interpolated displacement uh
% uh = zeros(ndof,1);
% for j=1:ndof
%     for i=1:nlink
%         uh(j,1) = uh(j,1) + Npar(i)*disp(ndof*(elem.nodes(i)-1)+j);
%     end
% end
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