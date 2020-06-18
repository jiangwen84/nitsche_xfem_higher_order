function[fe] = get_body_force(elem,NODE,VERT,PARAMS)

%% Get "interpolation" matrix
[Amat] = interpolationmatrix(elem,NODE,PARAMS);

%% Transform the physical part of elem into a cube
[quadNod] = getDegenQuadCoord(elem,NODE,VERT);
if(length(quadNod)~=4)
    err = MException('degeCubChk:numberofNodeChk',...
        'for ELEM %d, only %d for degenerated cube',...
        e,length(cubNod));
    throw(err);
end
if(~isempty([quadNod(:,2).X]))
    totSub = 2;
else
    totSub = 1;
end
%% Set element quantities to zero
volElm = 0.0; fe = zeros(PARAMS.nlink*PARAMS.ndof,1);
for iSub=1:totSub
    
    %% Loop over the Gauss points in the 3 spatial directions
    for iPsi=1:4
        for iEta=1:4
            
            [jcob,Pos,Weight] = get_4thOrderGaussQuad_info(quadNod(:,iSub),iPsi,iEta);
            
            % get the tet shape functions
            Npar = Amat\[1 Pos.X Pos.Y]';
            
            [fBar] = bodyforce(Pos.X,Pos.Y,Pos.Z,elem.domain,PARAMS);
            
            %% Element volume
            volElm = volElm + Weight.psi*Weight.eta*jcob;
            
            fe = fe + Weight.psi*Weight.eta*jcob*fBar*Npar;
            
        end
    end
end
if(abs(volElm-elem.volume)>1.0e-15)
    fprintf('ELEM(%d), abs(volElm-volEx) = %e\n',e,abs(volElm-volElmEx));
end

end