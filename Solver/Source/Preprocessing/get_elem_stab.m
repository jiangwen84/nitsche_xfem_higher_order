function [ELEM] = get_elem_stab(ELEM,VERT,csize,locGrn,PARAMS)
Diag = eye(3); %Diag(3,3) = 2; % To account for the factor of 2 in shear strains
for i=1:length(locGrn)
    Dmat = get_elas_tensor(ELEM(csize+i),PARAMS);
    Vol = ELEM(csize+i).volume;
    for iSib=1:length(ELEM(csize+i).siblings)
        xVal = [VERT(ELEM(csize+i).vertices(iSib,1:2)).X];
        yVal = [VERT(ELEM(csize+i).vertices(iSib,1:2)).Y];
        length_seg = sqrt((xVal(2)-xVal(1))^2 + (yVal(2)-yVal(1))^2);
        
        sibId = ELEM(csize+i).siblings(iSib);
        DmatSib = get_elas_tensor(ELEM(sibId),PARAMS);
        sibVol = ELEM(sibId).volume;
        
        if(PARAMS.gamFlag==2)
            %%%%%%%%% Robust Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%%%%
            wei_vol = norm(DmatSib*Diag)*Vol + norm(Dmat*Diag)*sibVol;
            ELEM(csize+i).gamma(iSib) = Vol*norm(DmatSib*Diag)/wei_vol;
            D_avg = norm(Dmat*Diag)*norm(DmatSib*Diag)/wei_vol;
            ELEM(csize+i).stab(iSib) = 2*length_seg*D_avg;
        elseif(PARAMS.gamFlag==1)
            %%%%%%%%% Hansbo's formulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ELEM(csize+i).gamma(iSib) = Vol/(Vol + sibVol);
            D_avg = (norm(DmatSib*Diag)*sibVol + norm(Dmat*Diag)*Vol)/(Vol + sibVol)^2;
            ELEM(csize+i).stab(iSib) = 2*length_seg*D_avg;
        else
            %%%%%%%%% Classical Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%            
            ELEM(csize+i).gamma(iSib) = 0.5;
            D_avg = norm(DmatSib*Diag)/sibVol + norm(Dmat*Diag)/Vol;
            ELEM(csize+i).stab(iSib) = (2/4)*length_seg*D_avg;
        end
        %%%%% If element is also a part of the external boundary, multiply
        % by a factor 1.5 (refer notes on stabilization parameter).
        if(length(locGrn)>2)
            % If the child element does not have any background nodes then it
            % is possible that the element itself is not a part of external
            % boundary but the sibling is. However, the stabilization parameter
            % for the same segment can't have two values from 2 different
            % elements. This might cause violation of traction equilibrium. For
            % this reason, check whether both Element and sibling are part of
            % the external boundary.
            if(~isempty(ELEM(csize+i).boundX) ||~isempty(ELEM(sibId).boundX)...
                    || ~isempty(ELEM(csize+i).boundY) || ~isempty(ELEM(sibId).boundY))
                ELEM(csize+i).stab(iSib) = (3/2)*ELEM(csize+i).stab(iSib);
            end
        end
    end
end
