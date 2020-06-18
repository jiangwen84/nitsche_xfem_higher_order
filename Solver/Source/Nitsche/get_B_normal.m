%%%%% Function to get the gradNdotn terms arising in a Nitsche formulation
function [gradXY,gradNT] = get_B_normal(elem,P,normal,PARAMS,gamma)

Dmat = get_elas_tensor(elem,PARAMS);
if(PARAMS.ndof==1)
    normal_voigt = normal;
elseif(PARAMS.ndof==2)
    %%%%%% Discretized stress tensor = normal_voigt*Dmat*B*d %%%%%%%
    normal_voigt = [normal(1),0,normal(2);
        0,normal(2),normal(1)];
end
gradXY = gamma*normal_voigt*Dmat;
gradNT = P*gradXY;

if(~PARAMS.nit_n), gradNT(1,:) = zeros(1,PARAMS.nlink*PARAMS.ndof); end
if(~PARAMS.nit_t && PARAMS.ndof==2), gradNT(2,:) = zeros(1,PARAMS.nlink*PARAMS.ndof); end
gradXY = P'*gradNT;

end