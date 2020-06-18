function [u] = solve_nitsche(ELEM,NODE,VERT,VERT_ELE,PARAMS)

fprintf('*************************************************************************** \n')
if(PARAMS.nit_t||PARAMS.nit_n)
    fprintf('                         Nitsche''s Method')
    if (PARAMS.gamFlag==2)
        fprintf(' : Robust \n');
    elseif (PARAMS.gamFlag==1)
        fprintf(' : Hansbo \n');
    else
        fprintf(' : Classical \n');
    end
else
    fprintf('                         Penalty Method \n')
end
fprintf('*************************************************************************** \n')

%%%%%%%% Prescribe boundary conditions %%%%%%%%%%%
[ifixu, u] = essentialbcs(NODE,PARAMS,zeros(PARAMS.ndof*length(NODE),1));
%%%%%%% Calculate global stiffness and external force %%%%%%%%%%%%
fprintf('Assembling linear algebraic system \n')
[bigk,fext] = assemble_nitsche(NODE,ELEM,VERT,VERT_ELE,PARAMS);
%%%%%%%% Apply boundary conditions %%%%%%%%%%%%%
[bigk,fext] = applybcs(ifixu,u,bigk,fext);

u = bigk\fext;
end