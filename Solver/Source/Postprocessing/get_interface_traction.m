function [int_traction,CentCoord] = get_interface_traction(NODE,ELEM,VERT,PARAMS,disp)

mat_param;
int_traction = zeros(length([ELEM.siblings]),1); CentCoord = zeros(length([ELEM.siblings]),1);
cutEleCount = 0;

for e=1:length(ELEM)
    if(~isempty(ELEM(e).vertices) && ELEM(e).domain==1)        
        [B] = sfderivatives(ELEM(e).nodes,NODE,PARAMS.ndof);
        [tracElm] = get_cut_element_traction(ELEM(e),ELEM(ELEM(e).siblings),...
            [NODE.X],[NODE.Y],B,VERT,PARAMS,disp);
        centEleCoord(1) = 0.5*(VERT(ELEM(e).vertices(1)).X + VERT(ELEM(e).vertices(2)).X);
        centEleCoord(2) = 0.5*(VERT(ELEM(e).vertices(1)).Y + VERT(ELEM(e).vertices(2)).Y);
        if(PARAMS.ndof==2)
            int_traction(2*cutEleCount+1:2*cutEleCount+2) = tracElm;
            CentCoord(2*cutEleCount+1:2*cutEleCount+2) = centEleCoord;
        else
            int_traction(cutEleCount+1) = tracElm;
        end
        cutEleCount = cutEleCount+1;
    end
end
end

function [tracElm] = get_cut_element_traction(elem,sibling,x,y,B,VERT,PARAMS,disp)

xe = x(elem.nodes); ye = y(elem.nodes);

mat_param;

nlink = PARAMS.nlink; ndof = PARAMS.ndof; tracElm = zeros(ndof,1);

%%%%%%%% Compute the stabilized or penalty part %%%%%%%%%%%%%%%%%%%%%%%%%%%
ueElm = get_local_solution(elem,PARAMS,disp); 
ueSib = get_local_solution(sibling,PARAMS,disp);

Amat = [ 1 1 1; xe(1) xe(2) xe(3); ye(1) ye(2) ye(3)];

%2-pt gauss rule
gauss = [-3^(-0.5), 3^(-0.5)];

xint = [VERT(elem.vertices).X];
yint = [VERT(elem.vertices).Y];

eta = [gauss(1) gauss(2)];

if(PARAMS.nit_n==1 || PARAMS.nit_t==1)
    alpha_e = elem.stab;
else
    alpha_e = PARAMS.pen_n;
end

for quadPoint=1:length(eta)
    NparGamma = get_shape_functions_gamma(xint,yint,eta(quadPoint),Amat,ndof);
    uhElm = NparGamma*ueElm; uhSib = NparGamma*ueSib;
    if(elem.domain==1)
        uhjump = (uhSib - uhElm);
    else
        uhjump = (uhElm - uhSib);
    end
    tracElm = tracElm + (1/length(eta))*alpha_e*uhjump;
end

%%%%%%%% Compute the Nitsche's part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(PARAMS.nit_n==1 || PARAMS.nit_t==1)
    if(PARAMS.ndof==2)
        DmatElm = get_elas_tensor(elem,PARAMS);
        DmatSib = get_elas_tensor(sibling,PARAMS);
        
        normal_voigt = [elem.normals(1),0,elem.normals(2);
            0,elem.normals(2),elem.normals(1)];
        sigmaAvg = elem.gamma*DmatElm*B*ueElm + sibling.gamma*DmatSib*B*ueSib;
        if(elem.domain==1)
            tracElm = tracElm +  normal_voigt*sigmaAvg;
        else
            tracElm = tracElm -  normal_voigt*sigmaAvg;
        end
    else
        fluxAvg = elem.gamma*E(elem.domain)*B*ueElm + sibling.gamma*E(sibling.domain)*B*ueSib;
        if(elem.domain==1)
            tracElm = tracElm + elem.normals*fluxAvg;
        else
            tracElm = tracElm - elem.normals*fluxAvg;
        end
    end
end
end
function [NparGamma] = get_shape_functions_gamma(xint,yint,eta,Amat,ndof)

Na_sub = 0.5*(1-eta); Nb_sub = 0.5*(1+eta);
xpos = Na_sub*xint(1) + Nb_sub*xint(2);
ypos = Na_sub*yint(1) + Nb_sub*yint(2);

Npar = Amat\[1 xpos ypos]';
if(ndof==2)
    NparGamma = [Npar(1),0,Npar(2),0,Npar(3),0;
        0,Npar(1),0,Npar(2),0,Npar(3)];
else
    NparGamma = Npar';
end
end