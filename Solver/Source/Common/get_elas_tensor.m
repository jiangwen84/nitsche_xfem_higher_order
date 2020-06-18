function[Dmat] = get_elas_tensor(elem,PARAMS)

if(PARAMS.ndof==1)
    Dmat = elem.MatParam(1);
elseif(PARAMS.ndof==2)
    E = elem.MatParam(1); nu = elem.MatParam(2);
    if (PARAMS.plane_strain)
        %%%%% Lame parameters for the two subdomains %%%%%%%%%%%%%%%%%%%%%%
        lambda = nu*E/((1+nu)*(1-2*nu));
        mu = E/(2*(1+nu));
        
        %%%%% D matrix for the Elasticity tensor %%%%%%%%%%%%%%%%%%%%%%%%%%
        Dmat = [(lambda+2*mu),        lambda,  0;
            lambda, (lambda+2*mu),  0;
            0,             0, mu];
    else
        %%%%% D matrix for the Elasticity tensor %%%%%%%%%%%%%%%%%%%%%%%%%%
        fac = E/(1. - nu^2);
        Dmat = fac*[  1, nu,        0;
            nu,  1,        0;
            0,  0, (1-nu)/2];
    end
end

end