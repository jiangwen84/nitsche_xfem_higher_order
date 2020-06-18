function [bigk,fext] = assemble_nitsche(NODE,ELEM,VERT,VERT_ELE,PARAMS)

num_bg_ele=PARAMS.num_bg_ele; % original fem mesh!
ndof = PARAMS.ndof; nlink=PARAMS.nlink;

elem_mat_size = (nlink*ndof)^2;
num_dof = ndof*length(NODE);
numele = length(ELEM);
n_cut_ele = PARAMS.numele-num_bg_ele;  % orginal fem mesh!
fext = zeros(num_dof,1);
% For a cut elem, we have 2 diag and 2 off-diag contributions to the stiffness, each of size elem_mat_size.
alloc_size = 4*elem_mat_size*n_cut_ele + elem_mat_size*(num_bg_ele-n_cut_ele);
I = zeros(alloc_size,1); J = zeros(alloc_size,1); V = zeros(alloc_size,1);
next = 0; 

gauss_six = [0.445948490915965   0.445948490915965   
             0.445948490915965   0.108103018168070   
             0.108103018168070   0.445948490915965   
             0.091576213509771   0.091576213509771  
             0.091576213509771   0.816847572980459   
             0.816847572980459   0.091576213509771];



for e=1:numele
    
    PARAMS.e = e;    
    
    Dmat = get_elas_tensor(ELEM(e),PARAMS);
    
    [weight] = get_weight(ELEM(e),NODE,VERT,VERT_ELE,PARAMS);
    
    %%%%%%%%%%Classic Stiffness Matrix%%%%%%%%%%%%
    
    [keclass]  = stiffness(ELEM(e),NODE,ndof,nlink,weight,Dmat,gauss_six);
    
                
    %%%%%%% Get body force %%%%%%%
    if(ndof==1)
       fe = get_bodyforce(ELEM(e),NODE,ndof,nlink,PARAMS,weight,gauss_six);
    else
        fe = zeros(ndof*nlink,1);
    end    
    
    if(~isempty(ELEM(e).vertices))          %%%% Check if element is cut
        %%%%%%%%%%%%% Get Nitsche and penalty terms for stiffness matrix %%
        [ke_diag,ke_offdiag,fe_gam] = nitsche_terms(ELEM(e),ELEM(ELEM(e).siblings),NODE,VERT,PARAMS);
    else
        ke_diag = zeros(nlink*ndof,nlink*ndof); ke_offdiag = zeros(nlink*ndof,nlink*ndof);
        fe_gam = zeros(nlink*ndof,1);
    end
    
    %%%%% Apply Dirichlet and Neumann conditions on the outer boundary
    %%%%% Handle Dirichlet conditions weakly using Nitsche's method
     if(PARAMS.weakBC)
        [ke_bdy,fe_bdy] = apply_boundary_conditions(ELEM(e),NODE,PARAMS);
    else
        ke_bdy = zeros(nlink*ndof,nlink*ndof);
        fe_bdy = zeros(nlink*ndof,1);
    end

    %%%%%%% (Block) Diagonal terms of global stiffness %%%%%%%%%%%%%
    %%%%%%% B matrix is a constant, therefore ke = B'*B*volume.
   % keclass = B'*Dmat*B*ELEM(e).volume;
   
   
    ke = keclass + ke_diag + ke_bdy;
        
    %%%%%%%%%% Assembly procedure %%%%%%%%%%%%%%
    for inod=1:nlink
        rg = (ndof*(ELEM(e).nodes(inod)-1) + 1);
        re = ndof*(inod-1)+1;
        %%%%%%% Assemble a row of the global stiffness %%%%%%%%%%%%
        for jnod=1:nlink
            %%%%%%% Get column indices for the block diagonal entries and
            %%%%%%% row indices for both diagonal and off-diagonal entries
            cbk_diag = ndof*(ELEM(e).nodes(jnod)-1) + 1;                       
            ce = ndof*(jnod-1)+1;
            for idof=1:ndof
                for jdof=1:ndof
                    I(next+1) = rg + (idof-1); J(next+1) = cbk_diag + (jdof-1);
                    V(next+1) = ke(re+(idof-1),ce+(jdof-1));
                    next = next+1;
                end
            end
            %%%%%%% For uncut elements there is no off-diagonal assembly %%
            if(~isempty(ELEM(e).siblings))
                for iSib=1:length(ELEM(e).siblings)
                    %%%%%%% Assemble cut elements %%%%%%%%%%%%%%%%
                    cbk_offdiag = ndof*(ELEM(ELEM(e).siblings(iSib)).nodes(jnod)-1) + 1;
                    for idof=1:ndof
                        for jdof=1:ndof
                            I(next+1) = rg + (idof-1); J(next+1) = cbk_offdiag + (jdof-1);
                            V(next+1) = ke_offdiag(re+(idof-1),ce+(jdof-1),iSib);
                            next = next+1;
                        end
                    end
                end
            end
        end
        %%%%%%%%% Assemble global force vector %%%%%%%%%%%%
        fext(rg:(rg+ndof-1)) = fext(rg:(rg+ndof-1)) + fe(re:(re+ndof-1)) +...
            + fe_gam(re:(re+ndof-1)) + fe_bdy(re:(re+ndof-1));
    end
end
bigk = sparse(I,J,V,num_dof,num_dof);
