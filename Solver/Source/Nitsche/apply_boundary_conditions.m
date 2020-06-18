function [ke_bdy,fe_bdy,PARAMS] = apply_boundary_conditions(elem,NODE,PARAMS)

ke_bdy = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
fe_bdy = zeros(PARAMS.nlink*PARAMS.ndof,1);

%%% Loop over outer boundaries of the domain.
% Currently looping order is: [Xmin Xmax Ymin Ymax]
for iBdy=1:PARAMS.num_out_bdy
    
    %%%% If current boundary is a DB and if current element is a part of it
    %%%% then apply DBCs weakly.
    if(PARAMS.DB(iBdy,1)==1 && ((~isempty(elem.boundX)>0 && elem.boundX==PARAMS.DB(iBdy,2) && (iBdy-1)*(iBdy-2)==0) ||...
            (~isempty(elem.boundY)>0 && elem.boundY==PARAMS.DB(iBdy,2)&& (iBdy-3)*(iBdy-4)==0 )))
        bdy_normal = get_bdy_normal(iBdy);
        if(iBdy==1||iBdy==2)
            xint = elem.bdySegXx; yint = elem.bdySegXy;
        else
            xint = elem.bdySegYx; yint = elem.bdySegYy;
        end
        
        %drawArrow = @(x,y,props) quiver( x(1),y(1),(x(2)-x(1))*0.1,(y(2)-y(1))*0.1,0, props{:} );     
        %drawArrow([0.5*(xint(1)+xint(2)) 0.5*(xint(1)+xint(2))+bdy_normal(1)],[0.5*(yint(1)+yint(2)) 0.5*(yint(1)+yint(2))+bdy_normal(2)],{'MaxHeadSize',0.8,'Color','r','LineWidth',1});
        
        
        if(~isempty(xint))
            [ke_dirichlet,fe_dirichlet] = nitsche_terms_bdy(elem,bdy_normal,NODE,xint,yint,PARAMS);
        else
            ke_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
            fe_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,1);
        end
        fe_neumann = zeros(PARAMS.nlink*PARAMS.ndof,1);
        fe_robin = zeros(PARAMS.nlink*PARAMS.ndof,1);
        ke_robin = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
        
        %%%% If current boundary is a NB and if current element is a part of it
        %%%% then build RHS to account for NBCs
    elseif(PARAMS.NB(iBdy,1)==1 && ((~isempty(elem.boundX)>0 && elem.boundX==PARAMS.NB(iBdy,2) && (iBdy-1)*(iBdy-2)==0 ) || ...
            (~isempty(elem.boundY)>0 && elem.boundY==PARAMS.NB(iBdy,2) && (iBdy-3)*(iBdy-4)==0 )))
        bdy_normal = get_bdy_normal(iBdy);
        if(iBdy==1||iBdy==2)
            xint = elem.bdySegXx; yint = elem.bdySegXy;
        else
            xint = elem.bdySegYx; yint = elem.bdySegYy;
        end
        
        %drawArrow = @(x,y,props) quiver( x(1),y(1),(x(2)-x(1))*0.1,(y(2)-y(1))*0.1,0, props{:} );     
        %drawArrow([0.5*(xint(1)+xint(2)) 0.5*(xint(1)+xint(2))+bdy_normal(1)],[0.5*(yint(1)+yint(2)) 0.5*(yint(1)+yint(2))+bdy_normal(2)],{'MaxHeadSize',0.8,'Color','b','LineWidth',1});
        
        if(~isempty(xint))
            [fe_neumann] = neumann_force_terms(elem,bdy_normal,NODE,xint,yint,PARAMS);
        else
            fe_neumann = zeros(PARAMS.nlink*PARAMS.ndof,1);
        end
       
        ke_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
        fe_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,1);
        fe_robin = zeros(PARAMS.nlink*PARAMS.ndof,1);
        ke_robin = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
    else
        ke_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
        fe_dirichlet = zeros(PARAMS.nlink*PARAMS.ndof,1);
        fe_neumann = zeros(PARAMS.nlink*PARAMS.ndof,1);
        fe_robin = zeros(PARAMS.nlink*PARAMS.ndof,1);
        ke_robin = zeros(PARAMS.nlink*PARAMS.ndof,PARAMS.nlink*PARAMS.ndof);
    end
    %%%% If current boundary is a RB and if current element is a part of it
    %%%% then build LHS and RHS to account for RBCs. Only built for scalar
    %%%% problems currently.
    if(PARAMS.ndof==1)
        if(PARAMS.RB(iBdy,1)==1 && ((~isempty(elem.boundX)>0 && elem.boundX==PARAMS.RB(iBdy,2) && (iBdy-1)*(iBdy-2)==0 ) || ...
                (~isempty(elem.boundY)>0 && elem.boundY==PARAMS.RB(iBdy,2) && (iBdy-3)*(iBdy-4)==0 )))
            %                 bdy_normal = get_bdy_normal(iBdy);
            if(iBdy==1||iBdy==2)
                xint = elem.bdySegXx; yint = elem.bdySegXy;
            else
                xint = elem.bdySegYx; yint = elem.bdySegYy;
            end
            if(~isempty(xint))
                [ke_robin,fe_robin] = robin_terms_bdy(elem,NODE,xint,yint,PARAMS);
            end
        end
    end
    
    ke_bdy = ke_bdy + ke_dirichlet + ke_robin;
    fe_bdy = fe_bdy+fe_dirichlet+fe_neumann + fe_robin;
end
end

%%%% Function to get normal:
% For a square/rectangular domain currently. So we know the normals
% apriori, if needed this can be changed to calculate the normals instead.
function [bdy_normal] = get_bdy_normal(iBdy)
if(iBdy==1), bdy_normal = [-1, 0];elseif(iBdy==2), bdy_normal = [1 0];
elseif(iBdy==3), bdy_normal = [0 -1];elseif(iBdy==4), bdy_normal = [0 1];
end
end