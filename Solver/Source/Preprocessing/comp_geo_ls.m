function[ELEM,NODE,VERT,VERT_ELE,PARAMS] = comp_geo_ls(x,y,node,ls,PARAMS)

fprintf('*************************************************************************** \n')
fprintf('      Pre-processing: Computational geometry for enriched formulation       \n')

numnod = size(x,1); numele = size(node,2);
nlink = 6; %6
nvert = 3;
nlset = 3;

global p;
global d;

p = 2; %polynomial order
d = 2; %dimension

% Computational domain
maxX=max(x); minX = min(x); lengthX=abs(maxX-minX);
maxY=max(y); minY = min(y); lengthY=abs(maxY-minY);

% For Nitsche's stabilization parameter
yieldParam = [];
[E,nu,yieldParam] = mat_param(PARAMS);

PARAMS.num_bg_nod = numnod;
PARAMS.num_bg_ele = numele;
PARAMS.sum_num_vv_e = 0;
PARAMS.nlink = nlink;
PARAMS.nlset = nlset;
PARAMS.nvert = nvert;
PARAMS.maxX = maxX; PARAMS.minX = minX; PARAMS.lengthX = lengthX;
PARAMS.maxY = maxY; PARAMS.minY = minY; PARAMS.lengthY = lengthY;

PARAMS.numele = 0;

%% Initialize number of virtual nodes
nb_virt_nod = 0;
nodepartochild = zeros(numnod,1);

%% Allocate memory in blocks
block_size = floor(numele + 0.1*numele);
alloc_size = block_size;
valloc_size = block_size;
valloc_ele_size = 3*block_size;

%% Structure ELEM
% ELEM(i).nodes:       modified connectivity of element i
% ELEM(i).domain:      index of the domain to which element i belongs
% ELEM(i).volume:      computational volume of element i
% ELEM(i).siblings:    other partial element(s) from the same parent
% ELEM(i).normals:     normal(s) of embedded interface(s) if any
% ELEM(i).vertices:    list of global vertex indices in element i
% ELEM(i).stab:        (list of) stabilization parameters for Nitsche's method at the interface
% ELEM(i).stab_dirich: (list of) stabilization parameter for Nitsche's method at Dirichlet boundaries
%                      ELEM(i).stab_dirich(1) corresponds to X boundary and
%                      ELEM(i).stab_dirich(2) corresponds to Y boundary. If
%                      either of these is empty then ELEM(i) does not
%                      belong to the corresponding boundary.
% ELEM(i).localVV:     list of vital vertices whose support are defined over elem i
% ELEM(i).activEP:     list of nodes of elem i which are active end-points
% ELEM(i).gamma:       weight for the bracket operator computed as a weighted average
% ELEM(i).boundX:      extremal value in x direction defining the outer boundaries of ELEM(i)
% ELEM(i).boundY:      extremal value in Y direction defining the outer boundaries of ELEM(i)
%                      rmk: if boundX and boundY are empty ELEM(i) doesn't
%                      belong to the corresponding boundary
% ELEM(i).uJumpPl:     Plastic jump at the interface for ELEM(i) at both gauss
%                      points and both load steps
% ELEM(i).MatParam:    Young's modulus and Poisson's ratio for ELEM(i)
% ELEM(i).YieldTrac:   Interfacial yield traction and Coulomb's friction
%                      coefficient for ELEM(i)
% ELEM(i).bdySegXx:    x-coordinates of the end-points of the segment of
%                      ELEM(i) that belongs to an extremal boundary in x-direction
% ELEM(i).bdySegXy:    y-coordinates of the end-points of the segment of
%                      ELEM(i) that belongs to an extremal boundary in x-direction
% ELEM(i).bdySegYx:    x-coordinates of the end-points of the segment of
%                      ELEM(i) that belongs to an extremal boundary in y-direction
% ELEM(i).bdySegYy:    y-coordinates of the end-points of the segment of
%                      ELEM(i) that belongs to an extremal boundary in y-direction

if(~PARAMS.debug && ~PARAMS.plotInt && ~PARAMS.plotIntTrac && ~PARAMS.refinement)
    ELEM(alloc_size) = struct('nodes', [], 'domain', [], 'volume', [], 'siblings', [],...
        'normals', [], 'vertices', [], 'stab', [], 'stab_dirich', [], 'localVV', [],...
        'activEP', [],'gamma', [], 'boundX', [], 'boundY', [], 'uJumpPl', [], 'MatParam', [],...
        'YieldTrac', [], 'bdySegXx',[], 'bdySegXy',[], 'bdySegYx',[], 'bdySegYy',[], 'stick', []);
else
    ELEM(alloc_size) = struct('nodes', [], 'domain', [], 'volume', [], 'siblings', [],...
        'normals', [], 'vertices', [],'vertices_ele', [], 'vertices_ele_virt',[], 'stab', [], 'stab_dirich',[], 'localVV', [],...
        'activEP', [], 'gamma', [], 'boundX', [], 'boundY', [], 'uJumpPl', [],...
        'MatParam', [], 'YieldTrac', [], 'tracInt', [], 'uJump', [], 'Pos', [],...
        'bdySegXx',[], 'bdySegXy',[], 'bdySegYx',[], 'bdySegYy',[], 'stick', []);
end

%% Structure VERT
% VERT(i).X:     x-coordinate of vertex i
% VERT(i).Y:     y-coordinate of vertex i
% VERT(i).vital: bool to indicate if vertex i is vital
VERT(valloc_size) = struct('X', [], 'Y', [], 'vital', []);
VERT_ELE(valloc_ele_size) = struct('X', [], 'Y', [], 'ls',[]);

csize = 0;
vsize = 0;
vsize_ele = 0;
cut_edge_conn = zeros(2,3*numele);
n_cut_ele = 0;
n_cut_edge = 0;

for e=1:numele
    
    %[int_flag] = intersection_test(ls,node,e,6);    
    [int_flag] = intersection_test2(node,e,x,y,PARAMS);
    
    %%%% Hack for Liu, Borja example to terminate a crack %%%%
    %%%% WARNING: Has to be made more general later if want to model cracks!
    if(PARAMS.LiuBorja == 1)
        if(~isempty(find(x(node(1:3,e))<PARAMS.CrackMin, 1)) || ~isempty(find(x(node(1:3,e))>PARAMS.CrackMax, 1)))
            int_flag = 0;
        end
    end
    
    if(int_flag==1)
        if(~PARAMS.refinement)
            PARAMS.numele = PARAMS.numele+2;
            n_cut_ele = n_cut_ele+1;
            %%%%% Create a map between parent and child nodes
            lse = ls(node(1:nlink,e));
            for i=1:nlink
                if(nodepartochild(node(i,e))==0)				%%%% check to identify if parent node is already associated with a child
                    nb_virt_nod = nb_virt_nod + 1;				%%%% if not, create a child
                    nodepartochild(node(i,e)) = numnod + nb_virt_nod;
                end
            end
            
            %%%%% Create element connectivity for cut elements
            for i=1:nlink
                if(lse(i)<0)
                    ELEM(csize+1).nodes(i) = node(i,e);			%%% if node belongs to Omega^-, retain it
                    ELEM(csize+2).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                else
                    ELEM(csize+2).nodes(i) = node(i,e);			%%% if node belongs to Omega^+, retain it
                    ELEM(csize+1).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                end
            end
            
            VERT_ELE(vsize_ele+1).X = x(node(1,e));
            VERT_ELE(vsize_ele+2).X = x(node(2,e));
            VERT_ELE(vsize_ele+3).X = x(node(3,e));
            
            VERT_ELE(vsize_ele+1).Y = y(node(1,e));
            VERT_ELE(vsize_ele+2).Y = y(node(2,e));
            VERT_ELE(vsize_ele+3).Y = y(node(3,e));
            
            ELEM(csize+1).vertices_ele(1) = vsize_ele+1;
            ELEM(csize+1).vertices_ele(2) = vsize_ele+2;
            ELEM(csize+1).vertices_ele(3) = vsize_ele+3;
            
            ELEM(csize+2).vertices_ele(1) = vsize_ele+1;
            ELEM(csize+2).vertices_ele(2) = vsize_ele+2;
            ELEM(csize+2).vertices_ele(3) = vsize_ele+3;
            
            %%%%% Create relationship between siblings
            ELEM(csize+1).siblings = csize+2;
            ELEM(csize+2).siblings = csize+1;
            
            %%%%% Get vertices and normals
            [xint,yint,coeff,nb_int,edge_id] = intersection_points(node,x,y,ls,e);
            [VERT,ELEM,cut_edge_conn,vsize,n_cut_edge] = get_vertices(VERT,ELEM,e,vsize,csize,node,cut_edge_conn,edge_id,nb_int,xint,yint,n_cut_edge);
            ELEM(csize+1).normals = -normal(coeff);
            ELEM(csize+2).normals = normal(coeff);
            
            %%%%% Get computational volume for children and domain
            [ELEM(csize+2).volume,ELEM(csize+1).volume] = sub_areas(node,e,x,y,ls,xint,yint);
            ELEM(csize+1).domain = 1;
            ELEM(csize+2).domain = 2;
            
            %%%% Store material parameters in struct ELEM
            ELEM(csize+1).MatParam(1) = E(ELEM(csize+1).domain);
            ELEM(csize+2).MatParam(1) = E(ELEM(csize+2).domain);
            if(PARAMS.ndof==2)
                ELEM(csize+1).MatParam(2) = nu(ELEM(csize+1).domain);
                ELEM(csize+2).MatParam(2) = nu(ELEM(csize+2).domain);
            end
            
            %%%% Get stabilization parameter for the interface:
            % \alpha_int = 2*C_int^2 regardless of whether a Dirichlet boundary
            % exists within element e.
            length_e = sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);
            [ELEM] = get_stab_param_int(PARAMS,ELEM,csize,length_e);
            
            %%%% Test whether domain boundary is a part of element boundary
            [ELEM] = test_boundary_elem(ELEM,PARAMS,VERT,VERT_ELE,x,y,node,e,int_flag,csize);
            
            %%%%% Get stabilization parameter for the Dirichlet surface
            %     For the first child
            [ELEM(csize+1)] = get_stab_param_dirich(ELEM(csize+1),PARAMS);
            %     For the second child
            [ELEM(csize+2)] = get_stab_param_dirich(ELEM(csize+2),PARAMS);
            
            
            %%%%% Store interface quantities: ujump and traction in struct ELEM:
            %%%%% uJump(l,m) & tracInt(l,m)
            %%%%% l--> Gauss points; m--> direction (n, tau);
            if(PARAMS.plotInt)
                ELEM(csize+1).uJump = zeros(2,2); ELEM(csize+2).uJump = zeros(2,2);
                ELEM(csize+1).tracInt = zeros(2,2); ELEM(csize+2).tracInt = zeros(2,2);
                ELEM(csize+1).Pos = zeros(2,2); ELEM(csize+2).Pos = zeros(2,2);
            end
            
            vsize_ele = vsize_ele + 3;
            csize = csize+2;
        else
            [node, x, y , ls, num_sub_elements] = hierarchical_tri( e, node, x, y , ls, PARAMS);
            cut_edge_conn = [cut_edge_conn, zeros(2,3*num_sub_elements)];
            %PARAMS.num_bg_ele = PARAMS.num_bg_ele + num_sub_elements-1;
            PARAMS.numele = PARAMS.numele+2;
            
            %%%%% Create a map between parent and child nodes
            for i=1:nlink
                if(nodepartochild(node(i,e))==0)				%%%% check to identify if parent node is already associated with a child
                    nb_virt_nod = nb_virt_nod + 1;				%%%% if not, create a child
                    nodepartochild(node(i,e)) = numnod + nb_virt_nod;
                    hold on;
                    %plot(x(node(i,e)),y(node(i,e)),'ro','MarkerSize',10,'MarkerFaceColor','r');
                end
            end
            
            volume_domain1 = 0.0;
            volume_domain2 = 0.0;
            length_interface = 0.0;
            temp_csize = csize;
            length_boundary1_x = 0.0;
            length_boundary1_y = 0.0;
            length_boundary2_x = 0.0;
            length_boundary2_y = 0.0;
            volume_domain1_bdry = 0.0;
            volume_domain2_bdry = 0.0;
            
            for ei = size(node,2)-num_sub_elements+1:size(node,2)
                
                [int_flag_sub] = intersection_test(ls,node,ei,6);
                %[int_flag_sub] = intersection_test2(node,ei,x,y,PARAMS);
                
                if(int_flag_sub==1)
                    n_cut_ele = n_cut_ele+1;
                    
                    lsei = ls(node(1:nlink,ei));
                    lse = ls(node(1:nlink,e));
                    %%%%% Create element connectivity for cut elements
                    for i=1:nlink
                        if(lsei(i)<0)                            
                            ELEM(csize+1).vertices_ele_virt(i) = 1; % 1 real node
                            ELEM(csize+2).vertices_ele_virt(i) = 2; % 2 virtural node               
                        else
                            ELEM(csize+1).vertices_ele_virt(i) = 2; %
                            ELEM(csize+2).vertices_ele_virt(i) = 1; %
                        end
                    end
                    
                    for i=1:nlink
                        if(lse(i)<0)
                            ELEM(csize+1).nodes(i) = node(i,e);			%%% if node belongs to Omega^-, retain it
                            ELEM(csize+2).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                     else
                            ELEM(csize+2).nodes(i) = node(i,e);			%%% if node belongs to Omega^+, retain it
                            ELEM(csize+1).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                        end
                    end
                                        
                    %%%%% Create relationship between siblings
                    ELEM(csize+1).siblings = csize+2;
                    ELEM(csize+2).siblings = csize+1;
                    
                    %%%%% Get vertices and normals
                    [xint,yint,coeff,nb_int,edge_id] = intersection_points(node,x,y,ls,ei);
          
                    [VERT,ELEM,cut_edge_conn,vsize,n_cut_edge] = get_vertices(VERT,ELEM,ei,vsize,csize,node,cut_edge_conn,edge_id,nb_int,xint,yint,n_cut_edge);
                    %plot([x(node(:,ei));x(node(1,ei))],[y(node(:,ei));y(node(1,ei))]);
                    %hold on;
                    %plot(x(cut_edge_conn(1:2,n_cut_edge-1)),y(cut_edge_conn(1:2,n_cut_edge-1)),'LineWidth',3);
                    %plot(x(cut_edge_conn(1:2,n_cut_edge)),y(cut_edge_conn(1:2,n_cut_edge)),'LineWidth',3);
                    
                    
                    ELEM(csize+1).normals = -normal(coeff);
                    ELEM(csize+2).normals = normal(coeff);
                    
                    %%%%% Get computational volume for children and domain
                    [ELEM(csize+2).volume,ELEM(csize+1).volume] = sub_areas(node,ei,x,y,ls,xint,yint);
                    ELEM(csize+1).domain = 1;
                    ELEM(csize+2).domain = 2;
                    
                    volume_domain1 = volume_domain1 + ELEM(csize+1).volume;
                    volume_domain2 = volume_domain2 + ELEM(csize+2).volume;
                    
                    %%%% Store material parameters in struct ELEM
                    ELEM(csize+1).MatParam(1) = E(ELEM(csize+1).domain);
                    ELEM(csize+2).MatParam(1) = E(ELEM(csize+2).domain);
                    if(PARAMS.ndof==2)
                        ELEM(csize+1).MatParam(2) = nu(ELEM(csize+1).domain);
                        ELEM(csize+2).MatParam(2) = nu(ELEM(csize+2).domain);
                    end
                    
                    VERT_ELE(vsize_ele+1).X = x(node(1,ei));
                    VERT_ELE(vsize_ele+2).X = x(node(2,ei));
                    VERT_ELE(vsize_ele+3).X = x(node(3,ei));
                    
                    VERT_ELE(vsize_ele+1).Y = y(node(1,ei));
                    VERT_ELE(vsize_ele+2).Y = y(node(2,ei));
                    VERT_ELE(vsize_ele+3).Y = y(node(3,ei));
                    
                    ELEM(csize+1).vertices_ele(1) = vsize_ele+1;
                    ELEM(csize+1).vertices_ele(2) = vsize_ele+2;
                    ELEM(csize+1).vertices_ele(3) = vsize_ele+3;
                    
                    ELEM(csize+2).vertices_ele(1) = vsize_ele+1;
                    ELEM(csize+2).vertices_ele(2) = vsize_ele+2;
                    ELEM(csize+2).vertices_ele(3) = vsize_ele+3;
                    
                    %%%% Get stabilization parameter for the interface:
                    % \alpha_int = 2*C_int^2 regardless of whether a Dirichlet boundary
                    % exists within element e.
                    length_e = sqrt((xint(2)-xint(1))^2 + (yint(2)-yint(1))^2);
                    length_interface = length_interface + length_e;
                   
                    %[ELEM] = get_stab_param_int(PARAMS,ELEM,csize,length_e);
                    
                    %%%% Test whether domain boundary is a part of element boundary
                    [ELEM] = test_boundary_elem(ELEM,PARAMS,VERT,VERT_ELE,x,y,node,ei,int_flag_sub,csize);
                    
                    %%%%% Get stabilization parameter for the Dirichlet surface
                    %     For the first child
                    %[ELEM(csize+1)] = get_stab_param_dirich(ELEM(csize+1),PARAMS);
                    %     For the second child
                    % [ELEM(csize+2)] = get_stab_param_dirich(ELEM(csize+2),PARAMS);
                    
                    
                    if(~isempty(ELEM(csize+1).boundX))
                        length_e = sqrt((ELEM(csize+1).bdySegXx(1)-ELEM(csize+1).bdySegXx(2))^2 + (ELEM(csize+1).bdySegXy(1)-ELEM(csize+1).bdySegXy(2))^2);
                        length_boundary1_x = length_boundary1_x + length_e;
                    end
                    if(~isempty(ELEM(csize+1).boundY))
                        length_e = sqrt((ELEM(csize+1).bdySegYx(1)-ELEM(csize+1).bdySegYx(2))^2 + (ELEM(csize+1).bdySegYy(1)-ELEM(csize+1).bdySegYy(2))^2);
                        length_boundary1_y = length_boundary1_y + length_e;
                    end
                    
                    if(~isempty(ELEM(csize+2).boundX))
                        length_e = sqrt((ELEM(csize+2).bdySegXx(1)-ELEM(csize+2).bdySegXx(2))^2 + (ELEM(csize+2).bdySegXy(1)-ELEM(csize+2).bdySegXy(2))^2);
                        length_boundary2_x = length_boundary2_x + length_e;
                    end
                    if(~isempty(ELEM(csize+2).boundY))
                        length_e = sqrt((ELEM(csize+2).bdySegYx(1)-ELEM(csize+2).bdySegYx(2))^2 + (ELEM(csize+2).bdySegYy(1)-ELEM(csize+2).bdySegYy(2))^2);
                        length_boundary2_y = length_boundary2_y + length_e;
                    end
                    
                    volume_domain1_bdry = volume_domain1_bdry + ELEM(csize+1).volume;
                    volume_domain2_bdry = volume_domain2_bdry + ELEM(csize+2).volume;
                    
                    %%%%% Store interface quantities: ujump and traction in struct ELEM:
                    %%%%% uJump(l,m) & tracInt(l,m)
                    %%%%% l--> Gauss points; m--> direction (n, tau);
                    if(PARAMS.plotInt)
                        ELEM(csize+1).uJump = zeros(2,2); ELEM(csize+2).uJump = zeros(2,2);
                        ELEM(csize+1).tracInt = zeros(2,2); ELEM(csize+2).tracInt = zeros(2,2);
                        ELEM(csize+1).Pos = zeros(2,2); ELEM(csize+2).Pos = zeros(2,2);
                    end
                    
                    vsize_ele = vsize_ele + 3;
                    
                    csize = csize+2;
                else
                    lsei_sum = sum(ls(node(1:nlink,ei)));
                    lse = ls(node(1:nlink,e));
                    if(lsei_sum < 0.0)
                        for i=1:nlink
                            if(lse(i)<0)
                                ELEM(csize+1).nodes(i) = node(i,e);			%%% if node belongs to Omega^-, retain it
                            else
                                ELEM(csize+1).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                            end
                        end
                    else
                        for i=1:nlink
                            if(lse(i)>0)
                                ELEM(csize+1).nodes(i) = node(i,e);			%%% if node belongs to Omega^-, retain it
                            else
                                ELEM(csize+1).nodes(i) = nodepartochild(node(i,e));	%%% and replace it with its child in the sibling
                            end
                        end
                    end
                    
                    ELEM(csize+1).vertices_ele_virt(1:3) = 1;
                    
                    ELEM(csize+1).volume = polyarea(x(node(1:3,ei)), y(node(1:3,ei))); %TO BE CHECKED polyarea(x(node(1:nlink,ei)), y(node(1:nlink,ei)))
                    set = find(ls(node(1:nlink,ei))<0);
                    if (set>0)
                        ELEM(csize+1).domain = 1;
                        volume_domain1 = volume_domain1 + ELEM(csize+1).volume;
                    else
                        ELEM(csize+1).domain = 2;
                        volume_domain2 = volume_domain2 + ELEM(csize+1).volume;
                    end
                    
                    %%%% Store material parameters in struct ELEM
                    ELEM(csize+1).MatParam(1) = E(ELEM(csize+1).domain);
                    if(PARAMS.ndof==2)
                        ELEM(csize+1).MatParam(2) = nu(ELEM(csize+1).domain);
                    end
                    
                    VERT_ELE(vsize_ele+1).X = x(node(1,ei));
                    VERT_ELE(vsize_ele+2).X = x(node(2,ei));
                    VERT_ELE(vsize_ele+3).X = x(node(3,ei));
                    
                    VERT_ELE(vsize_ele+1).Y = y(node(1,ei));
                    VERT_ELE(vsize_ele+2).Y = y(node(2,ei));
                    VERT_ELE(vsize_ele+3).Y = y(node(3,ei));
                    
                    ELEM(csize+1).vertices_ele(1) = vsize_ele+1;
                    ELEM(csize+1).vertices_ele(2) = vsize_ele+2;
                    ELEM(csize+1).vertices_ele(3) = vsize_ele+3;
                    
                    %%%% Test whether domain boundary is a part of element boundary
                    [ELEM] = test_boundary_elem(ELEM,PARAMS,VERT,VERT_ELE,x,y,node,ei,int_flag_sub,csize);
                    
                    if (ELEM(csize+1).domain == 1)
                        if(~isempty(ELEM(csize+1).boundX))
                            length_e = sqrt((ELEM(csize+1).bdySegXx(1)-ELEM(csize+1).bdySegXx(2))^2 + (ELEM(csize+1).bdySegXy(1)-ELEM(csize+1).bdySegXy(2))^2);
                            length_boundary1_x = length_boundary1_x + length_e;
                        end
                        if(~isempty(ELEM(csize+1).boundY))
                            length_e = sqrt((ELEM(csize+1).bdySegYx(1)-ELEM(csize+1).bdySegYx(2))^2 + (ELEM(csize+1).bdySegYy(1)-ELEM(csize+1).bdySegYy(2))^2);
                            length_boundary1_y = length_boundary1_y + length_e;
                        end
                        volume_domain1_bdry = volume_domain1_bdry + ELEM(csize+1).volume;
                    else
                        if(~isempty(ELEM(csize+1).boundX))
                            length_e = sqrt((ELEM(csize+1).bdySegXx(1)-ELEM(csize+1).bdySegXx(2))^2 + (ELEM(csize+1).bdySegXy(1)-ELEM(csize+1).bdySegXy(2))^2);
                            length_boundary2_x = length_boundary2_x + length_e;
                        end
                        if(~isempty(ELEM(csize+1).boundY))
                            length_e = sqrt((ELEM(csize+1).bdySegYx(1)-ELEM(csize+1).bdySegYx(2))^2 + (ELEM(csize+1).bdySegYy(1)-ELEM(csize+1).bdySegYy(2))^2);
                            length_boundary2_y = length_boundary2_y + length_e;
                        end
                        volume_domain2_bdry = volume_domain2_bdry + ELEM(csize+1).volume;
                    end
                    
                    %%%%% Get stabilization parameter for Dirichlet boundary
                    %[ELEM(csize+1)] = get_stab_param_dirich(ELEM(csize+1),PARAMS);
                    
                    vsize_ele = vsize_ele + 3;
                    
                    csize = csize+1;
                end
            end
                        
            for ei = size(node,2)-num_sub_elements+1:size(node,2)
                [int_flag_sub] = intersection_test(ls,node,ei,6);
                %[int_flag_sub] = intersection_test2(node,ei,x,y,PARAMS);
                if(int_flag_sub==1)
                    [ELEM] = get_stab_param_int_refine(PARAMS,ELEM,temp_csize,length_interface,volume_domain1,volume_domain2);
                    [ELEM] = get_stab_param_dirich_refine(ELEM,PARAMS, temp_csize+1, length_boundary1_x,length_boundary1_y, volume_domain1_bdry);
                    [ELEM] = get_stab_param_dirich_refine(ELEM,PARAMS, temp_csize+2, length_boundary2_x,length_boundary2_y, volume_domain2_bdry);
                    temp_csize = temp_csize + 2;
                else
                    if(ELEM(temp_csize).domain == 1)
                        [ELEM] = get_stab_param_dirich_refine(ELEM,PARAMS, temp_csize+1, length_boundary1_x,length_boundary1_y, volume_domain1_bdry);
                    else
                        [ELEM] = get_stab_param_dirich_refine(ELEM,PARAMS, temp_csize+1, length_boundary2_x,length_boundary2_y, volume_domain2_bdry);
                    end
                    temp_csize = temp_csize + 1;
                end
            end
            
        end
    else
        PARAMS.numele = PARAMS.numele+1;
        %%%%% Connectivity and computational volume for uncut elements
        ELEM(csize+1).nodes = node(1:nlink,e)';
        ELEM(csize+1).volume = polyarea(x(node(1:nlink,e)), y(node(1:nlink,e)));
        set = find(ls(node(1:nlink,e))<0);
        if (set>0)
            ELEM(csize+1).domain = 1;
        else
            ELEM(csize+1).domain = 2;
        end
        
        %%%% Store material parameters in struct ELEM
        ELEM(csize+1).MatParam(1) = E(ELEM(csize+1).domain);
        if(PARAMS.ndof==2)
            ELEM(csize+1).MatParam(2) = nu(ELEM(csize+1).domain);
        end
        
        ELEM(csize+1).vertices_ele_virt(1:3) = 1;
        
        VERT_ELE(vsize_ele+1).X = x(node(1,e));
        VERT_ELE(vsize_ele+2).X = x(node(2,e));
        VERT_ELE(vsize_ele+3).X = x(node(3,e));
        
        VERT_ELE(vsize_ele+1).Y = y(node(1,e));
        VERT_ELE(vsize_ele+2).Y = y(node(2,e));
        VERT_ELE(vsize_ele+3).Y = y(node(3,e));
        
        ELEM(csize+1).vertices_ele(1) = vsize_ele+1;
        ELEM(csize+1).vertices_ele(2) = vsize_ele+2;
        ELEM(csize+1).vertices_ele(3) = vsize_ele+3;
        
        %%%% Test whether domain boundary is a part of element boundary
        [ELEM] = test_boundary_elem(ELEM,PARAMS,VERT,VERT_ELE,x,y,node,e,int_flag,csize);
        
        %%%%% Get stabilization parameter for Dirichlet boundary
        [ELEM(csize+1)] = get_stab_param_dirich(ELEM(csize+1),PARAMS);
        
        vsize_ele = vsize_ele + 3;
        csize = csize+1;
    end
    
    %% Memory allocation for ELEM & VERT
    if(csize + 0.1*alloc_size>alloc_size)
        alloc_size = alloc_size + 4^(PARAMS.refine_levels)*block_size;
        if(~PARAMS.debug && ~PARAMS.plotInt && ~PARAMS.plotIntTrac && ~PARAMS.refinement)
            ELEM(alloc_size) = struct('nodes', [], 'domain', [], 'volume', [],...
                'siblings', [], 'normals', [], 'vertices', [], 'stab', [],...
                'stab_dirich', [], 'localVV', [], 'activEP', [], 'gamma', [],...
                'boundX', [], 'boundY', [], 'uJumpPl', [], 'MatParam', [],...
                'YieldTrac', [], 'bdySegXx',[], 'bdySegXy',[], 'bdySegYx',[], 'bdySegYy',[], 'stick', []);
        else
            ELEM(alloc_size) = struct('nodes', [], 'domain', [], 'volume', [], 'siblings', [],...
                'normals', [], 'vertices', [], 'vertices_ele',  [],'vertices_ele_virt',[], 'stab', [], 'stab_dirich',[], 'localVV', [],...
                'activEP', [], 'gamma', [], 'boundX', [], 'boundY', [], 'uJumpPl', [],...
                'MatParam', [], 'YieldTrac', [], 'tracInt', [], 'uJump', [], 'Pos', [],...
                'bdySegXx',[], 'bdySegXy',[], 'bdySegYx',[], 'bdySegYy',[], 'stick', []);
        end
    end
    if(vsize + 0.1*valloc_size>valloc_size)
        valloc_size = valloc_size + block_size*4^(PARAMS.refine_levels);
        VERT(valloc_size) = struct('X', [], 'Y', [], 'vital', []);
    end
    
    if(vsize_ele + 0.1*valloc_ele_size>valloc_ele_size)
        valloc_ele_size = valloc_ele_size + block_size*3*4^(PARAMS.refine_levels);
        VERT_ELE(valloc_ele_size) = struct('X', [], 'Y', [], 'ls', []);
    end
end

cut_edge_conn(:,(n_cut_edge+1):end) = [];
ELEM(csize+1:end) = [];
VERT(vsize+1:end) = [];
VERT_ELE(vsize_ele+1:end) = [];

%% update level set function
for i = 1:size(VERT_ELE,2)
    %VERT_ELE(i).ls = VERT_ELE(i).X - (PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin);
   % VERT_ELE(i).ls = (VERT_ELE(i).X)^2+ (VERT_ELE(i).Y)^2- 1/4;
       VERT_ELE(i).ls = get_level_set([VERT_ELE(i).X],[VERT_ELE(i).Y],PARAMS);

end


% hold on;
% plot_mesh([x y],node','T3','g.-',1.2);
% plot_mesh([x y],cut_edge_conn','L2','r.-',1.2);
% for i = 1:length(cut_edge_conn)
%     plot(x(cut_edge_conn(:,i)),y(cut_edge_conn(:,i)),'Color','b','LineWidth',3);
%     hold on;
% end
% 
% for i = 1:length(ELEM)
%     xe = [];
%     ye = [];
%     for j = 1:4
%         if (j~=4)
%             xe(j) = VERT_ELE(ELEM(i).vertices_ele(j)).X;
%             ye(j) = VERT_ELE(ELEM(i).vertices_ele(j)).Y;
%         else
%             xe(j) = VERT_ELE(ELEM(i).vertices_ele(1)).X;
%             ye(j) = VERT_ELE(ELEM(i).vertices_ele(1)).Y;
%         end
%     end
%     hold on;
%     plot(xe,ye);
% end

% theta = 0:0.01:2*pi;
% xo = 0.5*cos(theta) ;
% yo = 0.5*sin(theta) ;
% plot(xo,yo,'r-','Linewidth',1.2);

num_tot_nod = numnod + nb_virt_nod;

%% Structure NODE
% NODE(i).X:          x-coordinate of node i
% NODE(i).Y:          y-coordinate of node i
% NODE(i).domain:     index of the domain to which node i belongs
% NODE(i).ls:         value of the level-set
% NODE(i).childrenId: index of the chid(ren) of node i
% NODE(i).endPt:      bool to indicate whether node i is an end-point
% NODE(i).vvId:       list of vital vertices connected to end-point i

NODE(num_tot_nod) = struct('X', [], 'Y', [], 'domain', [], 'ls', [],...
    'childrenId', [], 'endPt', [], 'vvId', []);

for inod=1:numnod
    %% Get x-, y-coordinates and level set value
    NODE(inod).X = x(inod); NODE(inod).Y = y(inod); NODE(inod).ls = ls(inod);
    
    %% Identify domain index
    if (ls(inod)>0.0)
        NODE(inod).domain = 2;
    else
        NODE(inod).domain = 1;
    end
    if(nodepartochild(inod)~=0)
        NODE(inod).childrenId = nodepartochild(inod);
        NODE(nodepartochild(inod)).X = x(inod); NODE(nodepartochild(inod)).Y = y(inod);
        NODE(nodepartochild(inod)).ls = ls(inod);
        if(ls(inod)>0.0)
            NODE(nodepartochild(inod)).domain = 1;
        else
            NODE(nodepartochild(inod)).domain = 2;
        end
    end
    
    %% Identify end-points
    set = find(cut_edge_conn==inod);
    if(size(set,1)>1)
        NODE(inod).endPt = 1;
        children = NODE(inod).childrenId;
        if(nodepartochild(inod)~=0)
            NODE(children).endPt = 1;
        end
    else
        NODE(inod).endPt = 0;
        children = NODE(inod).childrenId;
        if(nodepartochild(inod)~=0)
            NODE(children).endPt = 0;
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     FUNCTION DEFINITIONS                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Function to test whether ELEM(csize+...) lies on the outer boundary:1
function [ELEM] = test_boundary_elem(ELEM,PARAMS,VERT,VERT_ELE,x,y,node,e,int_flag,csize)

for i=2:PARAMS.nlset
    for k=1:i-1
        if(abs(x(node(i,e))-PARAMS.maxX)<1.0e-12*PARAMS.lengthX && abs(x(node(k,e))-PARAMS.maxX)<1.0e-12*PARAMS.lengthX)
            [ELEM] = get_bdy_flag_and_bdy_seg_ep(ELEM,VERT,VERT_ELE,PARAMS,csize,x,y,int_flag,'X',PARAMS.maxX);
        elseif(abs(x(node(i,e))-PARAMS.minX)<1.0e-12*PARAMS.lengthX && abs(x(node(k,e))-PARAMS.minX)<1.0e-12*PARAMS.lengthX)
            [ELEM] = get_bdy_flag_and_bdy_seg_ep(ELEM,VERT,VERT_ELE,PARAMS,csize,x,y,int_flag,'X',PARAMS.minX);
        end
        if(abs(y(node(i,e))-PARAMS.maxY)<1.0e-12*PARAMS.lengthY && abs(y(node(k,e))-PARAMS.maxY)<1.0e-12*PARAMS.lengthY)
            [ELEM] = get_bdy_flag_and_bdy_seg_ep(ELEM,VERT,VERT_ELE,PARAMS,csize,x,y,int_flag,'Y',PARAMS.maxY);
        elseif(abs(y(node(i,e))-PARAMS.minY)<1.0e-12*PARAMS.lengthY && abs(y(node(k,e))-PARAMS.minY)<1.0e-12*PARAMS.lengthY)
            [ELEM] = get_bdy_flag_and_bdy_seg_ep(ELEM,VERT,VERT_ELE,PARAMS,csize,x,y,int_flag,'Y',PARAMS.minY);
        end
    end
end

end

%%%% Function to get a stabilization parameter for the interface
function [ELEM] = get_stab_param_int(PARAMS,ELEM,csize,length_e)
global p 
global d
DiagT = eye(3); DiagT(3,3) = 2; % Alex's trick to account for shear strains
Dmat_1 = get_elas_tensor(ELEM(csize+1),PARAMS);
Dmat_2 = get_elas_tensor(ELEM(csize+2),PARAMS);
if(PARAMS.gamFlag==2)
    %%%%%%%%% Robust Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%%%%
    wei_vol = norm(Dmat_2*DiagT)*ELEM(csize+1).volume + ...
        norm(Dmat_1*DiagT)*ELEM(csize+2).volume;
    ELEM(csize+1).gamma = ELEM(csize+1).volume*norm(Dmat_2*DiagT)/wei_vol;
    ELEM(csize+2).gamma = ELEM(csize+2).volume*norm(Dmat_1*DiagT)/wei_vol;
    D_avg = norm(Dmat_1*DiagT)*norm(Dmat_2*DiagT)/wei_vol;
    ELEM(csize+1).stab = 2*length_e*D_avg*(p*(p-1+d)/d);   %\alpha
    ELEM(csize+2).stab = 2*length_e*D_avg*(p*(p-1+d)/d);
elseif(PARAMS.gamFlag==1)
    %%%%%%%%% Hansbo Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%%%%
    vol_ele = ELEM(csize+1).volume + ELEM(csize+2).volume;
    ELEM(csize+1).gamma = ELEM(csize+1).volume/vol_ele;
    ELEM(csize+2).gamma = ELEM(csize+2).volume/vol_ele;
    D_avg = (norm(Dmat_2*DiagT)*ELEM(csize+2).volume...
        + norm(Dmat_1*DiagT)*ELEM(csize+1).volume)/vol_ele^2;
    ELEM(csize+1).stab = 2*length_e*D_avg;
    ELEM(csize+2).stab = 2*length_e*D_avg;
else
    %%%%%%%%% Classical Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%
    ELEM(csize+1).gamma = 0.5;
    ELEM(csize+2).gamma = 0.5;
    D_avg = norm(Dmat_2*DiagT)/ELEM(csize+2).volume + ...
        norm(Dmat_1*DiagT)/ELEM(csize+1).volume;
    ELEM(csize+1).stab = (2/4)*(length_e)*D_avg;
    ELEM(csize+2).stab = (2/4)*(length_e)*D_avg;
end

end


%%%% Function to get a stabilization parameter for the interface
function [ELEM] = get_stab_param_int_refine(PARAMS,ELEM,csize,length_e,volume1,volume2)
global p 
global d
DiagT = eye(3); DiagT(3,3) = 2; % Alex's trick to account for shear strains
Dmat_1 = get_elas_tensor(ELEM(csize+1),PARAMS);
Dmat_2 = get_elas_tensor(ELEM(csize+2),PARAMS);
%%%%%%%%% Robust Nitsche's formulation %%%%%%%%%%%%%%%%%%%%%%%%
wei_vol = norm(Dmat_2*DiagT)*volume1 + ...
    norm(Dmat_1*DiagT)*volume2;
ELEM(csize+1).gamma = volume1*norm(Dmat_2*DiagT)/wei_vol;
ELEM(csize+2).gamma = volume2*norm(Dmat_1*DiagT)/wei_vol;
D_avg = norm(Dmat_1*DiagT)*norm(Dmat_2*DiagT)/wei_vol;
ELEM(csize+1).stab = 2*length_e*D_avg*(p*(p-1+d)/d);   %\alpha
ELEM(csize+2).stab = 2*length_e*D_avg*(p*(p-1+d)/d);
end

%%%% Function to get stabilization parameter for a Dirichlet surface
function [elem] = get_stab_param_dirich(elem,PARAMS)
global p 
global d
DiagT = eye(3); DiagT(3,3) = 2; % Alex's trick to account for shear strains
Dmat = get_elas_tensor(elem,PARAMS);
if(~isempty(elem.boundX))
    length_e = sqrt((elem.bdySegXx(1)-elem.bdySegXx(2))^2 + (elem.bdySegXy(1)-elem.bdySegXy(2))^2);
    elem.stab_dirich(1) = 2*(length_e)*norm(Dmat*DiagT)/elem.volume*(p*(p-1+d)/d);
end
if(~isempty(elem.boundY))
    length_e = sqrt((elem.bdySegYx(1)-elem.bdySegYx(2))^2 + (elem.bdySegYy(1)-elem.bdySegYy(2))^2);
    elem.stab_dirich(2) = 2*(length_e)*norm(Dmat*DiagT)/elem.volume*(p*(p-1+d)/d);
end

end

function [ELEM] = get_stab_param_dirich_refine(ELEM,PARAMS, csize, length_x,length_y, volume)
global p
global d
DiagT = eye(3); DiagT(3,3) = 2; % Alex's trick to account for shear strains
Dmat = get_elas_tensor(ELEM(csize),PARAMS);
if(~isempty(ELEM(csize).boundX))
    ELEM(csize).stab_dirich(1) = 2*(length_x)*norm(Dmat*DiagT)/volume*(p*(p-1+d)/d);
end
if(~isempty(ELEM(csize).boundY))
    ELEM(csize).stab_dirich(2) = 2*(length_y)*norm(Dmat*DiagT)/volume*(p*(p-1+d)/d);
end

end

%% This function tests whether the current element belongs to any external
%% boundary and if it does what are the coordinates of the end-points of
%% the boundary segment. Basically, in the code it defines:
% ELEM.boundX, ELEM.boundY, ELEM.bdySegXx, ELEM.bdySegXy, ELEM.bdySegYx, ELEM.bdySegYy
function [ELEM] = get_bdy_flag_and_bdy_seg_ep(ELEM,VERT,VERT_ELE,PARAMS,csize,x,y,int_flag,iBdy,test_bdy)

if(int_flag==1),locDom = 2;else locDom=1;end

for iLocDom=1:locDom
    [xint, yint] = get_bdy_seg_ep(iBdy,ELEM(csize+iLocDom),test_bdy,x,y,VERT,VERT_ELE,PARAMS);
    if(iBdy=='X'), ELEM(csize+iLocDom).bdySegXx = xint; ELEM(csize+iLocDom).bdySegXy = yint; end
    if(iBdy=='Y'), ELEM(csize+iLocDom).bdySegYx = xint; ELEM(csize+iLocDom).bdySegYy = yint; end
    
    if(iBdy=='X' && ~isempty(ELEM(csize+iLocDom).bdySegXx)), ELEM(csize+iLocDom).boundX = test_bdy; end
    if(iBdy=='Y' && ~isempty(ELEM(csize+iLocDom).bdySegYx)), ELEM(csize+iLocDom).boundY = test_bdy; end
end

end

%%%% Function to get end points of the segment on the boundary
function [xint, yint] = get_bdy_seg_ep(iBdy,elem,Bdy,x,y,VERT,VERT_ELE,PARAMS)

if(~PARAMS.refinement)
    xe = x(elem.nodes(elem.nodes(1:3)<=PARAMS.num_bg_nod)); ye = y(elem.nodes(elem.nodes(1:3)<=PARAMS.num_bg_nod));
else
    
   inn = 1; 
   for in = 1:PARAMS.nvert
        if (elem.vertices_ele_virt(in)==1)
            xe(inn,1) = VERT_ELE(elem.vertices_ele(in)).X; ye(inn,1) = VERT_ELE(elem.vertices_ele(in)).Y;
            inn = inn+1;
        end
    end
end

if(~isempty(Bdy))
    if(iBdy=='X')
        BdyId = find(xe==Bdy);
    elseif(iBdy=='Y')
        BdyId = find(ye==Bdy);
    end
    
    if(isempty(elem.vertices))
        xint = xe(BdyId)';
        yint = ye(BdyId)';
    else
        yvert = [VERT(elem.vertices).Y]; xvert = [VERT(elem.vertices).X];
        if(iBdy=='X')
            VertBdyId = find(xvert==Bdy);
            BdyId = find(xe==Bdy);
        elseif(iBdy=='Y')
            VertBdyId = find(yvert==Bdy);
            BdyId = find(ye==Bdy);
        end
        xint = [xe(BdyId) xvert(VertBdyId)];
        yint = [ye(BdyId) yvert(VertBdyId)];
    end
else
    xint = []; yint = [];
end

if(~isempty(xint))
hold on;
plot(xint,yint);
end

end