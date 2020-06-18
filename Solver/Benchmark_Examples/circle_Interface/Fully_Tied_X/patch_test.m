%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Planar interface with perfect continuity at the interface        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

%%%%%%%%% Add path for the source codes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath;
absolutPath = '/Users/wenjiang/Documents/MATLAB/Nitsche_2nd_order/devel/tri_ele/';
%absolutPath = '/Users/wenjiang/Documents/MATLAB/Nitsche_2nd_order/Yingjie Liu/tri_ele/';
addpath(genpath([absolutPath 'Solver/Source/']))
addpath(genpath([absolutPath 'Solver/Benchmark_Examples/circle_Interface/Fully_Tied_X/']))

%%%%%%%%% Structure PARAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gamFlag:      ( 0 -> gamma = 1/2 )
%               ( 1 -> gamma = "Hansbo" : Weighted area)
%               ( 2 -> gamma = "Robust": Weighted area and stiffnesses ) 
% dirbdy_n:     flag for enforcing normal bcs on external bdy
% dirbdy_t:     flag for enforcing tangential bcs on external boundary
% nit_bdy:      flag for using Nitsche's method on external boundary
% plane_strain: flag to determine whether plane stress or plane strain
%               conditions exist
% plotInt:      flag for exporting interfacial fields to vtk format
% export:       flag for generating vtk files
% ndof:         (2 -> vector problems)
%               (1 -> scalar problems)
% tol:          Used at multiple places as a cut-off for zero, criteria for
%               convergence etc.
% num_out_bdy:  4 -> when the external boundary is a boxed region. Used in
%               sub-routines for weak enforcement of external bcs.
%               Obviously not general enough for other geometries.
% LiuBorja:     Used to terminate a crack as specified in the example
%               investigated in Liu and Borja, 2008. False for all other
%               benchmark problems
% jump:         flag to check for a non-zero jump constraint. If true,
%               additional terms to the right hand-side need to be added in 
%               the Nitsche routines
% weakBC:       flag for weak enforcement of Dirichlet bcs on external bdy

PARAMS = struct('gamFlag', 2, 'nit_n', true, 'nit_t', true, 'dirBdy_n', true, ...
    'dirBdy_t', true,'nit_bdy', true,'plane_strain', true, 'plotInt', true,...
    'export', true, 'ndof', 1, 'maxIter', 50,'tol', 1e-16,...
    'num_out_bdy', 4,'LiuBorja', false, 'jump', false, 'tjump',true,'gradient', false,...
    'weakBC', true, 'debug', false,'refinement',false);

[E,nu,yieldParam] = mat_param(PARAMS);
PARAMS.pen_n = 1e5*max(E); PARAMS.pen_t = 1e5*max(E);
PARAMS.pen_bdy = 1e5*mean(E);
%%%%%%%% Background mesh and level set information %%%%%%%%%%%%%%%%%%%%%%%
% squaremesh and rectmesh are internal meshers.
% xdiv=19; ydiv=19;
% [x,y,node] = squaremesh(xdiv,ydiv);

global L D
L = 2.2;
D = 2.2;

ndiv(1) = 64;
ndiv(2) = 64;
[ncoord,eco] = createmesh(ndiv,0.0,'T6');

%plot_mesh(ncoord,eco,'T6','g.-',1.2);

x = ncoord(:,1);
y = ncoord(:,2);
node = eco';

PARAMS.refine_levels = 2;
if(PARAMS.refine_levels==0)
    PARAMS.refinement = false;
else
    PARAMS.refinement = true;
end


% Also supports using external meshes as long as the mesh file can be read
% and the x and y coordinates of nodes and nodal connectivity can be
% determined. See the two lines below for an example. readmesh is included
% in the source folder and supports reading files from Gmsh

% mesh_file = [absolutPath,'../Pre-PostProcessor/Mesh_Files/unstructSquare.msh'];
%[x,y,node] = readmesh(mesh_file);

%%%%%%%%% Level set information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.Xmin = min(x); PARAMS.Xmax = max(x); PARAMS.xint = 0.5;
PARAMS.Ymin = min(y); PARAMS.Ymax = max(y); PARAMS.yint = 0.5;

%ls = x - (PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin);

xc = PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin;

yc = PARAMS.yint*(PARAMS.Ymax-PARAMS.Ymin)+PARAMS.Ymin;

%ls = x.^2 + y.^2 - 1/4; 

ls1 = (x-0.25).^2+ (y-0.25).^2- 0.1^2;
ls2 = (x-0.75).^2+ (y-0.75).^2- 0.1^2;
ls3 = (x-0.25).^2+ (y-0.75).^2- 0.1^2;
ls4 = (x-0.75).^2+ (y-0.25).^2- 0.1^2;

ls5 = (x+0.25).^2+ (y-0.25).^2- 0.1^2;
ls6 = (x+0.75).^2+ (y-0.75).^2- 0.1^2;
ls7 = (x+0.25).^2+ (y-0.75).^2- 0.1^2;
ls8 = (x+0.75).^2+ (y-0.25).^2- 0.1^2;

ls9 = (x-0.25).^2+ (y+0.25).^2- 0.1^2;
ls10 = (x-0.75).^2+ (y+0.75).^2- 0.1^2;
ls11 = (x-0.25).^2+ (y+0.75).^2- 0.1^2;
ls12 = (x-0.75).^2+ (y+0.25).^2- 0.1^2;

ls13 = (x+0.25).^2+ (y+0.25).^2- 0.1^2;
ls14 = (x+0.75).^2+ (y+0.75).^2- 0.1^2;
ls15 = (x+0.25).^2+ (y+0.75).^2- 0.1^2;
ls16 = (x+0.75).^2+ (y+0.25).^2- 0.1^2;

ls = min([ls1 ls2 ls3 ls4 ls5 ls6 ls7 ls8 ls9 ls10 ls11 ls12 ls13 ls14 ls15 ls16],[],2);

%ls = x - xc - (y - yc);

%ls = x - y + 0.25;

% = 2*(x-0.25) - y;
%ls = 10*(x-0.2)-y;

%%%%%%%%% Preprocessing to account for Hansbo formulation %%%%%%%%%%%%%%%%%
%[ELEM,NODE,VERT,PARAMS] = comp_geo_ls(x,y,node,ls,PARAMS);
[ELEM,NODE,VERT,VERT_ELE,PARAMS] = comp_geo_ls(x,y,node,ls,PARAMS);

% PARAMS.DB and PARAMS.NB contain information about Dirichlet and Neumann
% boundary conditions. Rows of PARAMS.DB/PARAMS.NB include a binary flag 
% indicating whether a boundary is a Dirichlet/Neumann boundary. Columns
% include information about which boundary.
PARAMS.DB = [1 PARAMS.Xmin;1 PARAMS.Xmax;1 PARAMS.Ymin;1 PARAMS.Ymax]; 
PARAMS.NB = [0 PARAMS.Xmin;0 PARAMS.Xmax;0 PARAMS.Ymin;0 PARAMS.Ymax];
PARAMS.RB = [0 PARAMS.Xmin;0 PARAMS.Xmax;0 PARAMS.Ymin;0 PARAMS.Ymax];

 
%%%%%%%%% Directory path for plotting files %%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS.plotPath = [absolutPath 'Pre-PostProcessor/Results/circle_Interface/Fully_Tied/Fully_Tied_X/'];
PARAMS.plotFileBase = ['FullyTied_Nit_gamma_' num2str(PARAMS.gamFlag)];


%%%%%%%%% Mechanical solver routines belong here %%%%%%%%%%%%%%%%%%%%%%%%%%
%[u] = solve_nitsche(ELEM,NODE,VERT,PARAMS);
[u] = solve_nitsche(ELEM,NODE,VERT,VERT_ELE,PARAMS);
[L2Err] = get_L2_error_bulk_scalar(NODE,ELEM,VERT,VERT_ELE,PARAMS,u);

%generateVTK(ELEM,NODE,VERT,PARAMS,u);
generateVTK_refine(ELEM,NODE,VERT,VERT_ELE,PARAMS,u)
