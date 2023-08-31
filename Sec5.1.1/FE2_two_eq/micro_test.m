% 
clear;
% clc


gradU_mat = [10 10]';
gradU_inc = [0 0]';
% u_mat     = 98.4189870004200;
% u_inc     = 109.040360024100;
% u_mat     = 98.565446;
% u_inc     = 109.194927;
u_mat     = 14.0939813064586;
u_inc     = 22.8986929065973;

% u_mat     = 20;
% u_inc     = 20;

% u_mat     = 1.86443737744952; 
% u_inc     = 9.46432280764660;

addpath('../utilities','./inputfiles', './results');

infile = 'inputfiles/RVE_inc';
% infile = 'inputfiles/RVEnbyn';
% infile = 'inputfiles/rvesizestudy';
elemType = 'Tri3';

% read gmsh input file

% infile = 'inputfiles/RVE';
% elemType = 'Quad4';
[nnode,coords,connect,interface,bound] = read_RVE_inc_mesh(infile,elemType);
% 
%% Assign dofs 
% each node has two dofs---potential and concentration
% an extra dof for each pair of interface nodes: lagrange multiplier
bc = "periodic";
% bc = "linear";
assignDofs_rve; 
% 
% dt = 0.4;
% 
input_paras_rve;
% 
% input_paras;
% 
% iteration steps and tolerance
NR.iter = 10;
NR.tol  = 1e-12;
%%
% Initialize solution vectors
% nlamb = length(bound.master);
% bounodes = [bound.master; bound.slave; bound.corner];
% bounodes = unique(bounodes);
% fixnodes = bounodes(2:end);
% nfix = length(fixnodes);
% initial guess value
u = zeros(ndofs,1);
u(:,1) = u_mat; 
% 
%% 
%   
% Newton-Raphson iterations
% 
%   err2 = 1.;
  nit = 0;
  resd = 1.0;
%   
%   K = zeros(ndofs+length(bound.master)+5,ndofs+length(bound.master)+5);
%   Res = zeros(ndofs+length(bound.master)+5,1);
%% 
while ((resd > NR.tol) && (nit < NR.iter)) % Newton Raphson loop
% 
  nit = nit + 1; % iteration number +1

% Initialize Global Tangent Matrix K and Global RHS
  [K0, Res0] = stiff_resi_rve(ncoord,ndof,nnode,nnodepair,coords,nelnodes,...
             connect,interface,dofArray,paras,u);
%          
  applyBCs_rve;
% 
% solve for the increment
  correc = - K \ Res;

% update solution at current time step
  u = u + correc;

% check convergence   
  if nit == 1
     correc1 = correc;
  end
      
  resd = norm( correc ) / norm( correc1 );
%   err1  = norm(correc) / norm(u);
  err2  = norm(Res)/ndofs;
  fprintf('Micro: ite num %d, res dis %e, res rhs %e\n', nit, resd, err2);
% 
end
%%
micro2macro;
%%
%   output to paraview
vtfile =['./results/', 'RVE_inc_twoeq'] ;
toParaview(coords,connect,elemType,vtfile,u(dofArray(1:nnode,1),end));
