% 
clear;clc
close all 


grad_ma = [0 0];
u_ma    = 101.755109365399;

addpath('../utilities','./inputfiles');

infile = 'inputfiles/RVE_inc';
elemType = 'Tri3';

% read gmsh input file

% infile = 'inputfiles/RVE';
% elemType = 'Quad4';
[nnode,coords,connect,elemTag,bound] = read_RVE_mesh(infile,elemType);
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
u  = zeros(ndofs,1);
% initial guess value
u(dofsC) = u_ma; 
% 
%% 
%   
% Newton-Raphson iterations
% 
  err2 = 1.;
  nit = 0;
  resd = 1.0;
%   
%   K = zeros(ndofs+length(bound.master)+5,ndofs+length(bound.master)+5);
%   Res = zeros(ndofs+length(bound.master)+5,1);
%% 
while ((err2 > NR.tol) && (nit < NR.iter)) % Newton Raphson loop
% 
  nit = nit + 1; % iteration number +1

% Initialize Global Tangent Matrix K and Global RHS
  [K0, Res0] = stiff_resi_rve(ncoord,ndof,nnode,coords,nelnodes,...
             connect,elemTag,dofArray,paras,u);
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
  err1  = norm(correc) / norm(u);
  err2  = norm(Res)/ndofs;
  fprintf('Micro: ite num %d, res dis %f, res rhs %f\n', nit, resd, err2);
% 
end
%%
micro2macro;
%%
vtfile =['results/', 'RVE_inc_oneq'] ;
toParaview(coords,connect,elemType,vtfile,u(1:nnode,end));