function [dFlux_mat,flux_mat,dFlux_inc,flux_inc,sour_mat,dSou_mat,sour_inc,dSou_inc] = micro2(gradU_mat,gradU_inc,u_mat,u_inc)
% 
% addpath('../utilities','./inputfiles');
% 
% gradU_mat = [2 0]';
% gradU_inc = [0 0]';
% u_mat     = 1;
% u_inc     = 2;

infile = 'inputfiles/RVE_inc';
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
%   fprintf('Micro: ite num %d, res dis %f, res rhs %f\n', nit, resd, err2);
% 
end
%%
% calculate homogenized flux and its tangent
% flux with interface nodes
volume = vol_mat + vol_inc;
node_int = interface.nodes;

bounodes_mat = [bounodes; node_int(:,2)];
F_bi = F_bond + F_infc;

cols_mat = dofArray(bounodes_mat,1);
D_b = coords(bounodes_mat,:);
flux_mat = - D_b' * F_bi(cols_mat)/volume;

% consistent matrix
delta = eye(ndofs,ndofs);
delta = delta([nnode+1:nnode+nlamb+3 nnode+nlamb+5+1:end],:);
rhs = [zeros(ndof*nnode,4); -C];
solu = K \ rhs;
% stiff = K1'*delta*solu;
stiff = [K_peri; K_corn; K_inte]'*delta*solu;
% 
dFlux_mat = D_b'*stiff(cols_mat,:)/volume; 


bounodes_inc = node_int(:,3);
cols_inc = dofArray(bounodes_inc,1);
D_b = coords(bounodes_inc,:);
flux_inc = - D_b' * F_bi(cols_inc)/volume;

delta = eye(ndofs,ndofs);
delta = delta([nnode+1:nnode+nlamb+3 nnode+nlamb+5+1:end],:);
rhs = [zeros(ndof*nnode,4); -C];
solu = K \ rhs;
stiff = [K_peri; K_corn; K_inte]'*delta*solu;
% 
dFlux_inc = D_b'*stiff(cols_inc,:)/volume; 

%% source terms and their tangent matrix
% 
delta = eye(ndofs,ndofs);
delta = delta(nnode+nlamb+5+1:end,:);
rhs = [zeros(ndof*nnode,4); -C];
solu = K \ rhs;
stiff = K_inte'*delta*solu;
%
sour_mat = sum(F_infc(node_int(:,2),:))/volume; 
dSou_mat = sum(stiff(node_int(:,2),:))/volume; 
% 
sour_inc = sum(F_infc(node_int(:,3),:))/volume; 
dSou_inc = sum(stiff(node_int(:,3),:))/volume; 

end

%%
%   output to paraview
% vtfile =['./results/', 'RVE_hole'] ;
% toParaview(coords,connect,elemType,vtfile,u(dofArray(1:nnode,1),end),F_int);