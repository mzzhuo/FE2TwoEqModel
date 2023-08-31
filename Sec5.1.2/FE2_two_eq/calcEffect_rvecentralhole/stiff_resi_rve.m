%====================== Assemble the global stiffness and residual =================
%
function [K, Res] = stiff_resi_rve(ncoord,ndof,nnode,coords,...
                                nelnodes,connect,dofArray,paras,u)
%
% ---------------------------------------------------------------
% get the coefficients for equations
% ---------------------------------------------------------------
%
%   ncoord = 2; nelnodes = 3; 
%   ndof = 2; % each node has 2 dofs: potential and concentration
% 
% ---------------------------------------------------------------
% initialize stiffness and residual matrix
% ---------------------------------------------------------------
  coefs = paras.coefs;
% no. of dofs for all nodes
  ndofs = ndof*nnode;% + nnodepair;
%   K   = sparse(ndofs,ndofs);
%   Res = sparse(ndofs,1);  
  K   = zeros(ndofs,ndofs);
  Res = zeros(ndofs,1);
%%  
%   u = zeros(ndofs);  
%   kel = zeros(ndof*nelnodes,ndof*nelnodes); 
%   rel = zeros(ndof*nelnodes); 
% 
% ---------------------------------------------------------------
% Assemble the stiffness matrix and residual for SPE
% ---------------------------------------------------------------
%
%   Loop over all the elements

for lmn = 1:size(connect,1)
%
%   Extract coords of nodes, DOF for the current element
    lmncoord = coords(connect(lmn,:),:);

% 
    % elu  = zeros(ndof*nelnodes,1);
% 
    elu = u(dofArray(connect(lmn,:),1));
% 
%   element stiffness and residual
    [kel, rel] = elstif_rve(ncoord,nelnodes,lmncoord,coefs,elu); 
%
%   Add the current element stiffness to the global stiffness
% 
%   For easy understanding: 
%   kpp = kel(1:nelnodes,1:nelnodes);
%   kpc = kel(1:nelnodes,nelnodes+1:end);
%   kcc = kel(nelnodes+1:end,nelnodes+1:end);
% 
    rows = dofArray(connect(lmn,:),1);
    cols = rows;

    K(rows,cols) = K(rows,cols) + kel;

    Res(rows) = Res(rows) + rel;
% 
end  
% F_int = Res;
% ---------------------------------------------------------------
% Apply Lagrangian multiplier to the stiffness matrix and residual
% ---------------------------------------------------------------

% interface node numbering: 1,2,3...
% int nodes; corresponding spe nodes; electrode nodes
% node_int = interface.nodes;
% 
% for i = 1:size(node_int,1)
% % left bottom
%     rw = nnode+i;
%     cl_mat = node_int(i,2); % dof of interface node of matrix
%     K(rw, cl_mat) = 1;  
%     cl_inc = node_int(i,3);
%     K(rw, cl_inc) = -1;
% % 
%     Res(rw) = Res(rw) + u(cl_mat) - u(cl_inc);
% %
% % right upper
% % 
%     K(cl_mat,rw) = 1;
%     Res(cl_mat) = Res(cl_mat) + u(rw);
% %     
%     K(cl_inc,rw) = -1; 
%     Res(cl_inc) = Res(cl_inc) - u(rw);
% end
% 
end
