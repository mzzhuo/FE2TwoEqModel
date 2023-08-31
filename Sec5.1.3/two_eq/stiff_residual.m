%====================== Assemble the global stiffness and residual =================
%
function [F_int,K, Res] = stiff_residual(ncoord,ndofs,coords,nnode,...
                                nelnodes,connect_0,interface,dofArray,paras,u,du)
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
  K   = sparse(ndofs,ndofs);
  Res = sparse(ndofs,1);  
%   K   = zeros(ndofs,ndofs);
%   Res = zeros(ndofs,1);
%% 
% ---------------------------------------------------------------
% Assemble the stiffness matrix and residual 
% ---------------------------------------------------------------
%
%   Loop over all the elements
%   
% 0 - matrix; 1 - inclusion
  elemTag = [zeros(size(connect_0.mat,1),1); ones(size(connect_0.inc,1),1)];
  connect = [ connect_0.mat; connect_0.inc ];
    
for lmn = 1:size(connect,1)
%
%   Extract coords of nodes, DOF for the current element
    lmncoord = coords(connect(lmn,:),:);
% 
    % elu  = zeros(ndof*nelnodes,1);
    % eldu = zeros(ndof*nelnodes,1);
% 
    elu   =  u(dofArray(connect(lmn,:),1));
    eldu  = du(dofArray(connect(lmn,:),1));
% 
%   element stiffness and residual
    elemtag = elemTag(lmn); % 0 - matrix; 1 - inclusion
    [kel, rel] = elstif(ncoord,nelnodes,lmncoord,coefs,elemtag,elu,eldu); 
%
%   Add the current element stiffness to the global stiffness
% 
    rows = dofArray(connect(lmn,:),1);
    cols = rows;

    K(rows,cols) = K(rows,cols) + kel;

    Res(rows) = Res(rows) + rel;

end  
%%
% ---------------------------------------------------------------
% Apply Lagrangian multiplier to the stiffness matrix and residual
% ---------------------------------------------------------------

F_int = Res;

% interface node numbering: 1,2,3...
% int nodes; corresponding spe nodes; electrode nodes
node_int = interface.nodes;

for i = 1:size(node_int,1)
% left bottom
    rw = dofArray(nnode+i,1);
    cl_mat = dofArray(node_int(i,2),1); % dof of interface node of matrix
    K(rw, cl_mat) = 1;  
    cl_inc = dofArray(node_int(i,3),1);
    K(rw, cl_inc) = -1;
% 
    Res(rw) = Res(rw) + u(cl_mat) - u(cl_inc);
%
% right upper
%     cl = dofArray(nnode+i,1); % cl = rw
%     rw_mat = dofArray(node_int(i,2),1); % dof of interface node of matrix
%     K(rw_mat,cl) = 1;
%     rw_inc = dofArray(node_int(i,3),1);
%     K(rw_inc,cl) = -1;
% 
    K(cl_mat,rw) = 1;
    Res(cl_mat) = Res(cl_mat) + u(rw);
%     
    K(cl_inc,rw) = -1; 
    Res(cl_inc) = Res(cl_inc) - u(rw);
end
% 
end
