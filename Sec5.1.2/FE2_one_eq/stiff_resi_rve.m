%====================== Assemble the global stiffness and residual =================
%
function [K, Res] = stiff_resi_rve(ncoord,ndof,nnode,coords,...
                                nelnodes,connect,elemTag,dofArray,paras,u)
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
  ndofs = ndof*nnode;
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
%
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
    elemtag = elemTag(lmn); % 0 - matrix; 1 - inclusion
    [kel, rel] = elstif_rve(ncoord,nelnodes,lmncoord,coefs,elemtag,elu); 
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

% 
end
