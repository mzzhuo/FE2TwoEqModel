%====================== Assemble the global stiffness and residual =================
%
function [K, Res] = stiff_resi_macro(ncoord,ndofs,coords,...
                                nelnodes,connect,dofArray,paras,u_pre,u,du)
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
    % eldu = zeros(ndof*nelnodes,1);
% 
    elu   =  u(dofArray(connect(lmn,:),1));
    eldu  = du(dofArray(connect(lmn,:),1));
    elu_pre = u_pre(dofArray(connect(lmn,:),1));

%   element stiffness and residual
    [kel, rel] = elstif_macro(ncoord,nelnodes,lmncoord,coefs,elu,eldu,elu_pre); 
%
%   Add the current element stiffness to the global stiffness
% 
    rows = dofArray(connect(lmn,:),1);
    cols = rows;

    K(rows,cols) = K(rows,cols) + kel;

    Res(rows) = Res(rows) + rel;
% 
end  

% 
end
