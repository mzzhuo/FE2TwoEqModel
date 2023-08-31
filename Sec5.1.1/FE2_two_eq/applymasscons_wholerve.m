function [K,Res,vol_mat,vol_inc] = ...
         applymasscons_wholerve(ncoord,coords,nnode,nelnodes,connect,elemTag,dofArray,u,u_mat,u_inc,coefs)
% 
% this func is to build the constraint of conserved c and phi
% c_bar = integral of c over whole domain / volume
% 
% initialize stiffness and residual matrix
% ---------------------------------------------------------------
  ndofs = nnode; 
%   K   = sparse(2,ndofs);
%   Res = sparse(2,1);  
  K   = zeros(1,ndofs);
  Res = zeros(1,1);
  vol_mat = 0.0;
  vol_inc = 0.0;
%%  
%   u = zeros(ndofs);  
%   kel = zeros(ndof*nelnodes,ndof*nelnodes); 
%   rel = zeros(ndof*nelnodes); 
% 
% ---------------------------------------------------------------
% Assemble the stiffness matrix and residual for SPE
% ---------------------------------------------------------------
%
%  Set up integration points && weights    
%   
  npoints = numberofintegrationpoints(ncoord,nelnodes);
  xilist = integrationpoints(ncoord,nelnodes,npoints);
  w = integrationweights(ncoord,nelnodes,npoints);
% 
%   Loop over all the elements
%
for lmn = 1:size(connect,1)
%
  kel = zeros(1,nelnodes);
  rel = zeros(1,1);
% 
% Extract coords of nodes, DOF for the current element
  lmncoord = coords(connect(lmn,:),:);
% 

  if (elemTag(lmn) == 0)
      hc = coefs.hc0;
      ismat = 1; isinc = 0;
  elseif (elemTag(lmn) == 1)
      hc = coefs.hc1;
      ismat = 0; isinc = 1;
  else
       error("element type is wrong!"); 
  end

  elu = u(dofArray(connect(lmn,:),1));
%   
  for intpt = 1:npoints
%
%   Compute shape functions && derivatives wrt local coords
%   
    xi = xilist(:,intpt);
  
    N = shapefunctions(nelnodes,ncoord,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,xi);
%
%   Compute the jacobian matrix && its determinant
%    
    dxdxi = lmncoord' * dNdxi;
%       
    dxidx = inv(dxdxi);
    dt = det(dxdxi);
%   Convert shape function derivatives:derivatives wrt global coords
    dNdx = dNdxi * dxidx;    
% 
    for a = 1:nelnodes
      kel(a) = kel(a) + hc * N(a) * w(intpt)*dt;
      rel(1) = rel(1) + hc * N(a)*elu(a) * w(intpt)*dt;
    end
    vol_mat  = vol_mat + ismat * w(intpt)*dt;
    vol_inc  = vol_inc + isinc * w(intpt)*dt;
  end
% 
% put elem value to global
  cl = dofArray(connect(lmn,:),1);
  K(1,cl) = K(1,cl) + kel;
%   
  Res = Res + rel;
end
% 
Res = Res - (coefs.hc0 * vol_mat * u_mat + coefs.hc1 * vol_inc * u_inc);

% volume = vol_mat + vol_inc;

end
