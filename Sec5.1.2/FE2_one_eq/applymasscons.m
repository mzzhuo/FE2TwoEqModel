function [K,Res,r_homo,hc_homo,volume] = applymasscons(ncoord,ndof,nnode,coords,nelnodes,connect,elemTag,dofArray,u,u_ma,coefs)
% 
% this func is to build the constraint of conserved c and phi
% c_bar = integral of c over whole domain / volume
% 
% initialize stiffness and residual matrix
% ---------------------------------------------------------------
  ndofs = ndof*nnode; 
%   K   = sparse(2,ndofs);
%   Res = sparse(2,1);  
  K   = zeros(1,ndofs);
  Res = zeros(1,1);
  hc_homo = 0.0;
   r_homo = 0.0;
  volume = 0.0;
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
% 
%     elu  = zeros(ndof*nelnodes,1);
% % 
  elu = u(dofArray(connect(lmn,:),1));
% 
  if (elemTag(lmn) == 0)
      hc = coefs.hc0;
      r  = coefs.r0;
  elseif (elemTag(lmn) == 1)
      hc = coefs.hc1;
      r  = coefs.r1;
  else
       error("element type is wrong!"); 
  end
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
    hc_homo = hc_homo + hc * w(intpt)*dt;
    r_homo  =  r_homo + r  * w(intpt)*dt;
    volume  = volume + w(intpt)*dt;

  end
% 
% put elem value to global
  cl = dofArray(connect(lmn,:),1);
  K(1,cl) = K(1,cl) + kel;
%   
  Res = Res + rel;
end
% 
hc_homo = hc_homo / volume;
 r_homo =  r_homo / volume;
Res = Res - volume * hc_homo * u_ma;

end
