%================= ELEMENT STIFFNESS MATRIX for SPE================================
%
function [kel, rel] = elstif_macro(ncoord,nelnodes,lmncoord,coefs,elu,eldu,elu_pre)
%
%    
   kel = zeros(nelnodes,nelnodes); 
   rel = zeros(nelnodes,1);
%
%  Set up integration points && weights    
%   
   npoints = numberofintegrationpoints(ncoord,nelnodes);
   xilist = integrationpoints(ncoord,nelnodes,npoints);
   w = integrationweights(ncoord,nelnodes,npoints);
%
%  Loop over the integration points
%%
for intpt = 1:npoints

%     Compute shape functions && derivatives wrt local coords
%   
    xi = xilist(:,intpt);
%     
    N = shapefunctions(nelnodes,ncoord,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,xi);
%
%   Compute the jacobian matrix && its determinant
%
    dxdxi = lmncoord' * dNdxi;
% 
%       
    dxidx = inv(dxdxi);
    dt = det(dxdxi);

%   Convert shape function derivatives:derivatives wrt global coords
%   dNdxi(nelnodes,ncoord): dNdx[i,j] = dN(i)/dx(j); 
    dNdx = dNdxi * dxidx;
%       
%   calculate potential and concentration and their gradient 
%   at this gauss point
    % phi_grad = dNdx' * elphi;
    u_grad   = dNdx' * elu;
    % phi = N' * elphi;
    u   = N' * elu;
%   call the micro FE computation
    [r_homo,hc_homo,flux,stiff] = micro(u_grad',u);
%
    coef1 = hc_homo / coefs.deltaT;

%   Compute the element stiffness matrix 
    kel = kel + coef1 * (N * N') * w(intpt)*dt ...
              - dNdx * stiff * [dNdx'; N'] * w(intpt)*dt;
% 
%   compute the residual of c: r_c
    rel = rel +  coef1 * N * (N' * eldu) * w(intpt)*dt ...
              - dNdx * flux *  w(intpt)*dt ...
              - r_homo * N *  w(intpt)*dt;          
          
end

% 

end