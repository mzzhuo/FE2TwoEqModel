%================= ELEMENT STIFFNESS MATRIX for SPE================================
% 
function [kel, rel] = elstif_rve(ncoord,nelnodes,lmncoord,coefs,elemtag,elu)
%
%  Assemble the element stiffness
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
%   for i = 1:ncoord
%     xi(i) = xilist(i,intpt);
%   end    
    N = shapefunctions(nelnodes,ncoord,xi);
    dNdxi = shapefunctionderivs(nelnodes,ncoord,xi);

      
%
%   Compute the jacobian matrix && its determinant
%
    dxdxi = lmncoord' * dNdxi;
%     for i = 1:ncoord
%       for j = 1:ncoord
%         dxdxi(i,j) = 0.;
%         for a = 1:nelnodes
%           dxdxi(i,j) = dxdxi(i,j) + lmncoord(a,i)*dNdxi(a,j);
%         end
%       end
%     end
%       
    dxidx = inv(dxdxi);
    dt = det(dxdxi);

%   Convert shape function derivatives:derivatives wrt global coords
    dNdx = dNdxi * dxidx;
 
    if (elemtag == 0)
        ku = coefs.k0;
%         hc = coefs.hc0;
    elseif (elemtag == 1)
        ku = coefs.k1;
%         hc = coefs.hc1;
    else
       error("element type is wrong!"); 
    end

%     if (ku < 0)
%         error("negative ku");
%     end

    kel = kel + ku * (dNdx*dNdx') *  w(intpt)*dt;
    
%   compute the residual of phi: r_c
    rel = rel + ku * dNdx * (dNdx'*elu) *  w(intpt)*dt;

end
% 
end