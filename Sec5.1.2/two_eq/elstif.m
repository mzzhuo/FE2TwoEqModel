%================= ELEMENT STIFFNESS MATRIX for SPE================================
%
function [kel, rel] = elstif(ncoord,nelnodes,lmncoord,coefs,elemtag,elu,eldu)
%
%  Assemble the element stiffness

   kel = zeros(nelnodes,nelnodes);

   rel = zeros(nelnodes,1);

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

    dxidx = inv(dxdxi);
    dt = det(dxdxi);

%   Convert shape function derivatives:derivatives wrt global coords
    dNdx = dNdxi * dxidx;

%   Compute the element stiffness matrix 

    if (elemtag == 0)
        ku = coefs.kc0_be + coefs.kc1_be*(N'*elu);
        kc1 = coefs.kc1_be;
        hc = coefs.hc0;
        r = 0; % no heat generation in the matrix
    elseif (elemtag == 1)
        ku = coefs.kc0_si + coefs.kc1_si*(N'*elu);
        kc1 = coefs.kc1_si;
        hc = coefs.hc1;
        r = 9e7;  % * 0.3318/0.3281; % no heat generation in the matrix
    else
       error("element type is wrong!"); 
    end

    coef1 = hc / coefs.deltaT;

%       coef1 = 1.0 / coefs.deltaT;
 
%     if (ku < 0)
%         error("negative ku");
%     end
%
    kel = kel + coef1 * (N * N') * w(intpt)*dt ...
                + ku * (dNdx*dNdx') *  w(intpt)*dt ...
                + dNdx * (dNdx'*elu) * kc1 * N' * w(intpt)*dt;

%   compute the residual of phi: r_c
    rel = rel +  coef1 * N * (N' * eldu) * w(intpt)*dt ...
              + ku * dNdx * (dNdx'*elu) *  w(intpt)*dt ...
              - r * N *  w(intpt)*dt;
end


end