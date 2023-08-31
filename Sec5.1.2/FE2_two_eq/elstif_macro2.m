%================= ELEMENT STIFFNESS MATRIX for SPE================================
%
function [kel, rel] = elstif_macro2(ncoord,nelnodes,lmncoord,coefs,elu_mat,eldu_mat,elu_inc,eldu_inc)
%
%    
   kel_mat = zeros(nelnodes,nelnodes); 
   rel_mat = zeros(nelnodes,1);
   kel_inc = zeros(nelnodes,nelnodes); 
   rel_inc = zeros(nelnodes,1);
   kel_mi = zeros(nelnodes,nelnodes); 
   kel_im = zeros(nelnodes,nelnodes); 
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
    gradU_mat   = dNdx' * elu_mat; % column vector
    gradU_inc   = dNdx' * elu_inc; % column vector
    
    u_mat       = N' * elu_mat;    % column vector
    u_inc       = N' * elu_inc;    % column vector
    
%   call the micro FE computation

%     dFlux_mat = [-0.011658 0 0; 0 -0.011658 0]; 
%     flux_mat = dFlux_mat * [gradU_mat; u_mat];

    % dFlux_mat = [-0.0533 0 3.041 -3.041; 0 -0.0533 3.041 -3.041];
    % flux_mat = dFlux_mat * [gradU_mat; u_mat];
    
    dFlux_inc = [0.0 0 0; 0 0.0 0]; %[-0.1064 0 0; 0 -0.1064 0]; 
    flux_inc = dFlux_inc * [gradU_inc; u_inc];

    [dFlux_mat, flux_mat, sour_mat, dSou_mat,sour_inc,dSou_inc]  = micro(gradU_mat,gradU_inc,u_mat,u_inc);
%     sour_mat = 0;
%     dSou_mat = [0 0 0 0]; 
%     sour_inc = 0; 
%     dSou_inc = [0 0 0 0]; 
%     

%%  
    eps_mat = 0.6682; eps_inc = 0.3318;
    r_inc = 9e7 * eps_inc;
    
    coef_mat = eps_mat * coefs.hc0 / coefs.deltaT;
    coef_inc = eps_inc * coefs.hc1 / coefs.deltaT;

%   For the matrix
%   Compute the element stiffness matrix 
    kel_mat = kel_mat + coef_mat * (N * N') * w(intpt)*dt ...
              - dNdx * dFlux_mat(:,1:3) * [dNdx'; N'] * w(intpt)*dt ...
              - N * (dSou_mat([1 2 3]) * [dNdx'; N']) * w(intpt)*dt;
% 
%   compute the residual 
    rel_mat = rel_mat +  coef_mat * N * (N' * eldu_mat) * w(intpt)*dt ...
              - dNdx * flux_mat *  w(intpt)*dt ...
              - sour_mat * N *  w(intpt)*dt;          
          
%   For the inclusion
%   Compute the element stiffness matrix 
    kel_inc = kel_inc + coef_inc * (N * N') * w(intpt)*dt ...
              - dNdx * dFlux_inc * [dNdx'; N'] * w(intpt)*dt ...
              - N * (dSou_inc(4) * N') *  w(intpt)*dt;
% 
%   compute the residual 
    rel_inc = rel_inc +  coef_inc * N * (N' * eldu_inc) * w(intpt)*dt ...
              - dNdx * flux_inc *  w(intpt)*dt ...
              - (r_inc + sour_inc) * N *  w(intpt)*dt;  
% 

    kel_mi = kel_mi - N * (dSou_mat(4) * N') * w(intpt)*dt ...
             - dNdx * dFlux_mat(:,4) * N' * w(intpt)*dt;
    kel_im = kel_im - N * (dSou_inc([1 2 3]) * [dNdx'; N']) * w(intpt)*dt;

    kel = [kel_mat, kel_mi; kel_im, kel_inc];     
    rel = [rel_mat; rel_inc];
    
end

% 

end