%================= ELEMENT STIFFNESS MATRIX for SPE================================
%
function [kel, rel] = elstif_macro(ncoord,nelnodes,lmncoord,coefs,elu_mat,eldu_mat,elu_inc,eldu_inc)
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

%     [flux_mat,dFlux_mat] = micro_mat(gradU_mat,u_mat);

%     dFlux_mat = [-0.011658 0 0; 0 -0.011658 0]; 
%     flux_mat = dFlux_mat * [gradU_mat; u_mat];

    dFlux_mat = [-0.0533 0 0; 0 -0.0533 0];
    flux_mat = dFlux_mat * [gradU_mat; u_mat];
%     
    dFlux_inc = [0.0 0 0; 0 0.0 0]; %[-0.1064 0 0; 0 -0.1064 0]; 
    flux_inc = dFlux_inc * [gradU_inc; u_inc];
    
    [sour_mat,dSou_mat,sour_inc,dSou_inc]  = micro_flux(gradU_mat,gradU_inc,u_mat,u_inc);
    
 
%     dSou_mat = [0 0 0 0]; 
%     dSou_inc = [0 0 0 0]; 
%     dSou_mat = [0, 0, -6708.16762837825, 6708.16762837834];
%     dSou_inc = [0, 0, 6708.16762837825, -6708.16762837834];
    
%     sour_mat = dSou_mat * [gradU_mat; u_mat; u_inc];
%     sour_inc = dSou_inc * [gradU_mat; u_mat; u_inc];
    
%%  
    eps_mat = 0.62; eps_inc = 0.38;
    
    coefref = (0.62*coefs.hc0 + 0.38*coefs.hc1) / (2.106*coefs.k0);
    
    coef_mat = eps_mat * coefs.hc0/ coefref / coefs.deltaT;
    coef_inc = eps_inc * coefs.hc1/ coefref / coefs.deltaT;

%   For the matrix
%   Compute the element stiffness matrix 
    kel_mat = kel_mat + coef_mat * (N * N') * w(intpt)*dt ...
              - dNdx * dFlux_mat * [dNdx'; N'] * w(intpt)*dt ...
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
              - sour_inc * N *  w(intpt)*dt;   
% 

    kel_mi = kel_mi - N * (dSou_mat(4) * N') * w(intpt)*dt;
    kel_im = kel_im - N * (dSou_inc([1 2 3]) * [dNdx'; N']) * w(intpt)*dt;

    kel = [kel_mat, kel_mi; kel_im, kel_inc];     
    rel = [rel_mat; rel_inc];
    
end

% 

end