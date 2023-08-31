%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% conductivity (W/mK)
coefs.k0 = 4e3;  % matrix beta phase
coefs.k1 = 1e-1;    % inclusion sigma phase

% % heat capacity Jm^3/K
% coefs.hc0 = 9e3 * 3.9e2;  % matrix    beta phase
% coefs.hc1 = 9 * 3.9e2;  % inclusion sigma phase

% heat capacity Jm^3/K
coefs.hc0 = 9e3 * 3.9e2*5;  % matrix    beta phase
coefs.hc1 = 9e3 * 3.9e2*5;  % inclusion sigma phase

% packing parameters
% coefs.deltaT = dt;

paras.coefs = coefs;
  