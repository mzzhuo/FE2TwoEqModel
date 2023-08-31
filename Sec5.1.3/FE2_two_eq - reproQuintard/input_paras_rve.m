%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% conductivity (W/mK)
coefs.k0 = 1.0;  % matrix beta phase
coefs.k1 = 1e-2;    % inclusion sigma phase

% heat capacity Jm^3/K
coefs.hc0 = 1.0e6;  % matrix    beta phase
coefs.hc1 = 1.0e6;  % inclusion sigma phase

% packing parameters
% coefs.deltaT = dt;

paras.coefs = coefs;
  