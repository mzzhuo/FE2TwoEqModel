%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% conductivity (W/mK)
coefs.k0 = 0.026;  % matrix beta phase
% coefs.k1 = 0.5;    % inclusion sigma phase

% heat capacity Jm^3/K
coefs.hc0 = 1202;  % matrix    beta phase
% coefs.hc1 = 1.7e6; % inclusion sigma phase

% packing parameters
% coefs.deltaT = dt;

paras.coefs = coefs;
  