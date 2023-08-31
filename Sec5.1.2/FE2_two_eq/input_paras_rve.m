%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% conductivity (W/mK)
coefs.kc0_be = 4e2;  % matrix beta phase
coefs.kc1_be = 4;

coefs.kc0_si = 1;    % inclusion sigma phase
coefs.kc1_si = 0.01;

% heat capacity Jm^3/K
coefs.hc0 = 9e3 * 3.9e2*5;  % matrix    beta phase
coefs.hc1 = 9e3 * 3.9e2*5;  % inclusion sigma phase

% % heat capacity Jm^3/K
% coefs.hc0 = 9e3 * 3.9e2*5;  % matrix    beta phase
% coefs.hc1 = 9e3 * 3.9e2*5;  % inclusion sigma phase

% packing parameters
% coefs.deltaT = dt;

paras.coefs = coefs;
  