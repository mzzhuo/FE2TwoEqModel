%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% conductivity (W/mK)
coefs.k0 = 1;  % matrix beta phase
coefs.k1 = 4e2;    % inclusion sigma phase

% heat capacity Jm^3/K
coefs.hc0 = 9e3 * 3.9e2;  % matrix    beta phase
coefs.hc1 = 9e3 * 3.9e2;  % inclusion sigma phase

% heat source density
coefs.r0 = 0;  % matrix    beta phase
coefs.r1 = 9e7;  % inclusion sigma phase

% 9e7*0.3318*20/coefs.hc0
% packing parameters
% coefs.deltaT = dt;

paras.coefs = coefs;
