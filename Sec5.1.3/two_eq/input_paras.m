%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
t0 = 0.0;      % Initial temp in matrix and inclusion

% conductivity (W/mK)
coefs.k0 = 4e2;  % matrix beta phase
coefs.k1 = 1;    % inclusion sigma phase

% heat capacity Jm^3/K
coefs.hc0 = 9e3 * 3.9e2*5;  % matrix    beta phase
coefs.hc1 = 9e3 * 3.9e2*5;  % inclusion sigma phase

% packing parameters
coefs.deltaT = dt;

paras.coefs = coefs;
% 
% coefs.spe1 = pro*spheat / dt;
% coefs.spe5 = 1.363636363636364;
% coefs.spe6 = 0.036037170892710;


% 0.0345*(0.62*1202+0.38*1.7e6)/(2.106*0.026)
% 0.0345*(2.14e-2)^2*6.45e6/(2.106*0.026)