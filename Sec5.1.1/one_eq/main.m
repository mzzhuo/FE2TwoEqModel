%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This codes solve diffusion equation
% fully resolved model and FE^2 
% 
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  
% 
clear;clc
% close all 

if ( ~exist('results', 'dir'))
   mkdir('results'); 
end

addpath('../utilities','./inputfiles','results');

% read gmsh input file

% infile = 'BatteryRef';
% [nnode,coords,nelem,connect,interface,nlamb,addLambda,bound]...
%     = read_input_file(infile);
% 
infile = 'fullorder';
elemType = 'Tri3';
% elemType = 'Quad4';
[nnode,coords,connect,elemTag,bound] = read_input_file_speT3(infile);
%% plot the mesh
% figure; hold on
% plot_mesh(coords,connect,'T3','k-')
% 
%% Assign dofs 
% each node has two dofs---potential and concentration
% an extra dof for each pair of interface nodes: lagrange multiplier

assignDofs; 
% 
%% Initialization
dt = 2;
n_time = 10; % total time = n_time * dt
timesteps = n_time + 1; % one more to consider the first step at t=0
% 
input_paras;
% 
% iteration steps and tolerance
NR.iter = 10;
NR.tol  = 1e-12;
%%
u  = zeros(ndofs, timesteps);
% 
u(:,1) = t0;  % initial c in domains
% u_pre = u(:,1);
% 
% time marching
% 
for itime = 2 : timesteps
% 
% Initialize solution at time step itime with previous step value
  if itime > 1
      u(:,itime) = u(:,itime-1);
%       u_pre      = u(:,itime-1);
  end
% 
% solution(concentration) increament of current time step
% difference bwt time steps: itime - (itime - 1)
  du = zeros(ndofs,1);
%   
% Newton-Raphson iterations
% 
  err2 = 1.;
  nit = 0;
  resd = 1.0;
%% 
  while ((resd > NR.tol) && (nit < NR.iter)) % Newton Raphson loop
% 
    nit = nit + 1; % iteration number +1
% 
%   Initialize Global Tangent Matrix K and Global RHS
    [K, Res] = stiff_residual(ncoord,ndofs,coords,nelnodes,...
                              connect,elemTag,dofArray,paras,u(:,itime),du);
    % F_int = Res(1:nnode);
%     
%   apply boundary conditions

%     applyBCs; % no heat flux

    correc = - K \ Res;
%
%   update difference between two timesteps
    du = du + correc;
%   update solution at current time step
    u(:,itime) = u(:,itime) + correc;
% 
%   check convergence   
    if nit == 1
       correc1 = correc;
    end
        
    resd  = norm( correc ) / norm( correc1 );
    err1  = norm(correc) / norm(u(:,itime));
    err2  = norm(Res)/ndofs;
    fprintf('step %d/%d, ite num %d, res dis %f, res rhs %f\n', ...
             itime-1, timesteps-1, nit, resd, err2);
% 
  end
%     u_pre = u(:,itime);
%   
% 
end
%%
%   output to paraview
vtfile =['results/', infile, '_oneeq'] ;
toParaview(coords,connect,elemType,vtfile,u(dofArray(1:nnode,1),end));

% save allresult.mat;
% 
% save solu.mat u;