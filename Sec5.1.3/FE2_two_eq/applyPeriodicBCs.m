%%  
% use Lagrange multiplier to enforce constraint 
% apply periodic boundary conditions 
% 
% 
K_peri = zeros(nlamb,nnode);
Res_peri = zeros(nlamb,1);
C_peri = zeros(nlamb,4);
% coefficient for macro variables downscaled
% [ grad_p(1) grad_p(2) phi grad_c(1) grad_c(2) c]'
% c_s - c_m - delta_c * (x_s - x_m) = 0
for i = 1:nlamb
  colm = dofArray(bound.master(i),1);
  cols = dofArray(bound.slave(i),1);
  K_peri(i,colm) = -1;
  K_peri(i,cols) = 1;
  dx = coords(bound.slave(i),:) - coords(bound.master(i),:);
  Res_peri(i) = u(cols) - u(colm) - dx * gradU_mat;
    
  C_peri(i,1) = - dx(1);
  C_peri(i,2) = - dx(2);
end
%% 
% apply periodic boundary conditions for corner nodes
% do it separately since it has 3 constraints
K_corn = zeros(3,nnode);
Res_corn = zeros(3,1);
C_corn = zeros(3,4);
% c_s - c_m - delta_c * (x_s - x_m) = 0
for i = 1:3
  colm = dofArray(bound.corner(1),1); % node 1
  cols = dofArray(bound.corner(i+1),1);
  K_corn(i,colm) = -1;
  K_corn(i,cols) = 1;
  dx = coords(bound.corner(i+1),:) - coords(bound.corner(1),:);
  Res_corn(i) = u(cols) - u(colm) - dx * gradU_mat;
    
  C_corn(i,1) = - dx(1);
  C_corn(i,2) = - dx(2);    
end
F_bond = - [K_peri; K_corn]' * u(nnode+1:nnode+nlamb+3);
%% 
% apply the conservation of heat energy
% 
elemtag = 0;
[K_consA,Res_consA,vol_mat] = ...
    applymasscons(ncoord,coords,nnode,nelnodes,connect.mat,elemtag,dofArray,u,u_mat,coefs);
% 
elemtag = 1;
[K_consB,Res_consB,vol_inc] = ...
    applymasscons(ncoord,coords,nnode,nelnodes,connect.inc,elemtag,dofArray,u,u_inc,coefs);

C_consA = zeros(1,4); 
C_consB = zeros(1,4);
C_consA(1,3) = coefs.hc0 * vol_mat; % this value may not be important.
C_consB(1,4) = coefs.hc1 * vol_inc; % this value may not be important.

K_cons = [K_consA; K_consB];
Res_cons = [Res_consA; Res_consB];
C_cons = [C_consA; C_consB];

% F_cons = - K_cons' * u(nnode+nlamb+3+1:nnode+nlamb+5);

% connect_all = [connect.mat; connect.inc];
% elemTag = [zeros(size(connect.mat,1),1); ones(size(connect.inc,1),1)];
% % 
% [K_cons,Res_cons,vol_mat,vol_inc] = ...
%     applymasscons_wholerve(ncoord,coords,nnode,nelnodes,connect_all,elemTag,dofArray,u,u_mat,u_inc,coefs);
% C_cons = zeros(1,4);
% C_cons(1,3) = coefs.hc0 * vol_mat;
% C_cons(1,4) = coefs.hc1 * vol_inc;
%%
% interface constraint
node_int = interface.nodes;
K_inte = zeros(nnodepair,nnode);
Res_inte = zeros(nnodepair,1);
C_inte = zeros(nnodepair,4);
for rw = 1:nnodepair
    cl_mat = node_int(rw,2); % dof of interface node of matrix
    K_inte(rw, cl_mat) = 1;  
    cl_inc = node_int(rw,3);
    K_inte(rw, cl_inc) = -1;
% 
    Res_inte(rw) = Res_inte(rw) + u(cl_mat) - u(cl_inc);
end
F_infc = - K_inte' * u(nnode+nlamb+5+1:end);
%%
% solve for the increment
% 
K1 = [K_peri; K_corn; K_cons; K_inte];
Res1 = [Res_peri; Res_corn; Res_cons; Res_inte];
C = [C_peri; C_corn; C_cons; C_inte];
% 
K = [K0, K1';  K1, zeros(nlamb+nnodepair+5,nlamb+nnodepair+5)];
% 
F_int = Res0;
% 
Res0 = Res0 + K1'*u(nnode+1:end);
Res = [Res0; Res1];

% try no interfacial condition
% K1 = [K_peri; K_corn; K_cons];
% Res1 = [Res_peri; Res_corn; Res_cons];
% C = [C_peri; C_corn; C_cons];
% K = [K0, K1';  K1, zeros(nlamb+5,nlamb+5)];
% F_int = Res0;
% Res0 = Res0 + K1'*u(nnode+1:end);
% Res = [Res0; Res1];

% interfacial condition enforced in assembling global matrix
% K1 = [K_peri; K_corn; K_cons];
% Res1 = [Res_peri; Res_corn; Res_cons];
% C = [C_peri; C_corn; C_cons];
% % 
% K1 = [K1, zeros(nlamb+5,nnodepair)];
% K = [K0, K1';  K1, zeros(nlamb+5,nlamb+5)];
% % F_int = Res0;
% Res0 = Res0 + K1'*u(nnode+nnodepair+1:end);
% Res = [Res0; Res1];