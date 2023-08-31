%%  
% use Lagrange multiplier to enforce constraint 
% apply periodic boundary conditions 
% 
K_peri = zeros(nlamb,nnode);
Res_peri = zeros(nlamb,1);
% coefficient for macro variables downscaled
C_peri = zeros(nlamb,3);

% [ grad_p(1) grad_p(2) phi grad_c(1) grad_c(2) c]'
% c_s - c_m - delta_c * (x_s - x_m) = 0
for i = 1:nlamb
  colm = dofArray(bound.master(i),1);
  cols = dofArray(bound.slave(i),1);
  K_peri(i,colm) = -1;
  K_peri(i,cols) = 1;
  dx = coords(bound.slave(i),:) - coords(bound.master(i),:);
  Res_peri(i) = u(cols) - u(colm) - grad_ma * dx';
    
  C_peri(i,1) = -dx(1);
  C_peri(i,2) = -dx(2);
end
%% 
% apply periodic boundary conditions for corner nodes
% do it separately since it has 3 constraints
K_corn = zeros(3,ndof*nnode);
Res_corn = zeros(3,1);
C_corn = zeros(3,3);
% c_s - c_m - delta_c * (x_s - x_m) = 0
for i = 1:3
  colm = dofArray(bound.corner(1),1); % node 1
  cols = dofArray(bound.corner(i+1),1);
  K_corn(i,colm) = -1;
  K_corn(i,cols) = 1;
  dx = coords(bound.corner(i+1),:) - coords(bound.corner(1),:);
  Res_corn(i) = u(cols) - u(colm) - grad_ma*dx';
    
  C_corn(i,1) = -dx(1);
  C_corn(i,2) = -dx(2);    
end
%% 
% apply the conservation of heat energy
% 
[K_cons,Res_cons,r_homo,hc_homo,volume] = applymasscons(ncoord,ndof,nnode,coords,nelnodes,connect,elemTag,dofArray,u,u_ma,coefs);
% 
C_cons = zeros(1,3);
C_cons(1,3) = volume * hc_homo; % this value may not be important because F_int/flux does not depend on u^M

F_cons = -K_cons'*u(ndof*nnode+nlamb+3+1:ndof*nnode+nlamb+4); 
%%
% solve for the increment
% K1 = [K_peri; K_corn];
% Res1 = [Res_peri; Res_corn];
% C = [C_peri; C_corn];
K1 = [K_peri; K_corn; K_cons];
Res1 = [Res_peri; Res_corn; Res_cons];
C = [C_peri; C_corn; C_cons];

K = [K0, K1';  K1, zeros(nlamb+4,nlamb+4)];
F_int = Res0; 
% F_int = -K1'*u(ndof*nnode+1:end); 
Res0 = Res0 + K1'*u(ndof*nnode+1:end);
Res = [Res0; Res1];

% F_int1 = - K_cons'*u(end); 
% F_int2 = - [K_peri; K_corn]'*u(ndof*nnode+1:end-1);

% fix u at node 1 to be 0;
% rw = dofArray(1,1);
% K(rw,:) = 0.0;
% K(rw,rw) = 1.0;
% Res(rw) = -(0 - u(rw));
