%%  
% use Lagrange multiplier to enforce constraint 
% apply fixed boundary conditions 
K_lamb = zeros(nfix,ndof*nnode);
Res_lamb = zeros(nfix,1);
C_lamb = zeros(nfix,3);
% 
for i = 1:nfix
  %
  cols = dofArray(fixnodes(i),1);
  colm = dofArray(1,1);
  K_lamb(i,cols) = 1;
  K_lamb(i,colm) = -1;
  dx = coords(fixnodes(i),:) - coords(1,:);
  Res_lamb(i) = u(cols) - u(colm) - grad_ma*dx';
  C_lamb(i,1) = -dx(1);
  C_lamb(i,2) = -dx(2);  
end

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

C_consA = zeros(1,4); C_consB = zeros(1,4);
C_consA(1,3) = coefs.hc0 * vol_mat; % this value may not be important.
C_consB(1,4) = coefs.hc1 * vol_inc; % this value may not be important.

K_cons = [K_consA; K_consB];
Res_cons = [Res_consA; Res_consB];
C_cons = [C_consA; C_consB];



%%
% solve for the increment
K1 = [K_lamb; K_cons];
Res1 = [Res_lamb; Res_cons];
C = [C_lamb; C_cons];

K = [K0, K1';  K1, zeros(nfix+1,nfix+1)];
F_int = Res0;
Res0 = Res0 + K1'*u(ndof*nnode+1:end);
Res = [Res0; Res1];
