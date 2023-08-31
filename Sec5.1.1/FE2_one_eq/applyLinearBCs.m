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
% apply the conservation of c and phi
[K_cons, Res_cons,volume] = applymasscons(ncoord,ndof,nnode,coords,nelnodes,connect,dofArray,u,u_ma);
C_cons = zeros(1,3);
C_cons(1,3) = -volume;
%%
% solve for the increment
K1 = [K_lamb; K_cons];
Res1 = [Res_lamb; Res_cons];
C = [C_lamb; C_cons];

K = [K0, K1';  K1, zeros(nfix+1,nfix+1)];
F_int = Res0;
Res0 = Res0 + K1'*u(ndof*nnode+1:end);
Res = [Res0; Res1];
