%% heat capacity
% in boundary conditions
volume = 1e-4;
%% flux
% calculate homegenized 
cols = dofArray(bounodes,1);
D_b = coords(bounodes,:);
flux = - D_b' * F_int(cols)/volume;


%% consistent matrix
delta = eye(ndofs,ndofs);
delta = delta(ndof*nnode+1:end,:);
rhs = [zeros(ndof*nnode,3); -C];
solu = K \ rhs;
stiff = K1'*delta*solu;
% 
D_b  = coords(bounodes,:);
cols = dofArray(bounodes,1);
stiff = D_b'*stiff(cols,:)/volume; 
