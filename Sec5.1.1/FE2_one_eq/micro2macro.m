%% heat capacity
% in boundary conditions
%% flux
% 
% calculate homegenized 
% F_int = K1'*u(ndof*nnode+1:end);
cols = dofArray(bounodes,1);
D_b = coords(bounodes,:);
flux = - D_b'*F_int(cols)/volume;

% flux = - coords'*F_int/1e-4;

%% consistent matrix
delta = eye(ndofs,ndofs);
delta2 = delta(ndof*nnode+1:end,:);
rhs = [zeros(ndof*nnode,3); -C];
solu = K \ rhs;
% 
Stiff_tem = K1'*delta2*solu;
D_b  = coords(bounodes,:);
cols = dofArray(bounodes,1);
stiff = D_b'*Stiff_tem(cols,:)/volume; 
%% the following means constraint force of mass conservation do not matter
% delta = eye(ndofs,ndofs);
% delta2 = delta(ndof*nnode+1:end-1,:);
% Stiff_tem = [K_peri; K_corn]'*delta2*solu;
% D_b  = coords(bounodes,:);
% cols = dofArray(bounodes,1);
% stiff = D_b'*Stiff_tem(cols,:)/volume; 


% D = coords';
% Stiff_p = D*Stiff(1:nnode,:)/volume;
% Stiff_c = D*Stiff(nnode+1:end,:)/volume;
% % 

% %% check internal nodes 
% allnodes = unique(connect);
% allnodes(bounodes) = [];
% internodes = allnodes;
% internodes = [internodes; internodes+nnode];
% Stiff(internodes,:);
% check = F_int(internodes,:);
% end