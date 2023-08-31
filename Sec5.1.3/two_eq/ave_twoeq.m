% analyze fully resolved solution

% parameters
L = 1; 
H = 1/140; 
nx = 140;  % circle number in x direction
ny = 1;   % circle number in y direction
r_ca = L/nx * 0.3250;
%%
T_mat = zeros(nx,2); % first colum count of node, sec column potential; third column conc 
% 
node_mat = unique(connect.mat);

for i = 1:length(node_mat)
    xcoord = coords(node_mat(i),1);
    for j = 1:nx
        if (xcoord >= (j-1)*L/nx) && (xcoord <= j*L/nx)
            T_mat(j,1) = T_mat(j,1) + 1; 
            T_mat(j,2) = T_mat(j,2) + u(dofArray(node_mat(i),1),7); % 1 means temp 
        end
    end
end
T_mat(:,2) = T_mat(:,2)./T_mat(:,1); 
T_mat(:,1) = L/nx/2 + (0:nx-1)*L/nx;  

%%
figure;
hold on;
plot([0;T_mat(:,1)],[1;T_mat(:,2)],'b-','Linewidth',3,'MarkerSize',3)

%%
% homogenized conc/potential in inclusion
%
%
T_inc = zeros(nx,2); % first colum count of node, sec column potential; third column conc 
% 
node_inc = unique(connect.inc);

for i = 1:length(node_inc)
    xcoord = coords(node_inc(i),1);
    for j = 1:nx
        if (xcoord >= (j-1)*L/nx) && (xcoord <= j*L/nx)
            T_inc(j,1) = T_inc(j,1) + 1; 
            T_inc(j,2) = T_inc(j,2) + u(dofArray(node_inc(i),1),7); % 1 means temp 
        end
    end
end
T_inc(:,2) = T_inc(:,2)./T_inc(:,1); 
T_inc(:,1) = L/nx/2 + (0:nx-1)*L/nx;  
%%
% figure;
hold on;
plot(T_inc(:,1), T_inc(:,2),'r-','Linewidth',2,'MarkerSize',3)

%%
% hold on;
% x = 0.01:.01:1;
% plot(x,ones(1,length(x))*T_inc(1,2),'-','Linewidth',2,'MarkerSize',3)
% plot(x,ones(1,length(x))*T_mat(1,2),'-','Linewidth',2,'MarkerSize',3)
%% write to files
data = T_mat;
fileID = fopen('eg2macroT_mat_72s_new.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);

data = T_inc;
fileID = fopen('eg2macroT_inc_72s_new.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);
