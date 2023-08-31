% analyze fully resolved solution

% parameters
L = 0.00357143; 
H = 0.00357143; 
nx = 1;  % circle number in x direction
ny = 1;   % circle number in y direction
r_ca = L/nx * 0.3250;


T_homo = zeros(timesteps,2); % first colum time, sec column temp.
% 
node = unique(connect);

for j = 1:timesteps
    count = 0;
    for i = 1:length(node)
        xcoord = coords(node(i),1);
        if (xcoord >= 0) && (xcoord <= L/nx)
            count = count + 1;
            T_homo(j,2) = T_homo(j,2) + u(dofArray(node(i),1),j); % 1 means temp 
        end
    end
    T_homo(j,2) = T_homo(j,2)./count; 
end
T_homo(:,1) = (0:n_time)*dt;

%%
figure;
hold on;
plot(T_homo(:,1),T_homo(:,2),'o','MarkerSize',5)
%% paraview data
data = T_homo;
fileID = fopen('aveTemp_oneq.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);