% output
bc = bound.right;
nodesright = unique(bc);
conc_right = mean(u( dofArray(nodesright,2), :));
%%
% output
bc = bound.right;
nodesright = unique(bc);
poten_right = mean(u( dofArray(nodesright,1), :));

% ---------------------------------------------------------------------------
%% plot the temp. at last step --- rve
figure;
hold on;
plot(coords(:,1),u(1:nnode,end),'<','MarkerSize',4);
plot(coords(:,1),u(nnode+1:end,end),'s','MarkerSize',4);
set(gca, 'Fontname', 'Times New Roman','FontSize',10)
xlabel('$x$','Interpreter', 'latex','fontsize',14);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',14);
% axis([0 10 -.05 1]);
box on;

%%
hold on;
plot(temp0(:,1),temp0(:,2),'r-','Linewidth',3);

%%
x = coords ((coords(:,2) == 0.0),1); 
y = u(coords(:,2) == 0.0,end);

data = [x y];
[~,idx] = unique(data(:,1));
data = data(idx,:);

figure;
plot(x,y,'o')

% ---------------------------------------------------------------------------
%%
macroT_mat = zeros(41,2);
macroT_mat(:,1) = data(:,1);
for i = 1:size(macroT_mat,1)
    dofs_x = dofArray( abs(coords(:,1)-macroT_mat(i,1))<1e-3, 1);
    macroT_mat(i,2) = mean(u(dofs_x,end));
end
macroT_mat(:,1) = macroT_mat(:,1) / 10;
%
fileID = fopen('quintard_macroT_mat_fe2_t2_notfixTsig.csv','w');
fprintf(fileID,'%12.6e, %12.6e \n', macroT_mat');
fclose(fileID);

%%
macroT_inc = zeros(41,2);
macroT_inc(:,1) = data(:,1);
for i = 1:size(macroT_inc,1)
    dofs_x = dofArray( abs(coords(:,1)-macroT_inc(i,1))<1e-3, 2);
    macroT_inc(i,2) = mean(u(dofs_x,end));
end
macroT_inc(:,1) = macroT_inc(:,1) / 10;
%
fileID = fopen('quintard_macroT_inc_fe2_t2_notfixTsig.csv','w');
fprintf(fileID,'%12.6f, %12.6f \n',macroT_inc');
fclose(fileID);