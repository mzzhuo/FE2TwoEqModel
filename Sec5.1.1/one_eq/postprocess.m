% output
bc = bound.right;
nodesright = unique(bc);
conc_right = mean(u( dofArray(nodesright,2), :));
%%
% output
bc = bound.right;
nodesright = unique(bc);
poten_right = mean(u( dofArray(nodesright,1), :));

%% plot the concentration at last step
figure
hold on;
plot(coords(:,1),u(:,end),'o','MarkerSize',3);
set(gca, 'Fontname', 'Times New Roman','FontSize',20)
xlabel('$x$','Interpreter', 'latex','fontsize',24);
ylabel('temp. ($^\circ$)','Interpreter', 'latex','fontsize',24);
% axis([0 4 0 1]);
box on
% ---------------------------------------------------------------------------
%%
hold on;
plot(temp0(:,1),temp0(:,2),'r-','Linewidth',3);

%%
x = coords ((coords(:,2) == 0.0),1); 
y = u_macro(coords(:,2) == 0.0);

data = [x y];
[~,idx] = unique(data(:,1));
data = data(idx,:);

figure;
plot(x,y,'o')

% ---------------------------------------------------------------------------
%%
fileID = fopen('ourmacro.txt','w');
fprintf(fileID,'%12.6f, %12.6f \n',data');
fclose(fileID);
%%
hold on;
plot(matrixref(:,1),matrixref(:,2),'rs')
%%
hold on;
plot(inclusion(:,1),inclusion(:,2),'gs')
%%
hold on;
plot(midlineour(:,4),midlineour(:,1),'gs')
%%
hold on;
plot(mixref(:,1),mixref(:,2),'rs')
%%