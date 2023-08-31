clear;
l = 0.00357143;  % edge length of a unit cell
r = l * 0.325; % inclusion radius
n_x = 1;      % number of unit cells/inclusion
n_y = 1;
L = l * n_x;  % edge length of RVE in x dir
H = l * n_y;  % edge length of RVE in y dir
n = n_x * n_y; % number of inclusions

Dist2Voids = 2.1 * r;
% 
offset = 1.08; % distance of point from border in terms of radius .... dire a pietro

i = 0; % initial number of anode particles
coords = zeros(n,2);

while i < n
% 
  i = i + 1;
% 
  u = rand (2,1);
  x = offset * r + (L - 2 * offset * r) * u(1);
  y = offset * r + (H - 2 * offset * r) * u(2);
% 
  skipVoid = 0;
%   
% check if voids intersect
  for j = 1:1:i-1
      if (j == 0)
          continue
      end
      
      % j is any previous particle center
      
      distance = norm([x y] - coords(j,:));
      
      if distance < Dist2Voids
%           disp(' - Skipping overlapping voids');
          i = i - 1;
          skipVoid = 1;
          break
      end
  end

  if (skipVoid == 0)
    coords(i,1) = x;
    coords(i,2) = y;
  end 
end 
%%
% figure 
% scatter(coords(:,1),coords(:,2),10);
% axis([0 L 0 H])
% box on
%% write geo file
% write geo file
oFile = fopen('rvesizestudy.geo','w');
% 
lcar1 = L/n_x/10;
lcar2 = L/n_x/10;
% write the external frame
fprintf(oFile,'\n\n');
fprintf(oFile,'Point(1) = {%.8e,%.8e,0.0,%.8e};\n',0.0,0.0,lcar1);
fprintf(oFile,'Point(2) = {%.8e,%.8e,0.0,%.8e};\n',L,0.0,lcar1);
fprintf(oFile,'Point(3) = {%.8e,%.8e,0.0,%.8e};\n',L,H,lcar1);
fprintf(oFile,'Point(4) = {%.8e,%.8e,0.0,%.8e};\n',0.0,H,lcar1);
% 
fprintf(oFile,'Line(1) = {1,2};\n');
fprintf(oFile,'Line(2) = {2,3};\n');
fprintf(oFile,'Line(3) = {4,3};\n');
fprintf(oFile,'Line(4) = {1,4};\n');

% write circles central points coords for cathode
fprintf(oFile,'\n\n');
for i = 1:size(coords,1)
  fprintf(oFile,'coordx[%6d] = %.8e; coordy[%6d] = %.8e;\n',...
                i,coords(i,1),i,coords(i,2));
end

% write Macro of circle
fprintf(oFile,'\n\n');
fprintf(oFile,'Macro CheeseHole\n');
fprintf(oFile,'p1 = newp; Point(p1) = {x,   y,   0, %.8e};\n',lcar2);
fprintf(oFile,'p2 = newp; Point(p2) = {x+r, y,   0, %.8e};\n',lcar2);
fprintf(oFile,'p3 = newp; Point(p3) = {x,   y+r, 0, %.8e};\n',lcar2);
fprintf(oFile,'p4 = newp; Point(p4) = {x-r, y,   0, %.8e};\n',lcar2);
fprintf(oFile,'p5 = newp; Point(p5) = {x,   y-r, 0, %.8e};\n',lcar2);

fprintf(oFile,'c1 = newreg; Circle(c1) = {p2,p1,p3};\n');
fprintf(oFile,'c2 = newreg; Circle(c2) = {p3,p1,p4};\n');
fprintf(oFile,'c3 = newreg; Circle(c3) = {p4,p1,p5};\n');
fprintf(oFile,'c4 = newreg; Circle(c4) = {p5,p1,p2};\n');

fprintf(oFile,'lineloops[nll] = newreg;\n');
fprintf(oFile,'Line Loop(lineloops[nll]) = {c1,c2,c3,c4};\n');
% fprintf(oFile,'Printf("Line Loop number: lineloops[%%g] = %%g.", nll, lineloops[nll]);\n');
fprintf(oFile,'Return\n');

% % --------------------------------------------------------------------------
% %                       generate anode and its interface
% % --------------------------------------------------------------------------
% fprintf(oFile,'\n\n');
% fprintf(oFile,'left = newreg;\n');
% fprintf(oFile,'Line Loop(left) = {1,2,-3,-4};\n');
% fprintf(oFile,'ps1 = newreg;\n');
% fprintf(oFile,'Plane Surface(ps1) = {left};\n');
% fprintf(oFile,'Physical Surface("anode") = {ps1};\n');
% fprintf(oFile,'Physical Line("anInt") = {2};\n');

% --------------------------------------------------------------------------
%                        generate matrix and its interface
% --------------------------------------------------------------------------
% write circles
fprintf(oFile,'\n\n');
fprintf(oFile,'r = %.8e;\n', r);
% fprintf(oFile,'lcar1 = %6.4f;\n', lcar2);
% 
fprintf(oFile,'n= %6d;\n', n_x*n_y);
% 
fprintf(oFile,'For i In {1:n}\n');
fprintf(oFile,'x = coordx[i];\n');
fprintf(oFile,'y = coordy[i];\n');
fprintf(oFile,'nll = i - 1;\n');
fprintf(oFile,'Call CheeseHole;\n');
fprintf(oFile,'intline[4*nll+0] = c1;\n');
fprintf(oFile,'intline[4*nll+1] = c2;\n');
fprintf(oFile,'intline[4*nll+2] = c3;\n');
fprintf(oFile,'intline[4*nll+3] = c4;\n');
fprintf(oFile,'EndFor\n');

fprintf(oFile,'square = newreg;\n');
fprintf(oFile,'Line Loop(square) = {1,2,-3,-4};\n');

fprintf(oFile,'mat = newreg;\n');
fprintf(oFile,'Plane Surface(mat) = {square,lineloops[]};\n');
fprintf(oFile,'Physical Surface("matrix") = {mat};\n');

fprintf(oFile,'Physical Line("int_mat") = {intline[]};\n');

% 
% --------------------------------------------------------------------------
%                       generate cathode and its interface
% --------------------------------------------------------------------------
fprintf(oFile,'\n\n');
% write particles in cathode
fprintf(oFile,'Delete lineloops;\n');
fprintf(oFile,'Delete intline;\n');

fprintf(oFile,'For i In {1:n}\n');
fprintf(oFile,'x = coordx[i];\n');
fprintf(oFile,'y = coordy[i];\n');
fprintf(oFile,'nll = i - 1;\n');
fprintf(oFile,'Call CheeseHole;\n');
fprintf(oFile,'psloops[nll] = newreg;\n');
fprintf(oFile,'Plane Surface(psloops[nll]) = {lineloops[nll]};\n');
fprintf(oFile,'intline[4*nll+0] = c1;\n');
fprintf(oFile,'intline[4*nll+1] = c2;\n');
fprintf(oFile,'intline[4*nll+2] = c3;\n');
fprintf(oFile,'intline[4*nll+3] = c4;\n');
fprintf(oFile,'EndFor\n');
fprintf(oFile,'Physical Surface("inclusion") = {psloops[]};\n');
fprintf(oFile,'Physical Line("int_inc") = {intline[]};\n');
% 

fprintf(oFile,'Physical Line("left")  = {4};\n');
fprintf(oFile,'Physical Line("right") = {2};\n');
fprintf(oFile,'Physical Line("lower") = {1};\n');
fprintf(oFile,'Physical Line("upper") = {3};\n');
fprintf(oFile,'Physical Point("corner") = {1,2,3,4};\n');

fprintf(oFile,'Mesh 2;\n// RefineMesh;\n');
fclose(oFile);
% 
% call Gmsh with cheese.geo as input and generate a 2d mesh (in cheese.msh)
% system("gmsh BatteryPars_test.geo -2 -o BatteryPars_test.msh");
% 
disp(' - Mesh generated');

