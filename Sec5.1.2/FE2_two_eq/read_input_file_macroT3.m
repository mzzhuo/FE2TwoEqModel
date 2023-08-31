%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==============================Function read_input_file =================
function [nnode,coords,connect,bound] = read_input_file_speT3(infile) 
% 
%   clear;clc;
%   addpath('./utilities','./inputfiles','results');
%   infile = 'SPET3';
%   
% read gmsh input file
% 
  infile = [infile, '.msh'];
  infile = fopen(infile,'r');
% 
%%
  line = fgetl(infile);
%   
% get physical line number for interface elems
  while( ~strcmp(line,'$PhysicalNames') )
      line = fgetl(infile);
  end
  line = fgetl(infile);
  nphy = sscanf(line, '%d'); % number of physical names
  %
  phyname = cell(nphy,2);
  %
  for i = 1 : nphy
      line = fgetl(infile);
      buf = textscan(line,'%d %d %s');
      phyname{i,1} = buf{2}; phyname{i,2} = buf{3};
  end
%%
% number of nodes 
  while( ~strcmp(line,'$Nodes') )
      line = fgetl(infile);
  end
  %
  line = fgetl(infile);
  nnode = sscanf(line, '%d');
 %
 % node coordinates stored in coords
  coords = zeros(nnode,3); % gmsh 2D/3D mesh both have 3 components
  for i = 1 : nnode
      line = fgetl(infile);
      buf = sscanf(line,'%f %f %f %f');
      coords(i,:) = buf(2:4)'; % not store node ID! make sure is 1,2,3....
  end
  coords = coords(:,1:2); % this is 2D problem
%%
%  No. elements && connectivity
  while( ~strcmp(line,'$Elements') )
      line = fgetl(infile);
  end
  line = fgetl(infile);
  nelem = sscanf(line, '%d');
%%
% gmsh file format:
% elem ID, elem type, no of tags, tags...,
% 1 2-node line.
% 2 3-node triangle.
% 3 4-node quadrangle.
%
  connect = zeros(nelem,8);
  for i = 1 : nelem
      line = fgetl(infile);
      buf = sscanf(line,'%f %f %f %f %f %f %f %f');
      connect(i,1:length(buf)) = buf; % elem ID not stored!
  end
% 
  connect(:,1:3) = [];
%
%% correct connectivity of elements with negative jacobian
%%
% Get nodes of domains: SPEC and electrodes
% 
  if ( string(phyname{3,2}) == '"spe"' )
      connect_spe = connect((connect(:,1) == phyname{3,1}), 3:5); % SPE 
  else
      fprintf('SPE physical number not found!');
  end
%   
%%
%
% node_BC0 = find( coords(:,1) == 0 );
% node_BC1 = find( coords(:,1) == max(coords(:,1)) );
%
  if ( string(phyname{1,2}) == '"left"' )
      bound.left = connect((connect(:,1) == phyname{1,1}),3:4);
  else
      fprintf('SPE left boundary physical number not found!');
  end 
  if ( string(phyname{2,2}) == '"right"' )
      bound.right = connect((connect(:,1) == phyname{2,1}),3:4);
  else
      fprintf('SPE right boundary physical number not found!');
  end 
%%
% combine individual domains
  clear connect;
  connect = connect_spe;
%   
%
  fclose(infile);
end