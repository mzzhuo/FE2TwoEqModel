%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% By: Mingzhao ZHUO (2018, Delft), m.zhuo@tudelft.nl
% (based on Lagrange multipler code by Davide Grazioli and Mohsen Goudarzi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==============================Function read_input_file =================
function [nnode,coords,connect,interface,bound] = read_RVE_inc_mesh(infile,elemType)
% 
%   clear;clc;
%   addpath('./utilities','./inputfiles','results');
%   infile = 'RVE';
%   
% read gmsh input file
% 
  infile = [infile, '.msh'];
  infile = fopen(infile,'r');
  if strcmp(elemType,'Tri3')
    nelnodes = 3;
  elseif strcmp(elemType,'Quad4') 
    nelnodes = 4;   % triangle element
  end  
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
  connect = zeros(nelem,5+nelnodes);
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

  if ( string(phyname{8,2}) == '"matrix"' )
      connect_mat = connect((connect(:,1) == phyname{8,1}), 3:2+nelnodes); % SPE 
  else
      fprintf('SPE physical number not found!');
  end
  if ( string(phyname{9,2}) == '"inclusion"' )
      connect_inc = connect((connect(:,1) == phyname{9,1}), 3:2+nelnodes); % SPE 
  else
      fprintf('SPE physical number not found!');
  end

%   
%%
%
% node_BC0 = find( coords(:,1) == 0 );
% node_BC1 = find( coords(:,1) == max(coords(:,1)) );
%
  if ( string(phyname{4,2}) == '"left"' )
      bound.left = connect((connect(:,1) == phyname{4,1}),3:4);
  else
      fprintf('SPE left boundary physical number not found!');
  end 
  if ( string(phyname{5,2}) == '"right"' )
      bound.right = connect((connect(:,1) == phyname{5,1}),3:4);
  else
      fprintf('SPE right boundary physical number not found!');
  end 
  if ( string(phyname{6,2}) == '"lower"' )
      bound.lower = connect((connect(:,1) == phyname{6,1}),3:4);
  else
      fprintf('SPE right boundary physical number not found!');
  end  
  if ( string(phyname{7,2}) == '"upper"' )
      bound.upper = connect((connect(:,1) == phyname{7,1}),3:4);
  else
      fprintf('SPE right boundary physical number not found!');
  end 
% Get interface nodes
% 
% left electrode interface with right spe 
  if ( string(phyname{2,2}) == '"int_mat"' )
      int_mat = connect((connect(:,1) == phyname{2,1}),3:4);
  else
      fprintf('Matrix interface physical number not found!');
  end 
  if ( string(phyname{3,2}) == '"int_inc"' )
      int_inc = connect((connect(:,1) == phyname{3,1}),3:4);
  else
      fprintf('Inclusion interface physical number not found!');
  end 
%  corner nodes
  if ( string(phyname{1,2}) == '"corner"' )
      bound.corner = connect((connect(:,1) == phyname{1,1}),3);
  else
      fprintf('SPE corner physical number not found!');
  end 
%%
% create matched master and slave node list
% master: left + lower; slave: right + upper
master1 = bound.left(2:end,1);
slave1  = bound.right(2:end,1);
if ( sum(abs(coords(master1,2)-coords(slave1,2))) ~= 0 )
    error('Wrong matched node ordering!')
end
master2 = bound.lower(2:end,1);
slave2  = bound.upper(2:end,1);
if ( sum(abs(coords(master2,1)-coords(slave2,1))) ~= 0 )
    error('Wrong matched node ordering!')
end    
bound.master = [master1;master2];    
bound.slave  = [slave1;slave2];
%%
% create interface nodes
int_mat_nodes = int_mat(:,1);
int_inc_nodes = int_inc(:,1);
%   
addLambda = [ int_mat_nodes int_inc_nodes ];
% check node pair at same coords
for ind = 1 : size(addLambda,1)
    nodecp = addLambda(ind,:);
    if (norm(coords(nodecp(1),:)-coords(nodecp(2),:)) ~= 0)
        error('Wrong interface node ordering!')
    end
end

nnodepair  = size(addLambda,1);    % node no of both interface
% define interface nodes and elements
int_node = (1:nnodepair)';       % new node number starting from 1
interface.nodes   = [ int_node addLambda ];

% combine individual domains
  clear connect;
  connect.mat = connect_mat; 
  connect.inc = connect_inc;

%   clear connect;
%   connect = [connect_mat; connect_inc];
% % 0 - matrix; 1 - inclusion
%   elemTag = [zeros(size(connect_mat,1),1); ones(size(connect_inc,1),1)];
%
  fclose(infile);
end