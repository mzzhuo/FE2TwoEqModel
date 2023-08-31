% Copyright (C) 2007 Garth N. Wells
%
% Write VTK post-processing files
%
% Modified by VP Nguyen for the IGA FEM code

function toParaview(coords,connect,etype,vtuFile,u)
% 
connect = [connect.mat; connect.inc];
dim   = size(coords,2);
nnode = size(coords,1);
ncell = size(connect,1);

% Output files

outfileVTU  = strcat(vtuFile, '.vtu');
results_vtu = fopen(outfileVTU, 'wt');

if(strcmp(etype, 'Quad4') || strcmp(etype, 'Quad8') || strcmp(etype, 'Quad9'))
    numVertexesPerCell = 4;
    VTKCellCode = 9;
elseif(strcmp(etype, 'B8'))
    numVertexesPerCell = 8;
    VTKCellCode = 12;
elseif(strcmp(etype, 'Tri3') || strcmp(etype, 'Tri6'))
    numVertexesPerCell = 3;
    VTKCellCode = 5;
elseif(strcmp(etype, 'Tet4'))
    numVertexesPerCell = 4;
    VTKCellCode = 10;
else
    error('Element type not known (VTKPostProcess)')
end


dof_per_vertex = 2;


%% Write headers
fprintf(results_vtu, '<VTKFile type="UnstructuredGrid"  version="0.1"   > \n');
fprintf(results_vtu, '<UnstructuredGrid> \n');
fprintf(results_vtu, '<Piece  NumberOfPoints="  %g" NumberOfCells=" %g"> \n', nnode, ncell);

%% Write point data
fprintf(results_vtu, '<Points> \n');
if( dof_per_vertex == 1)
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="1"  format="ascii" > \n');
else
    fprintf(results_vtu, '<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii" > \n');
end

for i=1:nnode
    if( dim == 3)
        fprintf(results_vtu, '%f ',  coords(i,1:3));
    elseif(dim == 2)
        fprintf(results_vtu, '%f ',  coords(i,1:2));
        fprintf(results_vtu, '0.0 ');
    end
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Points> \n');

%% Print cells
fprintf(results_vtu, '<Cells> \n');

%% Print cell connectivity
fprintf(results_vtu, '<DataArray  type="Int32"  Name="connectivity"  format="ascii"> \n');

for i=1:ncell
    fprintf(results_vtu, '%g ',  connect(i,1:numVertexesPerCell)-1 );
    fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell offsets
fprintf(results_vtu, '<DataArray  type="Int32"  Name="offsets"  format="ascii"> \n');

offset = 0;
for i=1:ncell
    offset = offset + numVertexesPerCell;
    fprintf(results_vtu,'%.20g ', offset);
%     fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');

% Print cell types
fprintf(results_vtu, '<DataArray  type="UInt8"  Name="types"  format="ascii"> \n');

for i=1:ncell
    fprintf(results_vtu, '%.20g ', VTKCellCode);
%     fprintf(results_vtu, '\n');
end

fprintf(results_vtu, '</DataArray> \n');
fprintf(results_vtu, '</Cells> \n');

%% Print result data
% print Concentration field
    
fprintf(results_vtu, '<PointData  Scalars="DataArray"> \n');


fprintf(results_vtu, '<DataArray  type="Float64"  Name="u" NumberOfComponents="1" format="ascii"> \n');
for i=1:nnode
    fprintf(results_vtu, '%12.8f   \n', u(i) );
end
fprintf(results_vtu, '</DataArray> \n');

fprintf(results_vtu, '</PointData> \n');

% end of VTK file

fprintf(results_vtu, '</Piece> \n');
fprintf(results_vtu, '</UnstructuredGrid> \n');
fprintf(results_vtu, '</VTKFile> \n');
% 
fclose(results_vtu);
end
