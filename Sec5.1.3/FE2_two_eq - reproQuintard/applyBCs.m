% 
% %   fixed temp. at left edge 
% %  
%     leftBC = bound.left;
%     nodesleft = unique(leftBC);
%     for inode = 1 : length(nodesleft)
%       rw = dofArray(nodesleft(inode),1);
%       for cl = 1 : ndofs
%         K(rw,cl) = 0;
%       end
%       K(rw,rw) = 1.0;
%       Res(rw) = -(1 - u(rw,itime));
%     end
%     
% %   fixed temp. at right edge 
% %  
%     rightBC = bound.right;
%     nodesright = unique(rightBC);
%     for inode = 1 : length(nodesright)
%       rw = dofArray(nodesright(inode),1);
%       for cl = 1 : ndofs
%         K(rw,cl) = 0;
%       end
%       K(rw,rw) = 1.0;
%       Res(rw) = -(0 - u(rw,itime));
%     end

    
    
    leftBC = bound.left;
    nodesleft = unique(leftBC);
%     rw = dofArray(nodesleft, 1:2);
    rw = dofArray(nodesleft, 1);
    rw = reshape(rw, [], 1);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(1 - u(rw, itime)); 
    
    rightBC = bound.right;
    nodesright = unique(rightBC);
%     rw = dofArray(nodesright, 1:2);
    rw = dofArray(nodesright, 1);
    rw = reshape(rw, [], 1);
    K(rw, :) = 0;
    K(rw, rw) = eye(length(rw), length(rw));
    Res(rw) = -(0 - u(rw, itime));     
    