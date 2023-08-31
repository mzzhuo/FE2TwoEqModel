
%   fixed temp. at left edge of electrode
%  
    leftBC = bound.left;
    nodesleft = unique(leftBC);
    for inode = 1 : length(nodesleft)
      rw = dofArray(nodesleft(inode),1);
      for cl = 1 : ndofs
        K_mat(rw,cl) = 0;
      end
      K_mat(rw,rw) = 1.0;
      Res_mat(rw) = -(1 - u_mat(rw,itime));
    end

%   fixed temp. at right edge of electrode
%  
    rightBC = bound.right;
    nodesright = unique(rightBC);
    for inode = 1 : length(nodesright)
      rw = dofArray(nodesright(inode),1);
      for cl = 1 : ndofs
        K_mat(rw,cl) = 0;
      end
      K_mat(rw,rw) = 1.0;
      Res_mat(rw) = -(0 - u_mat(rw,itime));
    end
