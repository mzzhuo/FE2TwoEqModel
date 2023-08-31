
%   fixed temp. at left edge of electrode
%  
    leftBC = bound.left;
    nodesleft = unique(leftBC);
    for inode = 1 : length(nodesleft)
      rw = dofArray(nodesleft(inode),1);
      for cl = 1 : ndofs
        K(rw,cl) = 0;
      end
      K(rw,rw) = 1.0;
      Res(rw) = -(1 - u(rw,itime));
    end

%   fixed temp. at right edge of electrode
%  
    rightBC = bound.right;
    nodesright = unique(rightBC);
    for inode = 1 : length(nodesright)
      rw = dofArray(nodesright(inode),1);
      for cl = 1 : ndofs
        K(rw,cl) = 0;
      end
      K(rw,rw) = 1.0;
      Res(rw) = -(0 - u(rw,itime));
    end
