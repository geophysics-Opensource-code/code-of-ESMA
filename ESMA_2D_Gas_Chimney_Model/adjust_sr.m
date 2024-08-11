function [isx,isz,igx,igz]=adjust_sr(sx,sz,gx,gz,dx,nbc)

% set and adjust the free surface position
isx=round(sx/dx)+1+nbc;isz=round(sz/dx)+1+nbc;
igx=round(gx/dx)+1+nbc;igz=round(gz/dx)+1+nbc;
if abs(sz) <0.5, isz=isz+1; end
igz=igz+(abs(gz)<0.5)*1;
end