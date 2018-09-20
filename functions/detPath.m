function targetPos = detPath(alpha,dir,step,design,visual)
%
% tangential apparent motion demo
%
% calculate apparent motion path:
%
% alpha = angular position of middle target
% dir = movement direction: 1= clockwise; -1 = counterclockwise
%
% return target positions matrix, dim [n,2]
%
% Matteo Lisi, 2013
%

tar = [design.displaySize*cos(alpha), design.displaySize*sin(alpha)];
gamma = pi/2-alpha;
targetPos = zeros(step,2);

for k=1:step
    v=dir*(k+floor(step/2)-step);
    targetPos(k,:)=visual.scrCenter(1:2)+[tar(1)+v*design.d*cos(gamma), tar(2)-v*design.d*sin(gamma)];
end