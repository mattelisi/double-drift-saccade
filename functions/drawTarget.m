function drawTarget(col,loc,size,scr)
%
% tangential apparent motion demo
%
% draw circular target
%

pu = size/2;
Screen(scr.main,'FillOval',col,repmat(loc,1,2)+[-pu -pu pu pu]);