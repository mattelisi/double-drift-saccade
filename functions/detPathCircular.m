function targetPos = detPathCircular(alpha, dir, stepBefore, stepAfter, beta, ecc, visual, rev)
%
% circular apparent motion - saccade task v1.0
%
% calculate apparent motion path:
%
% alpha = angular position of target (flash)
% dir = movement direction: 1= clockwise; -1 = counterclockwise
%
% return target positions matrix, dim [n,2]
%
% Matteo Lisi, 2013
%
% update: - added optional parameter 'rev' (flash at reversal) 
%         - require the angle and not the step in input (1/5/2013)
%

if nargin <= 7
    rev = 0;
end

ecc = ecc*visual.ppd;           % convert radius in pixel
targetPos = zeros((stepBefore+stepAfter+1),2);

if ~rev
    angles = 2*pi-linspace((alpha+dir*stepBefore*beta), (alpha-dir*stepAfter*beta), stepBefore+stepAfter+1);
    for k=1:(stepBefore+stepAfter+1)
        targetPos(k,:)=visual.scrCenter(1:2)+[ecc*cos(angles(k)), ecc*sin(angles(k))];
    end
else
    % direction reversal at flash location (alpha)
    angles_1 = 2*pi-linspace((alpha+dir*stepBefore*beta), (alpha+dir*beta), stepBefore);
    angles_2 = 2*pi-linspace(alpha, (alpha+dir*stepAfter*beta), stepAfter+1);
    angles = [angles_1, angles_2];
    
    for k=1:(stepBefore+stepAfter+1)
        targetPos(k,:)=visual.scrCenter(1:2)+[ecc*cos(angles(k)), ecc*sin(angles(k))];
    end
end
