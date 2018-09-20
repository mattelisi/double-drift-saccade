function [theta] = compDirection(x,y,type)
%
% Determine direction angle (radians) for the ordered position vector (x,y)
% with a linear fit to the data. 
%
% Matteo Lisi, 2013
%

if nargin<3
    type = 'average';
end

% determine hemi-field
q = sign(sum(diff(x)));

switch lower(type)
    
    case 'linear'
        
        coeff = polyfit(x, y, 1); % this include a non-zero intercept
        switch q
            case -1
                theta = pi + coeff(1);
            case 1
                theta = coeff(1);
        end
        
        
    case 'zerointercept'
        
        switch q
            case -1
                theta = pi + x\y;
            case 1
                theta = coeff(1);
        end
        
        
    case 'average'
        
        index = sort(repmat(1:3,1,floor(length(x)/3)));
        index = [index, repmat(3,1,mod(length(x),3))];
        x_start = mean(x(index==1));
        y_start = mean(y(index==1));
        x_end = mean(x(index==3));
        y_end = mean(y(index==3));
        theta = cart2pol(x_end-x_start, y_end-y_start);
        
end



% x = linspace(0,0.025,100) + randn(1,100)./10;
% y = linspace(0,7,100) + randn(1,100)./10;
% x=x'
% y=y'
% plot(x,y)


% a = [1 2 3 4 8 9 16 17 18];
% 
% b=a-(1:length(a));
% b = [true; diff(b(:)) ~= 0; true];
% split = mat2cell(a(:).', 1, diff(find(b)));
% 
% split{:}