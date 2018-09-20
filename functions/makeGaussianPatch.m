function[texture] = makeGaussianPatch(win,visual,sigma,width)
%
% make a simple white gaussian patch (texture, allows fast presentation)
%

% if is clipped on the sides, increase width.
halfWidthOfGrid = width / 2;
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  

% color codes
white = visual.white;
gray = visual.bgColor; 
  
absoluteDifference = abs(white - gray);
[x y] = meshgrid(widthArray, widthArray);
imageMatrix = exp(-((x .^ 2) + (y .^ 2)) / (sigma ^ 2));
grayscaleImageMatrix = gray + absoluteDifference * imageMatrix;
texture = Screen('MakeTexture', win, grayscaleImageMatrix); 