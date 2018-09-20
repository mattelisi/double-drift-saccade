function[texture] = makeColorPatch(win,visual,color,sigma,width)
%
% make a simple white gaussian patch (texture, fast presentation)
%

% if is clipped on the sides, increase width.
halfWidthOfGrid = width / 2;
widthArray = (-halfWidthOfGrid) : halfWidthOfGrid;  

% color codes
highlight = color;
background = visual.bgColor; 

if length(background)==1
    background=repmat(background,1,3);
end
if length(highlight)==1
    highlight=repmat(highlight,1,3);
end

difference = highlight - background;
[x y] = meshgrid(widthArray, widthArray);
imageMatrix = exp(-((x .^ 2) + (y .^ 2)) / (sigma ^ 2));
if length(difference)==1
    grayscaleImageMatrix = background + difference * imageMatrix;
else
    grayscaleImageMatrix(:,:,1) = background(1) + difference(1) * imageMatrix; 
    grayscaleImageMatrix(:,:,2) = background(2) + difference(2) * imageMatrix; 
    grayscaleImageMatrix(:,:,3) = background(3) + difference(3) * imageMatrix; 
end
texture = Screen('MakeTexture', win, grayscaleImageMatrix); 