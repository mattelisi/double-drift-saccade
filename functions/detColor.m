function colors = detColor(flash,step,visual)
% 
% determine color sequence (flash in the middle)
%
%
colors = repmat(visual.tarCol,step,1);
if flash
    colors(ceil(step/2),:) = visual.flashCol;
end