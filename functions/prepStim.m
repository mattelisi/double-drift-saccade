function[visual] = prepStim(scr, const)
%
% apparent motion demo
%
% Prepare display parameters & similar
%

visual.ppd = va2pix(1,scr);   % pixel per degree

visual.black = BlackIndex(scr.main);
visual.white = WhiteIndex(scr.main);

visual.bgColor = (visual.black + visual.white) / 2;     % background color
visual.fgColor = visual.black;

visual.scrCenter = [scr.centerX scr.centerY scr.centerX scr.centerY];

visual.fixCkRad = 2.5*visual.ppd;      % fixation check radius
visual.fixCkCol = visual.black;      % fixation check color
visual.fixCol = 50;

% target
visual.tarSize = 200;
visual.res = 1*[visual.tarSize visual.tarSize];

% gamma correction
if const.gammaLinear
    load(const.gamma);
    load(const.gammaRGB);
    
    % prepare and load lookup gamma table
    luminanceRamp = linspace(LR.LMin, LR.LMax, 256);
    invertedRamp = LR.LtoVfun(LR, luminanceRamp);
    invertedRamp = invertedRamp./255;
    % plot(invertedRamp)
    
    inverseCLUT = repmat(invertedRamp',1,3);
    % save gammaTable_greyscale.mat inverseCLUT
    
    Screen('LoadNormalizedGammaTable', scr.main, inverseCLUT);
    
    visual.bgColor = 14;
    visual.bgColorLuminance = LR.VtoLfun(LR, invertedRamp(visual.bgColor)*255);
end

% set priority of window activities to maximum
priorityLevel=MaxPriority(scr.main);
Priority(priorityLevel);


 