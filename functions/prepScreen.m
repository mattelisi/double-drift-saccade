function scr = prepScreen(const)
%
% apparent motion - saccade task v1
%
% Matteo Lisi, 2013
%

HideCursor;

scr.subDist = 100;   % subject distance (cm)
scr.colDept = 32;%16
scr.width   = 435;%410;  % monitor width (mm) (435 desk monitor)

% If there are multiple displays guess that one without the menu bar is the
% best choice.  Dislay 0 has the menu bar.
scr.allScreens = Screen('Screens');
scr.expScreen  = max(scr.allScreens);

% get rid of PsychtoolBox Welcome screen
% Screen('Preference', 'VisualDebugLevel',3);

% set resolution
%if ~const.saveMovie;
%    Screen('Resolution', scr.expScreen, 1024, 768);
%else
%    Screen('Resolution', scr.expScreen, 1344, 1008);
%end

% Open a fullscreen, onscreen window with gray background. Enable 32bpc
% floating point framebuffer via imaging pipeline on it, if this is possible
% on your hardware while alpha-blending is enabled. Otherwise use a 16bpc
% precision framebuffer together with alpha-blending. We need alpha-blending
% here to implement the nice superposition of overlapping gabors. The demo will
% abort if your graphics hardware is not capable of any of this.
%PsychImaging('PrepareConfiguration');
%PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
%[scr.main scr.rect] = PsychImaging('OpenWindow', scr.expScreen, 128);

% Open a window.  Note the new argument to OpenWindow with value 2, specifying the number of buffers to the onscreen window.
%imagingMode = kPsychNeed32BPCFloat;
%[scr.main,scr.rect] = Screen('OpenWindow',scr.expScreen, [0.5 0.5 0.5],[],scr.colDept,2,0,8,imagingMode);
[scr.main,scr.rect] = Screen('OpenWindow',scr.expScreen, [0.5 0.5 0.5],[],scr.colDept,2,0,2);

% get information about  screen
[scr.xres, scr.yres]    = Screen('WindowSize', scr.main);       % heigth and width of screen [pix]

% determine th main window's center
[scr.centerX, scr.centerY] = WindowCenter(scr.main);

% refresh duration
scr.fd = Screen('GetFlipInterval',scr.main);    % frame duration [s]

% Give the display a moment to recover from the change of display mode when
% opening a window. It takes some monitors and LCD scan converters a few seconds to resync.
WaitSecs(2);