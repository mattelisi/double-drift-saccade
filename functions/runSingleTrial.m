function data = runSingleTrial(td, scr, visual, const, design)
%
% saccade task - drifting gabor with envelope moving orthogonally
%
% td = trial design
%
% Matteo Lisi, 2013
% 


%% TRIAL PREP.

% clear keyboard buffer
FlushEvents('KeyDown');

if const.TEST == 1
    ShowCursor;
else
    HideCursor;
end

% predefine boundary information
cxm = round(td.fixLoc(1));
cym = round(td.fixLoc(2));
chk = round(visual.fixCkRad);

% draw trial information on EyeLink operator screen
Eyelink('command','draw_filled_box %d %d %d %d 15', round(cxm-chk/8), round(cym-chk/8), round(cxm+chk/8), round(cym+chk/8));    % fixation
Eyelink('command','draw_box %d %d %d %d 15', cxm-chk, cym-chk, cxm+chk, cym+chk);                   % fix check boundary
Eyelink('command','draw_cross %d %d', cxm, cym);

% make target and compute path
gabor = CreateProceduralGabor(scr.main, visual.tarSize, visual.tarSize, 0, [visual.bgColor/255 visual.bgColor/255 visual.bgColor/255 0.0], 1, 0.5);
flash = makeColorPatch(scr.main, visual, [0 255 0], td.flashSigma, visual.tarSize); % makeGaussianPatch(scr.main,visual,td.sigma, visual.tarSize)
pathRects = detRect(visual.ppd*td.tarPos(:,1) + scr.centerX + visual.ppd*td.ecc, visual.ppd*td.tarPos(:,2) + scr.centerY, visual.tarSize);
if design.tiltedMotion
    gaborAngle = - td.alpha;
else
    gaborAngle = 0;
end

% predefine time stamps
tBeg    = NaN;
tSac    = NaN;
tFlash  = NaN;
tEnd    = NaN;

eyePhase = 1;   % 1 = fixation phase (before cue), 2 = saccade phase (after cue)

eventStr = [];

% flags/counters
ex_fg = 0;      % 0 = ongoing; 1 = saccade OK; 2 = fix break; 3 = too slow
cycle = 0;
eventStr = 'EVENT_TargetOnset';

% draw fixation stimulus
drawFixation(visual.fixCol,td.fixLoc,scr,visual);
tFix = Screen('Flip', scr.main,0);
Eyelink('message', 'EVENT_FixationDot');
if const.TEST; fprintf(1,strcat('\n','EVENT_FixationDot')); end
if const.saveMovie
    eyePosmovie = [];
    Screen('AddFrameToMovie', scr.main, visual.imageRect, 'frontBuffer', const.moviePtr, round(td.fixDur/scr.fd)); 
end

tFlip = tFix + td.fixDur;
WaitSecs(td.fixDur - 2*design.preRelease);  

%% TRIAL

while ~ex_fg
    
    for i = 1:length(td.tarPos)
        
        % draw stimuli
        Screen('DrawTexture', scr.main, gabor, [], pathRects(i,:), gaborAngle, [], [], [], [], kPsychDontDoRotation, [td.tarPhase(i), td.tarFreq, td.sigma, td.contrast, 1, 0, 0, 0]);
        
        if ((cycle + i/length(td.tarPos)) <= td.flashTime)
            % draw fixation
            drawFixation(visual.fixCol,td.fixLoc,scr,visual);
            
        elseif eyePhase == 1
            eyePhase = 2;
            eventStr = 'EVENT_Flash';
        end
            
        
        
        % drawing finished, flip canvas
        Screen('DrawingFinished',scr.main);
        tFlip = Screen('Flip', scr.main, tFlip + scr.fd - design.preRelease);
        
        % send event triggers to eyelink
        if ~isempty(eventStr)
            Eyelink('message', eventStr);
            if const.TEST; fprintf(1,strcat('\n',eventStr)); end
            switch eventStr
                case 'EVENT_TargetOnset'
                    tBeg = tFlip;
                case 'EVENT_Flash'
                    tFlash = tFlip;
            end     
            eventStr = [];
        end
        
        % get eye position data
        [x,y] = getCoord(scr, const);   
        
        switch eyePhase
            case 1      % fixation phase
                if sqrt((mean(x)-cxm)^2+(mean(y)-cym)^2)>chk    % check fixation in a circular area
                    
                    ex_fg = 2;     % fixation break
                    
                    % blank screen after fixation break
                    drawFixation(visual.fixCol,td.fixLoc,scr,visual);
                    Screen(scr.main,'Flip');
                    break
                    
                end
                
            case 2      % saccade phase
                if sqrt((mean(x)-cxm)^2+(mean(y)-cym)^2)>chk 
                    
                    tSac = GetSecs;
                    
                    Eyelink('message', 'EVENT_Saccade1Started');
                    if const.TEST; fprintf(1,'\nEVENT_Saccade1Started'); end
                    ex_fg = 1;    % successful trial (with saccade)
                    
                    % blank screen after saccade
                    Screen(scr.main,'Flip');
                    break
                    
                elseif (GetSecs - tFlash) > design.maxSacRt
                    
                    ex_fg = 3;    % too slow
                    
                    % blank screen after fixation break
                    drawFixation(visual.fixCol,td.fixLoc,scr,visual);
                    Screen(scr.main,'Flip');
                    break
                    
                end
        end 

        if const.saveMovie; Screen('AddFrameToMovie', scr.main, visual.imageRect, 'frontBuffer', const.moviePtr, 1); end
        
    end % end of a stimulus motion cycle
    
    cycle = cycle + 1;
    
end

%% trial end

switch ex_fg
    
    case 2
        data = 'fixBreak';
        Eyelink('command','draw_text 100 100 15 Fixation break');
        
    case 3
        data = 'tooSlow';
        Eyelink('command','draw_text 100 100 15 Too slow');
        
    case 1
        
        WaitSecs(0.2);
        if const.saveMovie; Screen('AddFrameToMovie', scr.main, visual.imageRect, 'frontBuffer', const.moviePtr, round(0.2/scr.fd)); end

        % collect trial information
        trialData = sprintf('%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f',[td.alpha td.env_speed td.drift_speed td.trajLength td.initPos td.ecc td.flashTime]);
        
        % determine presentation times relative to 1st frame of motion
        timeData = sprintf('%i\t%i\t%i\t%i\t%i',round(1000*([tFix tBeg tFlash tSac tEnd]-tBeg)));
        
        % determine response data
        respData = sprintf('%i',round(1000*(tSac - tFlash)));
        
        % collect data for tab [6 x trialData, 5 x timeData, 1 x respData]
        data = sprintf('%s\t%s\t%s',trialData, timeData, respData);
end


% close active textures
Screen('Close', gabor)

