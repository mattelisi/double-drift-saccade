function design = genDesign(visual,scr,as)
%
% infinite regress #1 saccade task
%
% Matteo Lisi, 2013
% 


%% target parameters

% spatial
design.ecc = [10];                  % eccentricity of the trajectory midpoint (degree of visual angle)
design.alpha = [-30 30];  % tilt of the envelope's trajectory (degree) [NEGATIVE = left tilt]
design.tarFreq = 2;                 % spatial frequency of the carrier (cycles per degree)
design.sigma = 0.1;                 % sigma of the gaussian envelope (sigma = FWHM / 2.35)
design.constrast = 1;               % gabor's contrast [0 - 1]
design.lambda = (1/design.tarFreq); % wavelength of the carrier
design.flashSigma = 0.18;

% motion & temporal
design.envelope_speed = 2; %4;      % degree/sec.
design.drifting_speed = [1.5 1.5];%[0 3];  % (cycles of the carrier)/sec.
design.trajectoryLength = 3;    % degree of visual angle
design.motionPeriod = ((2*design.trajectoryLength)/design.envelope_speed);

design.fixDur = 0.400;          % minimum fixation duration [s]
design.fixDuJ = 0.200;          % additional fixation duration jitter [s]

design.motionType = 'triangular';    % allowed: triangular, sinusoidal

design.tiltedMotion = 1;        % put 0 here if you want vertical motion of the envelope
                                % in that case the sign of alpha gives the
                                % relation between upward-downward motion
                                % of the envelope and left-right drift

%% other parameter

% task design
design.removeTarget = 1;         % if you want to remove target after saccade
design.flashTime = [0.5:0.1:1];  % flash time (multiplier of motion period)
                                 % NB: qui corrisponde 



% task structure
design.nTrialsInBlock = 32;
design.nTrlsBreak = 200;    % number of trials between breaks, within a block
design.iti = 0.2;
design.totSession = 1;
rep = 16;

% timing
design.preRelease = scr.fd/3;           % must be half or less of the monitor refresh interval
design.fixDur = 0.400;                  % minimum fixation duration [s]
design.fixDuJ = 0.200;                  % additional fixation duration jitter [s]
design.maxSacRt = 10;

%% trials list
t = 0;
for r = 1:rep
for alpha = design.alpha
for tarFreq = design.tarFreq
for contrast = design.constrast
for sigma = design.sigma
for env_speed = design.envelope_speed
for drift_speed = design.drifting_speed
for trajLength = design.trajectoryLength
for initPos = [-1 1]
for ecc = design.ecc
for flashTime = design.flashTime
    
    t = t+1;
    
    % settings
    trial(t).alpha = alpha;
    
    trial(t).tarFreq = tarFreq / visual.ppd;    % in trial list the measures are in pixels
    trial(t).sigma = sigma * visual.ppd;        % 
    
    trial(t).contrast = contrast;
    trial(t).env_speed = env_speed;
    trial(t).drift_speed = drift_speed;
    trial(t).trajLength = trajLength;
    trial(t).initPos = initPos;
    trial(t).ecc = ecc;
    trial(t).flashTime = flashTime;
    trial(t).flashSigma = design.flashSigma * visual.ppd; 

    
    % compute motion parameters for a given trial
    
    % envelope motion
    trial(t).period = (trajLength*2) / env_speed;       % sec
    timeIndex = linspace(0, 1, round(trial(t).period/scr.fd)) + 0.25;
    switch design.motionType
        case 'triangular'
            tarRad = initPos * (trajLength/2) * sawtooth(2*pi*timeIndex, 0.5);  % target radius for each frame
        case 'sinusoidal'
            tarRad = initPos * (trajLength/2) * sin(2*pi*(timeIndex+0.25));
    end
    tarRad(end) = []; 
    if design.tiltedMotion
        tarAng = repmat(deg2rad(90 + alpha), 1, length(tarRad));
    else
        tarAng = repmat(deg2rad(90), 1, length(tarRad));
    end
    [x, y] = pol2cart(tarAng, tarRad);  % x-y coord. of the envelope for each frame
    trial(t).tarPos = [x' y'];
    
    % number of frames for a single cycle
    tarfra = length(tarRad);
    
    % drifting motion
    phaseRange = (drift_speed * 360) * (trial(t).period)/2;
    
    % grating phase at each frame
    initialRandomPhase = rand*180;
    switch design.motionType
        case 'triangular'
            tarPhase = initialRandomPhase + initPos*(-sign(alpha)) * phaseRange * sawtooth(2*pi*timeIndex, 0.5);
        case 'sinusoidal'
            tarPhase = initialRandomPhase + initPos*(-sign(alpha)) * phaseRange * -cos(2*pi*(timeIndex+0.25));
    end
                 
    tarPhase(end) = [];
    trial(t).tarPhase = tarPhase;

    %
    trial(t).fixDur = round((design.fixDur + design.fixDuJ*rand)/scr.fd)*scr.fd;
    trial(t).fixLoc = [scr.centerX scr.centerY];
    
%     % events
%     trial(t).event = zeros(design.nMaxMotionCycles*tarfra,1);
%     trial(t).event(ceil(tarfra*flashTime)) = 1;
    
end
end
end
end
end
end
end
end
end
end
end

design.totTrials = t;

% select trial for session
as = mod(as,design.totSession); 
if as==0; as=design.totSession; end

design.actualSession = as;

sessIndex = (repmat(1:design.totSession,1,ceil(design.totTrials/design.totSession)));
trial = trial(sessIndex==as);

% random order
r = randperm(length(trial));
trial = trial(r);

% generate blocks
design.nBlocks = length(trial)/design.nTrialsInBlock;
design.blockOrder = 1:design.nBlocks;

b=1; beginB=b; endB=design.nTrialsInBlock;

for i = 1:design.nBlocks
    design.b(i).trial = trial(beginB:endB);
    beginB  = beginB + design.nTrialsInBlock;
    endB    = endB   + design.nTrialsInBlock;
end



