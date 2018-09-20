% xanaEyeMovements.m
%
% Input:
% 1. vpcode.tabc
% 2. vpcode.dat
%
% Output:
% 2. vpcode.rea - saccadic responses
% 3. vpcode.btr - exclusion criteria 
%


clear;
home;

addpath('functions/');
addpath('../../Functions/');

% switches
welche    = 1;          % which files shall be analysed (1=all, 2=tmp)

plotData  = 0;          % show figure with raw data
printImages = 0;

SAMPRATE  = 1000;       % Eyetracker sampling rate 
maxRT     = 600;        % time window for corrective saccades
timBefCue = 700;        % time window to be analysed before saccade cue
velSD     = 5;          % lambda for microsaccade detectionc
minDur    = 8;          % threshold duration for microsaccades (ms)
VELTYPE   = 2;          % velocity type for saccade detection
maxMSAmp  = 1;          % maximum microsaccade amplitude
crit_cols = [2 3];      % critical gaps in dat files to find missings data
mergeInt  = 10;         % merge interval for subsequent saccadic events

% target radius for fixation and saccadic eye movements
fixRad = 2.5;
tarRad = 3.5;

% parameters for pursuit analysis
pursuitInterval = [10 80];     % start-end of the to be analysed interval after primary saccade (ms)

addFilter = 0;                  % logical: should additional filtering (below) be applied on the post-saccadic trace?
low_pass_pos=300;               % Position filter (Hz)
low_pass_vel=80;                % Velocity filter (Hz)
low_pass_acc=60;                % Acceleration filter (Hz)
order=16;                       % Order Nth for the filter
minPursuitDur = 30;             % how many ms of 'clean' (i.e., after excluding corrective saccades) pursuit are required for the analysis?


% Paris
MO_WIDE   = 1024;           % x resolution
ABSTAND   = 60;             % subject distance in cm
MO_PHYS   = 41.0;           % monitor width in cm
scrCen    = [1024, 768]/2;  % screen center (intial fixation position)

DPP = pix2deg(MO_WIDE,MO_PHYS,ABSTAND,1); % degrees per pixel
PPD = deg2pix(MO_WIDE,MO_PHYS,ABSTAND,1); % pixels per degree

% do analyses
anaEyeMovementsFilter_reloaded;

% 
% copyfile(reaAllFile,reaAllFile_R);

xcombineData;

fprintf(1,'\n\nOK!!\n');

% to save single traces
% exportfig(1,'ML02_1.9.eps','Color','rgb')
% exportfig(h3,'SB01_PS.eps','Color','rgb')

