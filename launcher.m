%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Code for Experiment 1 (saccade task) taken from
% "Dissociation Between Perceptual and Saccadic localization of moving objects"
% Matteo Lisi & Patrick Cavanagh, 2015, Current Biology
% DOI: 10.1016/j.cub.2015.08.021 
%

clear all;  clear mex;  clear functions;
addpath('functions/');

home;
tic;

% general parameters
const.TEST        = 1;      % 1 = test in dummy mode, 0 = test in eyelink mode
const.gammaLinear = 0;      % use monitor linearization
const.saveMovie   = 1;
const.nTrialMovie = 1;

%
const.gamma    = '../../gammaCalibration/EyelinkBoxcalData.mat';
const.gammaRGB = '../../gammaCalibration/EyelinkBoxcalDataRGB.mat';

% participant ID
newFile = 0;

%
analyzeNow = 0;

while ~newFile
    [vpcode] = getVpCode;

    % create data file
    datFile = sprintf('%s.mat',vpcode);
    
    % dir names
    subDir=substr(vpcode, 0, 4);
    sessionDir=substr(vpcode, 5, 2);
    resdir=sprintf('data/%s/%s',subDir,sessionDir);
    
    if exist(resdir,'file')==7
        o = input('      This directory exists already. Should I continue/overwrite it [y / n]? ','s');
        if strcmp(o,'y')
            newFile = 1;
            % delete files to be overwritten?
            if exist([resdir,'/',datFile])>0;                    delete([resdir,'/',datFile]); end
            if exist([resdir,'/',sprintf('%s.edf',vpcode)])>0;   delete([resdir,'/',sprintf('%s.edf',vpcode)]); end
            if exist([resdir,'/',sprintf('%s',vpcode)])>0;       delete([resdir,'/',sprintf('%s',vpcode)]); end
        end
    else
        newFile = 1;
        mkdir(resdir);
    end
end

currentDir = cd;

% how many consecutive session?
nsess = getTaskInfo;

for sess = 1:str2double(nsess)
    
    cd(currentDir);
    
    % update session number e vpcode, create directory
    actualSess = str2double(sessionDir) + sess -1;
    actualSessStr = num2str(actualSess);
    
    if length(actualSessStr)==1
        actualSessStr = strcat('0',actualSessStr);
    end
    
    if sess > 1
        vpcode = sprintf('%s%s',subDir,actualSessStr);
        
        % create data file
        datFile = sprintf('%s.mat',vpcode);
    
        % dir names
        resdir=sprintf('data/%s/%s',subDir,actualSessStr);
        
        % keep the control to avoid potential deleting of good data
        if exist(resdir,'file')==7
            o = input('      This directory exists already. Should I continue/overwrite it [y / n]? ','s');
            if strcmp(o,'y')
                newFile = 1;
                % delete files to be overwritten?
                if exist([resdir,'/',datFile])>0;                    delete([resdir,'/',datFile]); end
                if exist([resdir,'/',sprintf('%s.edf',vpcode)])>0;   delete([resdir,'/',sprintf('%s.edf',vpcode)]); end
                if exist([resdir,'/',sprintf('%s',vpcode)])>0;       delete([resdir,'/',sprintf('%s',vpcode)]); end
            end
        else
            newFile = 1;
            mkdir(resdir);
        end
        
    end
    
    % prepare screens
    scr = prepScreen(const);
    
    % prepare stimuli
    visual = prepStim(scr, const);
    
    % generate design
    design = genDesign(visual, scr, actualSess);
    
    % prepare movie
    if const.saveMovie
        movieName = sprintf('%s.mov',vpcode);
        % use GSstreamer
        Screen('Preference', 'DefaultVideocaptureEngine', 3)
        const.moviePtr = Screen('CreateMovie', scr.main, movieName, 800, 600, 30, 'CodecSettings= Keyframe=10 Videoquality=1');
        visual.imageRect =  [(scr.centerX-100) (scr.centerY-300) (scr.centerX+700) (scr.centerY+300)];
        
%        const.moviePtr = Screen('CreateMovie', scr.main, movieName, 800, 600, 60, 'CodecSettings= Keyframe=10 Videoquality=1');
%        visual.imageRect =  [scr.centerX+(round(design.ecc*visual.ppd)-400) (scr.centerY-300) (scr.centerX+round(design.ecc*visual.ppd)+400) (scr.centerY+300)];
    end
    
    % initialize eyelink-connection
    [el, err]=initEyelink(vpcode,visual,const,scr);
    if err==el.TERMINATE_KEY
        return
    end
    
    as = mod(actualSess,design.totSession);
    if as==0; as=design.totSession; end
    
    % instructions
    systemFont = 'Arial'; % 'Courier';
    systemFontSize = 19;
    GeneralInstructions = ['Welcome to our experiment. \n\n',...
        'Session ',actualSessStr,' (',num2str(as),' of ',num2str(design.totSession),').\n\n',...
        'At the beginning of each trial, you will see a black dot: look at it to start.\n\n',...
        'A Gabor patch will appear and start moving either on the left or rigth side of the screen. \n\n',...
        'Your task is to make a fast and accurate eye movement to the green flash, as soon as you see it.\n\n'...
        'Press any key when ready to begin.'];
    Screen('TextSize', scr.main, systemFontSize);
    Screen('TextFont', scr.main, systemFont);
    Screen('FillRect', scr.main, visual.bgColor);
    
    DrawFormattedText(scr.main, GeneralInstructions, 'center', 'center', visual.fgColor,70);
    
    % DrawFormattedText(scr.main, 'circular pursuit #1', 20, 20, [0.3 0.3 0.3],70);
    
    Screen('Flip', scr.main);
    
    SitNWait;
    
    try
        % runtrials
        design = runTrials(design,vpcode,el,scr,visual,const);
    catch
        reddUp;
        rethrow(lasterror);
    end
    
    % finalize
    if const.saveMovie
        Screen('FinalizeMovie', const.moviePtr);
    end
    
    % shut down everything
    reddUp;
    
    % save updated design information
    save(sprintf('%s.mat',vpcode),'design','visual','scr','const');
    
    % sposto i risultati nella cartella corrispondente
    movefile(datFile,resdir);
    movefile(vpcode,resdir);
    if ~const.TEST; movefile(sprintf('%s.edf',vpcode),resdir); end
    
    fprintf(1,'\nThis part of the experiment took %.0f min.',(toc)/60);
    fprintf(1,'\n\nOK!\n');
    
    % copy also edf files for analysis
    copyfile(sprintf('%s/%s.edf',resdir,vpcode),'analysis/edf/')
    
    %% ANALYSIS
    
    if analyzeNow
    
        % add vpcode, convert msg2tab and change working directory
        cd('analysis/prepare')
        fid = fopen('subjects.tmp','w');
        fprintf(fid,'%s\n',vpcode);
        fclose(fid);
        
        fid = fopen('subjects.all','a');
        fprintf(fid,'\n%s',vpcode);
        fclose(fid);
        
        cd('../edf')
        fid = fopen('edf2prep','w');
        fprintf(fid,'%s\n',vpcode);
        fclose(fid);
        
        % convert edf2asc
        system('sh prepare2.sh')
    
        % convert edf2asc
        % cd('../edf')
        % system('sh prepare.sh')
        
        % analysis
        cd('../prepare')
        xmsg2tab;
        xanaEyeMovements;
        
        cd('../analysis')
        plot_single_1;
        
    end
    
end



