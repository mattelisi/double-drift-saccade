%% Drifting Gabor - moving envelope - eye movement analysis
%
% Matteo Lisi, 2013
%

if welche == 1
    reaAllFile = sprintf('../rea/AllSubjects_V%iD%iT%iS%i.rea',velSD,minDur,VELTYPE,SAMPRATE);
    reaAllFile_R = sprintf('../R/AllSubjects_V%iD%iT%iS%i.txt',velSD,minDur,VELTYPE,SAMPRATE);
    btrAllFile = sprintf('../btr/AllSubjects_V%iD%iT%iS%i.btr',velSD,minDur,VELTYPE,SAMPRATE);
    
    freaAll = fopen(reaAllFile,'w');
    fbtrAll = fopen(btrAllFile,'w');
    fRAll = fopen(reaAllFile_R,'w');
    
    fprintf(fRAll,  'code\tsub\tsess\tblock\ttrial\talpha\tenvSpeed\tdriftSpeed\tpathLength\tinitPos\tecc\tfixOff\ttFix\ttBeg\tFixoff\ttSac\tEnd\tsacRT\tother1\tother2\tother3\ttEdfFix\ttEdfTOn\ttEdfFixOff\tEdfSac\t');
    fprintf(fRAll,  'sacRT\treaSacNumber\tsacType\tsacOnset\tsacOffset\tsacDur\tsacVPeak\tsacDist\tsacAngle1\tsacAmp\tsacAngle2\tsacxOnset\tsacyOnset\tsacxOffset\tsacyOffset\t');
    fprintf(fRAll,  'driftDirLand\tenvDirAng\tdriftDirAng\tpursuitDirAng\tpos_p_X\tpos_p_Y\tvel_p_X\tvel_p_Y\tacc_p_X\tacc_p_Y\ttarPosLandX\ttarPosLandY\ttarPosSacOnX\ttarPosSacOnY\n');
    
    
    sfid = fopen('subjects.all','r');
else
    sfid = fopen('subjects.tmp','r');
end

%% prepare figures and axes (if plotting is required)
if plotData
    
    cout = [0.3 0.3 0.3];
    cbac = [1.0 1.0 1.0];
    
    close all;
    h1 = figure;
    set(gcf,'pos',[50 50 1000 500],'color',cbac);
    ax(1) = axes('pos',[0 0 0.5 1]);    % left bottom width height
    ax(2) = axes('pos',[0.55 0.5 0.5 0.5]);
    ax(3) = axes('pos',[0.55 0.0 0.5 0.5]);
     
    
    % plot 2: pursuit
    h2 = figure;
    set(gcf,'pos',[50 50 600 600],'color',cbac);


end


% TAB file conventions
%
% 1     2      3     4         5           6          7       8   9
% block trial3 alpha env_speed drift_speed trajLength initPos ecc fixOff 
% 
% 10   11   12      13   14   15    16     17     18     19      20      21         22
% tFix tBeg tFixOff tSac tEnd sacRT other1 other2 other3 tedfFix tedfTOn tedfFixOff tedfSac


%% prepare folder to store pursuit traces
pursuitFolder = ('../pursuit');
mkdir(pursuitFolder)


%% analysis

ntGoodAll = 0;
ntBadAll = 0;
vp = 0;
cnt = 1;
while cnt ~= 0  % per tutti i soggetti
    [vpcode, cnt] = fscanf(sfid,'%s',1);
    if cnt ~= 0
        vp = vp + 1;
        
        % create file strings
        tabfile = sprintf('../tab/%s.tab',vpcode);
        datfile = sprintf('../raw/%s.dat',vpcode);
        reafile = sprintf('../rea/%s_V%iD%iT%iS%i.rea',vpcode,velSD,minDur,VELTYPE,SAMPRATE);
        btrfile = sprintf('../btr/%s_V%iD%iT%iS%i.btr',vpcode,velSD,minDur,VELTYPE,SAMPRATE);
        
        if ~exist(tabfile,'file')
            fprintf(1,'\n\ttab file %s not found!',tabfile);
        end
        if ~exist(datfile,'file')
            fprintf(1,'\n\ttab file %s not found!',datfile);
        end
        
        fprintf(1,'\n\n\tloading %s ...',vpcode);
        
        tab = load(tabfile);
        dat = load(datfile);
        
        % load experimental design specifications
        subDir=substr(vpcode, 0, 4);
        sessionDir=substr(vpcode, 5, 2);
        designfile = sprintf('../../data/%s/%s/%s.mat',subDir,sessionDir,vpcode);
        load(designfile);
        
        fprintf(1,' preparing\n');
        frea = fopen(reafile,'w');
        fbtr = fopen(btrfile,'w');
        
        % counting variables
        nt   = 0;   % number of trials
        ntGood=0;   % number of good trials
        
        % exclusion criteria (for up to 100 response saccades)
        ex.sbrs = zeros(1,100);   % saccades > 1 deg before response saccade
        ex.mbrs = zeros(1,100);   % Missing before saccade (blinks)
        ex.nors = zeros(1,100);   % no response saccade
        ex.samp = zeros(1,  1);   % sampling rate error
        
        % empty pursuit traces
%         clear pursuitTraces
%         pursuitTraces.x = [];
%         pursuitTraces.y = [];
%         pursuitTracesFileNames = sprintf('%s/%s_pTrace.mat',pursuitFolder, vpcode);
%         ptcount = 0;
        
        for t = 1:size(tab,1) % per ogni trial
            nt = nt + 1;
            
            % load relevant variables
            block = tab(t,1);
            trial = tab(t,2);
            
            alpha = tab(t,3);
            ecc = tab(t,8);
            env_speed = tab(t,4);
            drift_speed = tab(t,5);
            trajLength = tab(t,6);
            initPos = tab(t,7);
            
            % load trial spec. and target trajectory
            td = design.b(block).trial(trial);

            tTarOnEDF = tab(t,20);
            tFixOnEDF = tab(t,19);
            tCueOnEDF = tab(t,21);    % saccadic cue (fixation offset)
            
            % nRS = size(tarPos,1)-1;     % nRS is 0 if no saccade required
            nRS = 1;
            
            % reset all rejection criteria
            if nRS > 0
                sbrs = zeros(1,nRS);
                mbrs = zeros(1,nRS);
                resp = zeros(1,nRS);
                samp = zeros(1,1);
            else
                sbrs = zeros(1,1);
                mbrs = zeros(1,1);
                resp = zeros(1,1);
                samp = zeros(1,1);
            end
            
            %% Primary response saccade analysis
            
            % get data in response time interval
            idxrs = find(dat(:,1)>=tCueOnEDF & dat(:,1)<=tCueOnEDF+maxRT);
            timers = dat(idxrs,1);	% time stamp
            
            % determine sampling rate (differs across
            % trials if something went wrong)
            if ~mod(length(timers),2)   % delete last sample if even number
                timers(end) = [];
                idxrs(end)  = [];
            end
            samrat = round(1000/mean(diff(timers)));
            minsam = minDur*samrat/1000;
            if 0 %samrat<SAMPRATE
                samp = 1;
                ex.samp = ex.samp+1;
            end
            
            xrsf = DPP*([dat(idxrs,2)-scrCen(1) -(dat(idxrs,3)-scrCen(2))]);    % positions
            
            % filter eye movement data
            clear xrs;
            clear mrs;
            xrs(:,1) = filtfilt(fir1(35,0.05*SAMPRATE/samrat),1,xrsf(:,1));
            xrs(:,2) = filtfilt(fir1(35,0.05*SAMPRATE/samrat),1,xrsf(:,2));
            
            vrs = vecvel(xrs, samrat, VELTYPE);    % velocities
            vrsf= vecvel(xrsf, samrat, VELTYPE);   % velocities
            mrs = microsaccMerge(xrsf,vrsf,velSD,minsam,mergeInt);  % saccades
            mrs = saccpar(mrs);
            if size(mrs,1)>0
                amp = mrs(:,7);
                mrs = mrs(amp>maxMSAmp,:);
            end
            nSac = size(mrs,1);
            
            %% target path and position at the time of the saccade
            if ~isempty(mrs)
            
                tarRadiusAtLanding = NaN;
                fixPos = [0 0];
                period = ((trajLength*2) / env_speed)*1000;
                travelTime  = timers(mrs(1,2))-tTarOnEDF;
                timeIndex_1 = travelTime./period; % mod(travelTime,(period*1000))/1000;
                timeIndex_2 = mod(travelTime+(1000/120),(period*1000))/1000;    % target position one refresh after saccade landing
                
                travelTime2  = timers(mrs(1,1))-tTarOnEDF; % this is saccade onset
                timeIndex_3 = travelTime2./period; % mod(travelTime2,(period*1000))/1000;
                timeIndex_4 = mod(travelTime2+(1000/120),(period*1000))/1000;    % target position one refresh after saccade landing
                
                switch design.motionType
                    case 'triangular'
                        tarRadiusAtLanding = initPos * (trajLength/2) * sawtooth(2*pi*timeIndex_1, 0.5);  % target radius for each frame
                        tarRadiusAtLanding_2 = initPos * (trajLength/2) * sawtooth(2*pi*timeIndex_2, 0.5);
                        tarRadiusAtSacOnset = initPos * (trajLength/2) * sawtooth(2*pi*timeIndex_3, 0.5);
                        tarRadiusAtSacOnset_2 = initPos * (trajLength/2) * sawtooth(2*pi*timeIndex_4, 0.5);
                    case 'sinusoidal'
                        tarRadiusAtLanding = initPos * (trajLength/2) * sin(2*pi*(timeIndex_1+0.25));
                        tarRadiusAtLanding_2 = initPos * (trajLength/2) * sin(2*pi*(timeIndex_2+0.25));
                        tarRadiusAtSacOnset = initPos * (trajLength/2) * sin(2*pi*(timeIndex_3+0.25));
                        tarRadiusAtSacOnset_2 = initPos * (trajLength/2) * sin(2*pi*(timeIndex_4+0.25));
                end
                
                [tarLand_x, tarLand_y] = pol2cart(deg2rad(90 + alpha), tarRadiusAtLanding); % target position at landing
                [tarLand_x_2, tarLand_y_2] = pol2cart(deg2rad(90 + alpha), tarRadiusAtLanding_2);
                [tarSacOn_x, tarSacOn_y] = pol2cart(deg2rad(90 + alpha), tarRadiusAtSacOnset);
                [tarSacOn_x_2, tarSacOn_y_2] = pol2cart(deg2rad(90 + alpha), tarRadiusAtSacOnset_2);
                
                %%% ADJIUST ALL POSITIONS FOR MONITOR COORDINATES (reverse y)
                tarSacOn_y = - tarSacOn_y;
                tarSacOn_y_2 = - tarSacOn_y_2;
                tarLand_y = - tarLand_y;
                tarLand_y_2 = - tarLand_y_2;
                
                % compute target envelope and internal motion at saccade onset, to be compared with post landing pursuit direction
                envDirAng = cart2pol(tarSacOn_x - tarLand_x, tarSacOn_y - tarLand_y);
                driftDirAng = envDirAng + sign(alpha)*(pi/2); % cambiato segno 20-2-2015
                
                tarPosLand = [tarLand_x + ecc, tarLand_y];
                tarPosSacOn = [tarSacOn_x + ecc, tarSacOn_y];
                tarPos = [fixPos; [tarLand_x + ecc, tarLand_y]];
                tarLand_x = tarLand_x + ecc;
                tarSacOn_x = tarSacOn_x + ecc;
                
                % determine also the direction of carrier drift within the
                % gabor
                driftDir = NaN;
                timeToSaccOnset = timers(mrs(1,1))-tTarOnEDF;
                timeIndex_son_1 = mod(timeToSaccOnset-(1000/120),(period*1000))/1000;
                timeIndex_son_2 = mod(timeToSaccOnset,(period*1000))/1000;
                switch design.motionType
                    % here we are interested only in the sign of target
                    % phase variation at the time of saccade onset
                    case 'triangular'
                        tarPhase_1 = initPos*(-sign(alpha)) * sawtooth(2*pi*timeIndex_son_1, 0.5);
                        tarPhase_2 = initPos*(-sign(alpha)) * sawtooth(2*pi*timeIndex_son_2, 0.5);
                    case 'sinusoidal'
                        tarPhase_1 = initPos*(-sign(alpha)) * -cos(2*pi*(timeIndex_son_1+0.25));
                        tarPhase_2 = initPos*(-sign(alpha)) * -cos(2*pi*(timeIndex_son_2+0.25));
                end
                % note that the direction of phase changes can be zero, es.
                % at the point of reversal or if target has no internal
                % motion
                driftDir = sign(tarPhase_2 - tarPhase_1);
                
            else
                tarPos = [0 0; ecc 0];
            end
            
            % target trajectory (for plotting)
            [tar1end_x, tar1end_y] = pol2cart(deg2rad(90 + alpha), -(trajLength/2));
            [tar2end_x, tar2end_y] = pol2cart(deg2rad(90 + alpha),  (trajLength/2));
            tar1end_x = tar1end_x + ecc;
            tar2end_x = tar2end_x + ecc;
         
            tar2end_y = - tar2end_y; % 20-2-2015
            tar1end_y = - tar1end_y;
            
            %% Analysis of post-landing pursuit
            
            % define output variables
            pursuitDirAng = NaN;
            pos_p = [NaN NaN];
            vel_p = [NaN NaN];
            acc_p = [NaN NaN];
            clear posF;
            clear velF;
            clear accF;
            
            if size(mrs,1)>=1   % if there has been a saccade (here i look only after the first one)
                
                pursuitOK = 0;
                
                if (mrs(1,2)+pursuitInterval(2)) < size(xrsf,1)  % check if pursuit trace is enough for analysis
                    
                    pursuitOK = 1;
                    
                    if addFilter %% zero-phase filtering ---------------------------------------------
                        
                        if low_pass_pos ~= 0    % --- Filter position
                            cutv = low_pass_pos/(SAMPRATE/2);
                            PXf = filtfilt(fir1(35,cutv),1,xrsf(:,1));
                            PYf = filtfilt(fir1(35,cutv),1,xrsf(:,2));
                            posF = [PXf PYf];
                        else
                            posF = xrsf;
                        end
                        
                        % determine velocity
                        vel = vecvel(posF, samrat, VELTYPE);
                        
                        if low_pass_vel ~= 0    % --- Filter velocity
                            cutv = low_pass_vel/(SAMPRATE/2);
                            VXf = filtfilt(fir1(35,cutv),1,vel(:,1));
                            VYf = filtfilt(fir1(35,cutv),1,vel(:,2));
                            velF = [VXf VYf];
                        else
                            velF = vel;
                        end
                        
                        % determine acc
                        acc = vecacc(velF); % acc=[diff(velF(:,1)) diff(velF(:,1))];
                        
                        if low_pass_acc ~= 0    % --- Filter accelleration
                            cutv = low_pass_acc/(SAMPRATE/2);
                            accXf = filtfilt(fir1(35,cutv),1,acc(:,1));
                            accYf = filtfilt(fir1(35,cutv),1,acc(:,2));
                            accF = [accXf accYf];
                        else
                            accF = acc;
                        end
                        
                    else
                        posF = xrsf;
                        velF=[vrsf(:,1) vrsf(:,2)];
                        accF = vecacc(velF);    % accF = [diff(vrsf(:,1)) diff(vrsf(:,2))];
                    end
                    
                end %% end filter part ---------------------------------------------------------------
                
                % check available data
                if (mrs(1,2) + pursuitInterval(2)) < size(xrsf,1)
                    actPursuitInterval = [mrs(1,2)+pursuitInterval(1) : mrs(1,2)+pursuitInterval(2)];
                else
                    actPursuitInterval = [mrs(1,2)+pursuitInterval(1) : size(xrsf,1)-1];
                end
                
                % corrective saccade to NaN values
                if size(mrs,1)>1 && (actPursuitInterval(end)) >= mrs(2,1)
                    corrSacc = mrs(2,1):mrs(2,2);
                end
                
                if length(actPursuitInterval) > minPursuitDur && exist('posF','var')
                    
                    % analysis: NB. FILTERED DATA
                    posPursuit = posF(actPursuitInterval, :);
                    velPursuit = velF(actPursuitInterval, :);
                    accPursuit = accF(actPursuitInterval(1:end-1), :);
                    
                    if exist('coorSacc','var')
                        posPursuit(corrSacc) = NaN;
                        velPursuit(corrSacc) = NaN;
                        accPursuit(corrSacc) = NaN;
                    end
                    
                    %rel_p_x = posPursuit(:,1) - posPursuit(1,1);
                    %rel_p_y = posPursuit(:,2) - posPursuit(1,2);
                    
                    rel_p_x = posPursuit(:,1) - xrs(mrs(1,2),1);
                    rel_p_y = posPursuit(:,2) - xrs(mrs(1,2),2);
                    
                    % transform in polar coordinates for plotting
                    [theta_p, rho_p] = cart2pol(rel_p_x, rel_p_y);
                    
                    % direction angle
                    % pursuitDirAng = compDirection(rel_p_x, rel_p_y);  % default: average between the first and third part
                    % sacxOffset = xrs(mrs(reaSacNumber,2),1);
                    % sacyOffset = xrs(mrs(reaSacNumber,2),2);
                    pursuitDirAng = cart2pol(mean(posPursuit(:,1))-xrs(mrs(1,2),1), mean(posPursuit(:,2))-xrs(mrs(1,2),2));
                    
                    % horizontal and vertical components (position)
                    pos_p = [(posPursuit(end,1) - posPursuit(1,1)), (posPursuit(end,2) - posPursuit(1,2))];
                    
                    % peak horizontal and vertical velocities
                    vel_p = [absMax(velPursuit(:,1)), absMax(velPursuit(:,2))];
                    
                    % peak accellerations
                    acc_p = [absMax(accPursuit(:,1)), absMax(accPursuit(:,2))];
            
                end
                
            end
            
            pursuitRes = [envDirAng driftDirAng pursuitDirAng pos_p vel_p acc_p];
            % format string: %.2f/t%.2f/t%.2f/t%.2f/t%.2f/t%.2f/t%.2f/t%.2f/t%.2f/t
            
            
            %% PLOT TRACES
            
            if plotData
                
                figure(h1);
                
                axes(ax(1));
                hold off;
                
                plot(xrsf(:,1),xrsf(:,2),'k-','color',[0.5 0.5 0.5]);
                hold on;
                
                % fixation area
                plot(tarPos(1,1),tarPos(1,2),'ko','color',[0.5 0.5 0.5],'MarkerFaceColor',[1 1 1],'markersize',(2.6*PPD)/2,'LineWidth',1);
                
                % trajectory
                if drift_speed ~= 0
                    line([ecc ecc],[-(trajLength/2) (trajLength/2)],'color',[0.1 0.1 0.1],'LineStyle','--','LineWidth',2)
                end
                line([tar1end_x tar2end_x],[tar1end_y tar2end_y],'color',[0.9 0.9 0.9],'LineWidth',9)
                
                if ~isempty(mrs)
                    plot(tarLand_x,tarLand_y,'.','color',[0 0 0],'markersize',40);
                    plot(tarSacOn_x,tarSacOn_y,'.','color',[0.8 0.8 0.8],'markersize',20);
                end
                
                for i = 1:size(mrs,1)
                    plot(xrsf(mrs(i,1):mrs(i,2),1),xrsf(mrs(i,1):mrs(i,2),2),'r-','linewidth',2);
                end
                if ~isempty(mrs)
                    plot([xrsf(mrs(:,1),1) xrsf(mrs(:,2),1)]',[xrsf(mrs(:,1),2) xrsf(mrs(:,2),2)]','r.');
                end
                
                axes(ax(2));
                plot(timers-tCueOnEDF,vrsf(:,1),'k-','color',[0.6 0.6 0.6]);
                hold on;
                for i = 1:size(mrs,1)
                    plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrsf(mrs(i,1):mrs(i,2),1)','r-','linewidth',2);
                end
                hold off;
                xlim(timers([1 end])-tCueOnEDF);
                ylim([-1000 1000]);
                
                axes(ax(3));
                plot(timers-tCueOnEDF,vrsf(:,2),'k-','color',[0.6 0.6 0.6]);
                hold on;
                for i = 1:size(mrs,1)
                    plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrsf(mrs(i,1):mrs(i,2),2)','r-','linewidth',2);
                end
                hold off;
                xlim(timers([1 end])-tCueOnEDF);
                ylim([-1000 1000]);
                
                
                if pursuitOK % plot post saccadic trace --------------------------------------------------------
                
                    figure(h2);
                        
                    polar(0, 0.4);
                    
                    hold on;
                    
                    % pursuit trace
                    tp = polar(theta_p, rho_p);
                    set(tp,'Color',[0.5 0.5 0.5],'LineWidth',3);
                    
                    % linear fit
                    tp1 = polar([pursuitDirAng pursuitDirAng], [0 0.4],'--');
                    set(tp1,'Color',[0.5 0.5 0.5],'LineWidth',2);
                    
                    tp2 = polar([envDirAng envDirAng], [0 0.4],'-.');
                    set(tp2,'Color',[0 0.5 0],'LineWidth',2);
                    
                    if drift_speed ~= 0
                        tp3 = polar([driftDirAng driftDirAng], [0 0.4],'-.');
                        set(tp3,'Color',[0 0 0.5],'LineWidth',2);
                        legend([tp1,tp2,tp3],'gaze','envelope','carrier','Location','Best');
                    else
                        legend([tp1,tp2],'gaze','envelope','Location','Best');
                    end
                    
                    
                    hold off;
                    
                end % end plot post saccadic trace -------------------------------------------------------------
                
            end
            
            %%
            % loop over all necessary response saccades
            s = 0;
            for rs = 1:nRS
                fixRec = repmat(tarPos(rs  ,:),1,2)+[-fixRad -fixRad fixRad fixRad];
                tarRec = repmat(tarPos(rs+1,:),1,2)+[-tarRad -tarRad tarRad tarRad];
                
                % check for response saccade
                reaSacNumber = 0;
                if s<nSac
                    while ~resp(rs) && s<nSac
                        s = s+1;
                        onset = timers(mrs(s,1))-tCueOnEDF;
                        xBeg  = xrs(mrs(s,1),1);    % initial eye position x
                        yBeg  = xrs(mrs(s,1),2);	% initial eye position y
                        xEnd  = xrs(mrs(s,2),1);    % final eye position x
                        yEnd  = xrs(mrs(s,2),2);	% final eye position y
                        
                        % is saccade out of fixation area?
                        if ~resp(rs)
                            fixedFix = isincircle(xBeg,yBeg,fixRec);
                            fixedTar = isincircle(xEnd,yEnd,tarRec);
                            if fixedTar && fixedFix
                                reaSacNumber = s;   % which saccade after cue went to target
                                reaLoc = fixedTar;  % which target location did saccade go to (if multiple)
                                resp(rs) = 1;
                            end
                        end
                        
                        %%
                        if plotData
                            axes(ax(1));
                            if resp(rs)
                                plot([xBeg xEnd],[yBeg yEnd],'ko','color',[0 0 0],'markersize',8,'linewidth',2);
                            end
                        end
                        %%
                    end
                end
                
                if ~resp(rs)
                    ex.nors(rs) = ex.nors(rs) + 1;
                    
                    sacType    = 0;    % 0 = no response; 1 = microsaccade; 2 = large saccade; 3 = no saccade task
                    sacOnset   = NaN;
                    sacOffset  = NaN;
                    sacDur     = NaN;
                    sacVPeak   = NaN;
                    sacDist    = NaN;
                    sacAngle1  = NaN;
                    sacAmp     = NaN;
                    sacAngle2  = NaN;
                    sacxOnset  = NaN;
                    sacyOnset  = NaN;
                    sacxOffset = NaN;
                    sacyOffset = NaN;
                    sacRT      = NaN;
                else
                    % compile reaSac data
                    sacType    = rs+1; % 0 = no response; 1 = microsaccade; 2 = large saccade; 3 = second large saccade
                    sacOnset   = timers(mrs(reaSacNumber,1))-tCueOnEDF; % start-screen-geloggt
                    sacOffset  = timers(mrs(reaSacNumber,2))-tCueOnEDF; % start-screen-geloggt
                    sacDur     = mrs(reaSacNumber,3)*1000/samrat;
                    sacVPeak   = mrs(reaSacNumber,4);
                    sacDist    = mrs(reaSacNumber,5);
                    sacAngle1  = mrs(reaSacNumber,6);
                    sacAmp     = mrs(reaSacNumber,7);
                    sacAngle2  = mrs(reaSacNumber,8);
                    sacxOnset  = xrs(mrs(reaSacNumber,1),1);
                    sacyOnset  = xrs(mrs(reaSacNumber,1),2);
                    sacxOffset = xrs(mrs(reaSacNumber,2),1);
                    sacyOffset = xrs(mrs(reaSacNumber,2),2);
                    
                    if rs>1
                        sacRT = sacOnset-reaSac(rs-1,5);
                    else
                        sacRT = sacOnset;
                    end
                end
                
                reaSac(rs,:) = [sacRT reaSacNumber sacType sacOnset sacOffset sacDur, ...
                    sacVPeak sacDist sacAngle1 sacAmp sacAngle2 sacxOnset, ...
                    sacyOnset sacxOffset sacyOffset];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % saccades before current response saccade %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if reaSacNumber>rs
                    sbrs(rs) = 1;
                    ex.sbrs(rs) = ex.sbrs(rs)+1;
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % missings before current response saccade %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if resp(rs)
                    rt = sacRT;
                    if rs == 1
                        t1 = tCueOnEDF-timBefCue;
                        t2 = tCueOnEDF+sacRT+sacDur;
                    else
                        t1 = tCueOnEDF+reaSac(rs-1,5);
                        t2 = t1+sacRT+sacDur;
                    end
                    
                    % get index for current data period
                    idx  = find(dat(:,1)>=t1 & dat(:,1)<t2);
                    
                    % check for missing data
                    idxmbrs = find(dat(idx,crit_cols)==-1);
                    if ~isempty(idxmbrs)
                        mbrs(rs) = 1;
                        ex.mbrs(rs) = ex.mbrs(rs) + 1;
                    end
                end
            end
            
            %%
            if plotData
                axes(ax(1));
                hold on;
                xlim([-3 13]);
                ylim([-8 8]);
                text(-2,6,sprintf('%s',substr(vpcode, 3, 2)),'Fontsize',14,'hor','left','ver','bottom');
                text(-2,-6,sprintf('Trial %i (b. %i)',trial,block),'Fontsize',14,'hor','left','ver','bottom');
                
                condStr = 'Saccade';
                
                if sum(sbrs)==0 && sum(mbrs)==0 && sum(resp)==nRS % && samp==0
                    text(4,-6,sprintf('%s: good',condStr),'Fontsize',14,'hor','center','ver','bottom');
                else
                    text(4,-6,sprintf('%s: bad',condStr),'Fontsize',14,'hor','center','ver','bottom');
                end
                
                text(10,-6,sprintf('alpha= %.2f',alpha),'Fontsize',14,'hor','center','ver','bottom');
                text(10, 6,sprintf('drift = %.2f',drift_speed),'Fontsize',14,'hor','center','ver','bottom');
                
                set(gca,'visible','off');
                
                % export image
                if printImages
                    o = input('\n   Save image? [y / n]? ','s');
                    if strcmp(o,'y')
                        exportfig(h1,sprintf('%s_trial%i.%i.eps',substr(vpcode, 3, 2), block, trial),'bounds','tight','color','rgb','LockAxes',0);
                    end
                end
                
                if waitforbuttonpress
                    plotData = 1;
                    % pause;
                    
                end
                cla;
            end
            %% end % this would be the end for if ~samp
            
            % evaluate rejection criteria
            if sum(sbrs)==0 && sum(mbrs)==0 && sum(resp)==nRS % && samp==0
                
                % additional parameters for analyses
                % \t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f
                added = [driftDir pursuitRes tarPosLand tarPosSacOn];
                
                % output vector
                rea = [repmat(tab(t,:),nRS,1) reaSac added];
                
                % tab format
                % %i\t%i\t.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\n
                
                % this is what you need in addition to the tab format:
                % %i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\n
                              %i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n
                fprintf(frea,'%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',rea);
                                                                                                                                                                                                                   
                if welche == 1
                    
                    fprintf(fRAll,  '%s\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',vpcode,vp,str2num(sessionDir),rea);
                    
                    rea = [vp*ones(size(rea,1),1), rea added];
                                     %i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n
                    fprintf(freaAll,'%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n',rea);
                end
                ntGood = ntGood + 1;
                ntGoodAll = ntGoodAll + 1;
            else
                ntBadAll = ntBadAll + 1;
            end
            fprintf(fbtr,'%i\t%i\t%i\t%i\t%i\t%i\n',tab(t,1),tab(t,2),sbrs,mbrs,~resp,samp);
        end
        fclose(frea);
        fclose(fbtr);
        
        % save pursuit traces matrix for a single subject
%         save(pursuitTracesFileNames,'-struct','pursuitTraces');
        
        fprintf(1,'\t Tot. subject trials:\t\t%i\n\t GoodTrials :\t\t%i\n\t Samp. error:\t\t%i\n',nt,ntGood,ex.samp);
        for rs = 1:nRS
            fprintf(1,'\t Response saccade %i:\n\t  SacBefCue :\t\t%i\n\t  MisBefSac :\t\t%i\n\t  NoSaccade :\t\t%i\n',rs,ex.sbrs(rs),ex.mbrs(rs),ex.nors(rs));
        end
        
        % save pursuit image ----------------------------------------------
%         exportfig(h3,sprintf('%s.eps',vpcode),'Color','rgb');
%         figure(h3); 
%         subplot(1,2,1); hold off
%         subplot(1,2,2); hold off
        % -----------------------------------------------------------------
        
        if welche == 1
            fprintf(fbtrAll,'%s\t%i\t%i\t%i',vpcode,nt,ntGood,ex.samp);
            for rs = 1:nRS
                fprintf(fbtrAll,'\t%i\t%i\t%i',ex.sbrs(rs),ex.mbrs(rs),ex.nors(rs));
            end
            fprintf(fbtrAll,'\n');
        end
    end
end
fclose(sfid);

if welche == 1
    fclose(fbtrAll);
    fclose(freaAll);
    fclose(fRAll);
end

fprintf(1,'\n\nTotal came %i good trials out of %i total trials!',ntGoodAll,ntGoodAll+ntBadAll);
