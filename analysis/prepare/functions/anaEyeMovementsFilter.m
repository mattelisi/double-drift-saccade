%%
%
% analysis - saccade task (saccade to moving target - sampled motion)
% 
% [in questa versione il target diventa verde, e scompare appena il
% soggetto fa la saccade]
%
% - analizza anche il pursuit post-saccadico
% 
% Matteo Lisi, 2013
%
% updates: 13/06/2013 - updated exclusion of corrective saccades from post-landing pursuit response.
%                       (avoiding error of different column numbers in the output in case of shorter analyzed intervals)
%
%



if welche == 1
    reaAllFile = sprintf('../rea/AllSubjects_V%iD%iT%iS%i.rea',velSD,minDur,VELTYPE,SAMPRATE);
    reaAllFile_R = sprintf('../R/AllSubjects_V%iD%iT%iS%i.txt',velSD,minDur,VELTYPE,SAMPRATE);
    btrAllFile = sprintf('../btr/AllSubjects_V%iD%iT%iS%i.btr',velSD,minDur,VELTYPE,SAMPRATE);
    
    freaAll = fopen(reaAllFile,'w');
    fbtrAll = fopen(btrAllFile,'w');
                    
%      fprintf(freaAll, 'subject\tblock\ttrial\ttargetAngle\ttarDir\tnStepPre\tnStepPost\tsigma\tecc\td\tangVel\ttanVel\t'); % 1 - 12
%      fprintf(freaAll, 'tFix\ttBeg\ttChange\ttSac\ttEnd\ttClr\tsacRTno\ttedfFix\ttedfTOn\ttedfSac\ttedfChange\ttedfTOf\ttedfClr\t'); % 13 - 25
%      fprintf(freaAll, 'sacRT\treaSacnumber\tsacType\tsacOnset\tsacOffset\tsacDur\tsacVPeak\tsacDist\tsacAngle1\tsacAmp\tsacAngle2\t'); % 26 - 36
%      fprintf(freaAll, 'sacxOnset\tsacyOnset\tsacxOffset\tsacyOffset\t'); % 37 - 40
%      fprintf(freaAll, 'tarDirAngle\tvelDev\tpursuitVel\tpursuitVelAngle\tposDiffDev\tpursuitPosDiff\tpursuitPosDiffAngle\n'); % 41 - 47
%      fprintf(freaAll, 'tarLandDirAngle\ttarLand\ttarLandPosDirAngle\tpursuitAcc'); % 48 - 51
%      fprintf(freaAll, 'dT'); % 52
    
    sfid = fopen('subjects.all','r');
else
    sfid = fopen('subjects.tmp','r');
end

%%
if plotData
    sl = 0.5;   % box side length
    fp = 0.6;   % fixation point width
    ec = 15.0;  % eccentricity
    cl = 0.2;   % cue length
    
    cout = [0.3 0.3 0.3];
    cbac = [1.0 1.0 1.0];
    
    close all;
    h1 = figure;
%     set(gcf,'pos',[50 50 1000 500],'color',cbac);
%     ax(1) = axes('pos',[0 0 0.5 1]);
%     ax(2) = axes('pos',[0.5 0.5 0.5 0.5]);
%     ax(3) = axes('pos',[0.5 0.0 0.5 0.5]);
    set(gcf,'pos',[50 50 1500 500],'color',cbac);
    % left bottom width height
    ax(1) = axes('pos',[0 0 0.33333 1]);
    ax(2) = axes('pos',[0.35 0.5 0.3 0.5]);
    ax(3) = axes('pos',[0.35 0.0 0.3 0.5]);
    ax(4) = axes('pos',[0.65 0.0 0.35 1]);
    
    
    % plot 2: pursuit
    h2 = figure;
    set(gcf,'pos',[50 50 600 600],'color',cbac);
    
    % plot 3: angles
%     h3 = figure;
%     set(gcf,'pos',[50 50 1000 500],'color',cbac);

    % plot deviations
     h5 = figure;
     set(gcf,'pos',[50 50 400 400],'color',cbac);
    
% else
%     % plot 3: angles
%     h3 = figure;
%     cbac = [1.0 1.0 1.0];
%     set(gcf,'pos',[50 50 1000 500],'color',cbac);

end

%
% TAB file
%
% 1     2      3   4   5        6         7     8   9 10     11     12
% block trial3 tar dir nStepPre nStepPost sigma ecc d angVel tanVel tFix
%
% 13   14      15   16   17   18    19      20      21      22         23      24      25
% tBeg tChange tSac tEnd tClr sacRT tedfFix tedfTOn tedfSac tedfChange tedfTOf tedfClr lastStep
%

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
        
        for t = 1:size(tab,1) % per ogni trial
            nt = nt + 1;
            
            % load relevant variables
            block = tab(t,1);
            trial = tab(t,2);
            dir = tab(t,4);
            ecc = tab(t,8);
            angVel = tab(t,10);
            lastStep = tab(t,25);
            
            % load trial spec. and target trajectory
            td = design.b(block).trial(trial);
            
            path =  detPathCircular(td.tar, td.dir, td.nStepPre, td.nStepPost, td.beta, td.ecc, visual);
            
            for pathStep = 1:size(path,1)
                path(pathStep,:)=DPP*[path(pathStep,1)-scrCen(1),-(path(pathStep,2)-scrCen(2))];
            end
            tarIndex = find(td.events==2);  
            
            if ~isnan(lastStep)
                path = path(1:lastStep,:);
            end
            
            tFixOnEDF = tab(t,19);
            tCueOnEDF = tFixOnEDF+tab(t,22);    % saccadic cue (flash or color change)
            
            % fixation and target position
            fixPos = [0 0];
            tarPos = path(tarIndex,:);
            
            tar = tab(t,3);
            
            tarPos = [fixPos; tarPos];
            
            nRS = size(tarPos,1)-1;     % nRS is 0 if no saccade required
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Analysis of saccadic response times %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            
            
            if size(mrs,1)>=1
            
                % pre-define pursuit variables
                tarDirAngle     = NaN;
                velDev          = NaN;
                pursuitVel      = NaN;
                pursuitVelAngle = NaN;
                posDiffDev      = NaN;
                pursuitPosDiff  = NaN;
                pursuitPosDiffAngle = NaN;
                tarLandDirAngle = NaN;
                tarLand         = NaN;
                tarLandPosDirAngle = NaN;
                pursuitAcc = NaN;
                
                
                %% PURSUIT ANALYSIS
                %
                % this check pursuit ONLY AFTER THE FIRST RESPONSE SACCADE
                % prima di tutto aggiungere un check per veder se c'� il
                % tracciato o mancano dei dati (es. se chiude gli occhi subito dopo la saccade)
                %

                nsamples = size(xrsf,1);

                % calculate target direction at landing time...
                % then linspace tra onset e landing
                % [tar tarTime landtime velocity (angular)]
                
                travelTime = timers(mrs(1,2))-tCueOnEDF; % dopo la prima saccade
                tarLand = tar - dir * deg2rad(angVel) *(travelTime/1000);
                
                % target [x, y] position at landing
                tarPosLand = [ecc*cos(tarLand), ecc*sin(tarLand)];
                
                % target direction angle
                tarDirAngle = tar - dir*(pi/2);
                
                % target direction angle at landing time
                tarLandDirAngle = tarLand - dir*(pi/2);
                
                % CALCOLO ANCHE DIREZIONE AT LANDING POSITION
                sacAngle = mean([mrs(1,6), mrs(1,8)]);
                if sacAngle<0; sacAngle = 2*pi + sacAngle; end
                tarLandPosDirAngle = sacAngle - dir*(pi/2);
                
                % check if pursuit trace is enough for analysis
                pursuitOK = 0;
                
                if (mrs(1,2)+pursuitInterval(2))<=size(xrsf,1)
                    
                    pursuitOK = 1;
                    
                    if addFilter
                        
                        % new: added zero-phase filtering
                        
                        % --- Filter position
                        if low_pass_pos ~= 0
                            % [PXf PYf]=FiltrSacc2D(xrsf',nsamples,order,low_pass_pos, SAMPRATE);
                            cutv = low_pass_pos/(SAMPRATE/2);
                            PXf = filtfilt(fir1(35,cutv),1,xrsf(:,1));
                            PYf = filtfilt(fir1(35,cutv),1,xrsf(:,2));
                            posF = [PXf PYf];
                        else
                            posF = xrsf;
                        end
                        
                        % determine velocity
                        vel = vecvel(posF, samrat, VELTYPE); 
                        
                        % --- Filter velocity
                        if low_pass_vel ~= 0
                            % [VXf VYf]=FiltrSacc2D(vrsf',nsamples,order,low_pass_vel, SAMPRATE);
                            cutv = low_pass_vel/(SAMPRATE/2);
                            VXf = filtfilt(fir1(35,cutv),1,vel(:,1));
                            VYf = filtfilt(fir1(35,cutv),1,vel(:,2));
                            velF = [VXf VYf];
                        else
                            % velF=[vrsf(:,1) vrsf(:,2)];
                            velF = vel;
                        end
                        
                        % determine acc
                        acc = vecacc(velF); 
                        % acc=[diff(velF(:,1)) diff(velF(:,1))];
                        
                        % --- Filter accelleration
                        if low_pass_acc ~= 0
                            % [accXf accYf]=FiltrSacc2D(acc',nsamples-1,order,low_pass_acc, SAMPRATE);
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
                        % accF = [diff(vrsf(:,1)) diff(vrsf(:,2))];
                        accF = vecacc(velF); 
                    end
                    
                    % exclude corrective saccades from pursuit analysis
%                     actPursuitInterval = [mrs(1,2)+pursuitInterval(1) : mrs(1,2)+pursuitInterval(2)];
%                     if size(mrs,1)>1 && (mrs(1,2)+pursuitInterval(2)) >= mrs(2,1)
%                         if mrs(2,2) < (mrs(1,2)+pursuitInterval(2))
%                             corrSacc = mrs(2,1):mrs(2,2);
%                         else
%                             corrSacc = mrs(2,1):(mrs(1,2)+pursuitInterval(2));
%                         end
%                         
%                         actPursuitInterval = setdiff(actPursuitInterval,corrSacc);
%                     end
                    
                    % corrective saccade to NaN values
                    actPursuitInterval = [mrs(1,2)+pursuitInterval(1) : mrs(1,2)+pursuitInterval(2)];
                    if size(mrs,1)>1 && (mrs(1,2)+pursuitInterval(2)) >= mrs(2,1)
                        corrSacc = mrs(2,1):mrs(2,2);
                        % actPursuitInterval(corrSacc) = NaN;
                    end

                    
                    if length(actPursuitInterval) > minPursuitDur
                    
                        % analysis: NB. FILTERED DATA
                        posPursuit = posF(actPursuitInterval, :);
                        velPursuit = velF(actPursuitInterval, :);
                        accPursuit = accF(actPursuitInterval(1:end-1), :);
                        
                        if exist('coorSacc')
                            posPursuit(corrSacc) = NaN;
                            velPursuit(corrSacc) = NaN;
                            accPursuit(corrSacc) = NaN;
                        end
                            
                        rel_p_x = posPursuit(:,1) - posPursuit(1,1);
                        rel_p_y = posPursuit(:,2) - posPursuit(1,2);
                        [theta_p, rho_p] = cart2pol(rel_p_x, rel_p_y);
                        
                        % 1) velocities vector
                        [theta_v, rho_v] = cart2pol(velPursuit(:,1), velPursuit(:,2));
                        pursuitVel = mean(rho_v);
                        % pursuitVel = max(rho_v);
                        % pursuitVel = mean(rho_v(round(length(rho_v)/2):end)); % se vuoi considerare solo la seconda met�
                        pursuitVelAngle = mean(theta_v);
                        
                        % 2) position difference vector
                        pursuitPosDiff = sum(diff(rho_p));
                        % pursuitPosDiff = rho_p(end);
                        pursuitPosDiffAngle = mean(theta_p);
                        
                        % 3) accelleration
                        [theta_acc, rho_acc] = cart2pol(accPursuit(:,1), accPursuit(:,2));
                        pursuitAcc = max(rho_acc);
                        % pursuitAcc = (mean(accPursuit(:,1))+mean(accPursuit(:,1)))/2;
                        
                        % calcolo deviazioni dalla direzione del target
                        posDiffDev = angDiff(tarLandDirAngle, pursuitPosDiffAngle);
                        velDev = angDiff(tarLandDirAngle, pursuitVelAngle);
                    
                        pursuitRes = [tarDirAngle, velDev, pursuitVel, pursuitVelAngle, posDiffDev, pursuitPosDiff, pursuitPosDiffAngle, tarLandDirAngle, tarLand, tarLandPosDirAngle, pursuitAcc];
                        
                        if plotData
%                             figure(h3);
%                             
%                             LWD = 1;
%                             if plotData
%                                 LWD = 2;
%                             end
%                             
%                             subplot(1,2,1); if plotData; hold off; end
%                             % polar(theta,rho); % theta is the angle
%                             polar(0, 28.5); hold on
%                             p = polar([velDev velDev], [0 pursuitVel]);
%                             set(p,'Color','red','LineWidth',LWD);
%                             title('Post-saccadic velocity vectors (deg./sec)')
%                             
%                             subplot(1,2,2); if plotData; hold off; end
%                             polar(0, 2.85); hold on
%                             p = polar([posDiffDev posDiffDev], [0 pursuitPosDiff],'-r');
%                             set(p,'Color','red','LineWidth',LWD);
%                             title('Post-saccadic position-drift vectors (deg.)')
                            
                            
                            % plot pursuit trace
                            figure(h2);
                            
                            tarTanVel = angVel/ecc;
                            xRange = [1 pursuitInterval(2)];
                            tarTrace = [xRange(1):xRange(2)]*(tarTanVel/1000);
                            
                            % calculate expected target path
                            timeInd = xRange(1):xRange(2);
                            tarPathAngles = tarLand - dir .* deg2rad(angVel) .*(timeInd./1000);
                            tarPath = [ecc.*cos(tarPathAngles); ecc.*sin(tarPathAngles)]';
                            % tarVel = 1000*vecacc(tarPath);
                            % tarAcc = 1000*vecacc(tarVel);
                            tarVel = 1000*[diff(tarPath(:,1)) diff(tarPath(:,2))];
                            tarAcc = 1000*[diff(tarVel(:,1)) diff(tarVel(:,2))];
                            
                            subplot(2,1,1)
                            plot(tarPath,'--','LineWidth',1); hold on
                            pl = plot(posPursuit,'LineWidth',2); hold off
                            xlim(xRange); 
                            ylabel('distance [deg]')
                            legend(pl,'horizontal','vertical','Location','NorthEast');
                            
                            subplot(2,1,2)
                            plot(tarVel,'--','LineWidth',1); hold on
                            plot(velPursuit,'LineWidth',2); hold off
                            xlim(xRange); 
                            ylabel('velocity [deg/sec]')
                            % legend(pl,'horizontal','vertical','Location','NorthEastOutside');
                            
%                             subplot(3,1,3)
%                             plot(tarAcc,'--','LineWidth',1); hold on
%                             plot(accPursuit,'LineWidth',2); hold off
%                             xlim(xRange);
%                             ylabel('acceleration [deg/sec/sec]')
%                             legend(pl,'horizontal','vertical','Location','NorthEastOutside');
                            
                        end
                    
                    else
                        pursuitRes = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
                    end
                    
                else
                    pursuitRes = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
                    
                end
            end
            
            %% Saccade deviation
            if size(mrs,1)>=1
                X=xrsf(mrs(1,1):mrs(1,2),1);
                Y=xrsf(mrs(1,1):mrs(1,2),2);
                if ~plotData
                    [curvature, R2] = saccDeviation(X,Y);
                    curvatureData = [curvature, R2];
                else
                    [curvature, R2, curvatureFit] = saccDeviation(X,Y);
                    curvatureData = [curvature, R2];
                    
                    % plot curvature
                    figure(h5);
                    
                    plot(curvatureFit.Xh_n, curvatureFit.Yh,'color',[0.6 0.6 0.6],'linewidth',2);    hold on
                    plot(curvatureFit.Xh_n, curvatureFit.predicted,'color',[0 0 0],'linewidth',1);  hold off
                    xlim([-1.5 1.5]);
                    ylim([-7 7]);
                    ylabel('curvature [deg]')
                    xlabel('normalized amplitude')
                    
                    % legend(pl,num2str(dTs),'Location','NorthEastOutside');
                    
                    if dir == 1 % 1 = CLOCKWISE; -1 = COUNTERCLOCKWISE
                        text(-1.3,5,'target motion: clockwise','Fontsize',14,'hor','left','ver','bottom');
                    else
                        text(-1.3,5,'target motion: counterclockwise','Fontsize',14,'hor','left','ver','bottom');
                    end
                    text(-1.2,-6,sprintf('quadratic curvature: %.2f \t\t (R-squared: %.2f)',curvature, R2),'Fontsize',10,'hor','left','ver','bottom');
                    
                end
            else
                curvatureData = [NaN NaN];
            end
            
            %% PLOT TRACES
            
            if plotData
                
                figure(h1);
                
                axes(ax(1));
                hold off;
                
                % sxrs = smoothdata(xrs);
                % plot(sxrs(:,1),sxrs(:,2),'k-','color',[0.5 0.5 0.5]);
                plot(xrsf(:,1),xrsf(:,2),'k-','color',[0.5 0.5 0.5]);
                hold on;
                
                % big circle
                % plot(0,0,'ko','color',[0.2 0.2 0.2],'markersize',(10*PPD),'linewidth',1);
                
                plot(tarPos(1,1),tarPos(1,2),'ko','color',[0.5 0.5 0.5],'MarkerFaceColor',[0.9 0.9 0.9],'markersize',(2.6*PPD)/2,'linewidth',1);
                
                % trajectory
                for pathStep = 1:size(path,1)
                    if pathStep<= tarIndex
                        plot(path(pathStep,1),path(pathStep,2),'.','color',[0.9 0.9 0.9],'markersize',20);
                    else
                        plot(path(pathStep,1),path(pathStep,2),'.','color',[0 1 0],'markersize',20);
                    end
                end
                
                
                plot(tarPosLand(1),tarPosLand(2),'ko','color',[0.6 0 0],'markersize',(1.5*PPD)/2,'linewidth',1.5);
                plot(tarPosLand(1),tarPosLand(2),'.','color',[0.6 0 0],'markersize',30);
                
                for i = 1:size(mrs,1)
                    plot(xrsf(mrs(i,1):mrs(i,2),1),xrsf(mrs(i,1):mrs(i,2),2),'r-','linewidth',2);
                    % plot(sxrs(mrs(i,1):mrs(i,2),1),sxrs(mrs(i,1):mrs(i,2),2),'r-','linewidth',2);
                end
                if ~isempty(mrs)
                    plot([xrsf(mrs(:,1),1) xrsf(mrs(:,2),1)]',[xrsf(mrs(:,1),2) xrsf(mrs(:,2),2)]','r.');
                    % plot([sxrs(mrs(:,1),1) sxrs(mrs(:,2),1)]',[sxrs(mrs(:,1),2) sxrs(mrs(:,2),2)]','r.');
                end
                
                axes(ax(2));
                plot(timers-tCueOnEDF,vrsf(:,1),'k-','color',[0.6 0.6 0.6]);
                hold on;
                % plot(timers-tCueOnEDF,vrs(:,1),'g-');
                for i = 1:size(mrs,1)
                    %plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrs(mrs(i,1):mrs(i,2),1)','r-','linewidth',2);
                    plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrsf(mrs(i,1):mrs(i,2),1)','r-','linewidth',2);
                end
                hold off;
                xlim(timers([1 end])-tCueOnEDF);
                ylim([-1000 1000]);
                
                axes(ax(3));
                plot(timers-tCueOnEDF,vrsf(:,2),'k-','color',[0.6 0.6 0.6]);
                hold on;
                % plot(timers-tCueOnEDF,vrs(:,2),'g-');
                for i = 1:size(mrs,1)
                    % plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrs(mrs(i,1):mrs(i,2),2)','r-','linewidth',2);
                    plot(timers(mrs(i,1):mrs(i,2))-tCueOnEDF,vrsf(mrs(i,1):mrs(i,2),2)','r-','linewidth',2);
                end
                hold off;
                xlim(timers([1 end])-tCueOnEDF);
                ylim([-1000 1000]);
                
                % post saccadic trace
                axes(ax(4));
                
                if pursuitOK
                
                    polar(0, max(rho_p));
                    
                    hold on;
                    
                    tp = polar([tarLandPosDirAngle tarLandPosDirAngle], [0 max(rho_p)+0.1],'--');
                    set(tp,'Color',[0.5 0.5 0.5],'LineWidth',2);
                    
                    tp = polar([tarLandDirAngle tarLandDirAngle], [0 max(rho_p)+0.1],'--');
                    set(tp,'Color','black','LineWidth',2);
                    
                    tp = polar(theta_p, rho_p);
                    set(tp,'Color','red','LineWidth',2);
                    
                    hold off;
                
                else
                    plot([0,1],[0,1],'black')
                end
                
            end
            %%
            % loop over all necessary response saccades
            s = 0;
            for rs = 1:nRS
                fixRec = repmat(tarPos(rs  ,:),1,2)+[-tarRad -tarRad tarRad tarRad];
                tarRec = repmat(tarPos(rs+1,:),1,2)+10*[-tarRad -tarRad tarRad tarRad];
                
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
                xlim([-18 18]);
                ylim([-18 18]);
                text(-14,14,sprintf('%s',substr(vpcode, 3, 2)),'Fontsize',14,'hor','left','ver','bottom');
                text(-14,-14,sprintf('Trial %i (b. %i)',trial,block),'Fontsize',14,'hor','left','ver','bottom');
                
                if dir == 1 % 1 = CLOCKWISE; -1 = COUNTERCLOCKWISE
                    text(-14,-16,'target motion: clockwise','Fontsize',14,'hor','left','ver','bottom');
                else
                    text(-14,-16,'target motion: counterclockwise','Fontsize',14,'hor','left','ver','bottom');
                end
                
                condStr = 'Saccade';
                
                if sum(sbrs)==0 && sum(mbrs)==0 && sum(resp)==nRS % && samp==0
                    text(0,-14,sprintf('%s: good',condStr),'Fontsize',14,'hor','center','ver','bottom');
                else
                    text(0,-14,sprintf('%s: bad',condStr),'Fontsize',14,'hor','center','ver','bottom');
                end
                
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
                rea = [repmat(tab(t,:),nRS,1) reaSac pursuitRes td.dT curvatureData];
                
                % pursuit res.
                % %.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n
                % [tarDirAngle, velDev, pursuitVel, pursuitVelAngle, posDiffDev, pursuitPosDiff, pursuitPosDiffAngle, tarLandDirAngle, tarLand, tarLandPosDirAngle, pursuitAcc]
                
                % tab format
                % %i\t%i\t%.4f\t%i\t%i\t%i\t%.4f\t%i\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n
                
                % this is what you need in addition to the tab format:
                % %i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n
                              %i\t%i\t%.4f\t%i\t%i\t%i\t%.4f\t%i\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n
                fprintf(frea,'%i\t%i\t%.4f\t%i\t%i\t%i\t%.4f\t%i\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',rea);
                                                                                                                                                                                                                   
                if welche == 1
                    rea = [vp*ones(size(rea,1),1), rea];
                                     %i\t%i\t%i\t%.4f\t%i\t%i\t%i\t%.4f\t%i\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n
                    fprintf(freaAll,'%i\t%i\t%i\t%.4f\t%i\t%i\t%i\t%.4f\t%i\t%.4f\t%.4f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n',rea);
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
end

fprintf(1,'\n\nTotal came %i good trials out of %i total trials!',ntGoodAll,ntGoodAll+ntBadAll);
