%% PURSUIT ANALYSIS
%
% esp.drifting-moving gabor


if size(mrs,1)>=1
    
    % pre-define pursuit variables
    tarDirAngle = NaN;          velDev = NaN;           pursuitVel = NaN;
    pursuitVelAngle = NaN;      posDiffDev = NaN;       pursuitPosDiff = NaN;
    pursuitPosDiffAngle = NaN;  tarLandDirAngle = NaN;  tarLand = NaN;
    tarLandPosDirAngle = NaN;   pursuitAcc = NaN;
    
    nsamples = size(xrsf,1);
    
    
    % Here I calculate target parameters that will be compared
    % to pursuit parameter. These will also be used for
    % plotting (if plotting is enabled)
    
    % calculate target position at the landing time of the
    % saccade landing
    travelTime = timers(mrs(1,2))-tCueOnEDF; % dopo la prima saccade
    tarLand = tar - dir * deg2rad(angVel) *(travelTime/1000);
    
    % target [x, y] position at landing
    tarPosLand = [ecc*cos(tarLand), ecc*sin(tarLand)];
    
    % target direction angle
    tarDirAngle = tar - dir*(pi/2);
    
    % target direction angle at landing time
    tarLandDirAngle = tarLand - dir*(pi/2);
    
    % target direction at landing position
    sacAngle = cart2pol(xrs(mrs(1,2),1), xrs(mrs(1,2),2));
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
            
            % Here I calculate pursuit parameters (direction,
            % peak velocity, accelleration, ecc.). Direction is
            % calculated by fitting a line through the pursuit
            % gaze points
            
            rel_p_x = posPursuit(:,1) - posPursuit(1,1);
            rel_p_y = posPursuit(:,2) - posPursuit(1,2);
            
            % transform in polar coordinates for plotting
            [theta_p, rho_p] = cart2pol(rel_p_x, rel_p_y);
            
            % direction angle
            theta_linearFit = compDirection(rel_p_x, rel_p_y);  % default: average between the first and third part
            
            % horizontal and vertical components (position)
            pos_p = [(posPursuit(end,1) - posPursuit(1,1)), (posPursuit(end,2) - posPursuit(1,2))];
            
            % peak horizontal and vertical velocities
            vel_p = [absMax(velPursuit(:,1)), absMax(velPursuit(:,2))];
            
            % peak accellerations
            acc_p = [absMax(accPursuit(:,1)), absMax(accPursuit(:,2))];
            
            
            % calcolo deviazioni dalla direzione del target
            posDiffDev = angDiff(tarLandDirAngle, pursuitPosDiffAngle);
            velDev = angDiff(tarLandDirAngle, pursuitVelAngle);
            
            tarTanVel = angVel/ecc;
            
            % target trajectory within pursuit interval
            xRange = [1 length(actPursuitInterval)];
            
            % calculate expected target path (during pursuit interval) (i.e., TIME)
            timeInd = xRange(1):xRange(2);
            tarPathAngles = tarLand - dir .* deg2rad(angVel) .*(timeInd./1000); % tarLand is at landing time
            tarPath = [ecc.*cos(tarPathAngles); ecc.*sin(tarPathAngles)]';      % so this is target path in pursuit interval
            
            tarVel = 1000*[diff(tarPath(:,1)) diff(tarPath(:,2))];
            tarAcc = 1000*[diff(tarVel(:,1)) diff(tarVel(:,2))];
            
            % target expected direction in pursuit arc/chord (i.e., POSITION)
            tarPathAngles_pos = sacAngle - dir .* deg2rad(angVel) .*(timeInd./1000);
            tarPath_pos = [ecc.*cos(tarPathAngles_pos); ecc.*sin(tarPathAngles_pos)]';
            [tarPursuitDirAngle_pos] = cart2pol(tarPath_pos(end,1)-tarPath_pos(1,1), tarPath_pos(end,2)-tarPath_pos(1,2));
            
            % target direction angle in pursuit interval (it is
            % the angle of the chord, not the tnagential direction) TIME
            [tarPursuitDirAngle] = cart2pol(tarPath(end,1)-tarPath(1,1), tarPath(end,2)-tarPath(1,2));
            
            % add target direction at SACCADE ONSET
            travelTime2 = timers(mrs(1,1))-tCueOnEDF; % dopo la prima saccade
            tarSacOnset = tar - dir * deg2rad(angVel) *(travelTime2/1000);
            tarPathAngles_sacOnset = tarSacOnset - dir .* deg2rad(angVel) .*(timeInd./1000);
            tarPath_sacOnset = [ecc.*cos(tarPathAngles_sacOnset); ecc.*sin(tarPathAngles_sacOnset)]';
            [tarPursuitDirAngle_sacOnset] = cart2pol(tarPath_sacOnset(end,1)-tarPath_sacOnset(1,1), tarPath_sacOnset(end,2)-tarPath_sacOnset(1,2));
            
            % output of pusuit analysis
            pursuitRes = [tarPursuitDirAngle, tarLand, tarPursuitDirAngle_sacOnset, tarPursuitDirAngle_pos, pos_p, theta_linearFit, vel_p, acc_p]; %11
            
            % save pursuit trace
            ptcount = ptcount +1;
            pursuitTraces.x = alignTrial(pursuitTraces.x, posPursuit(:,1));
            pursuitTraces.y = alignTrial(pursuitTraces.y, posPursuit(:,2));
            pursuitTraces.trialInfo(ptcount) = td;
            pursuitTraces.sacAngle(ptcount) = sacAngle;
            
            if plotData
                
                % plot pursuit trace
                figure(h2);
                
                subplot(2,1,1)
                plot(tarPath,'--','LineWidth',1); hold on
                pl = plot(posPursuit,'LineWidth',2); hold off
                xlim(xRange);
                ylabel('position [deg]')
                legend(pl,'horizontal','vertical','Location','NorthEast');
                
                subplot(2,1,2)
                plot(tarVel,'--','LineWidth',1); hold on
                plot(velPursuit,'LineWidth',2); hold off
                xlim(xRange);
                ylabel('velocity [deg/sec]')
                % legend(pl,'horizontal','vertical','Location','NorthEastOutside');
                
            end
            
        else
            pursuitRes = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        end
        
    else
        pursuitRes = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];
        
    end
end