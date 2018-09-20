% xmsg2tab.m
%
% Experiment: Drifting gabor - moving envelope. saccade task
%
% Creates tab-File containing information specified
% for a certain trial.
%
% tab-Files have the following columns
% (predefined by defNaNVariables.m):

clear;

addpath('functions/');

msgpath = '../raw/';
tabpath = '../tab/';

subfid = fopen('subjects.tmp','r');
cnt = 1;
while cnt ~= 0
    [vpcode, cnt] = fscanf(subfid,'%s',1);
    if cnt ~= 0
        msgstr = sprintf('%s%s.msg',msgpath,vpcode);
        msgfid = fopen(msgstr,'r');

        fprintf(1,'\nprocessing ... %s.msg',vpcode);
        stillTheSameSubject = 1;
        tab = [];
        while stillTheSameSubject
            % predefine critical variables
            defNaNVariables;

            stillTheSameTrial = 1;
            while stillTheSameTrial

                line = fgetl(msgfid);
                if ~ischar(line)    % end of file
                    stillTheSameSubject = 0;
                    break;
                end

                if ~isempty(line) && stillTheSameSubject   % skip empty lines
                    la = strread(line,'%s');    % array of strings in line

                    if length(la) >= 3
                        switch char(la(3))
                            case 'TRIAL_START'
                                trial = str2double(char(la(4)));
                            case 'EVENT_FixationDot'
                                tedfFix = str2double(char(la(2)));
                            case 'EVENT_TargetOnset'
                                tedfTOn = str2double(char(la(2)));
                            case 'EVENT_Flash'
                                tedfFixOff = str2double(char(la(2)));
                            case 'EVENT_Saccade1Started'
                                tedfSac = str2double(char(la(2)));
                            case 'TRIAL_ENDE'
                                trial2 = str2double(char(la(4)));
                            case 'TrialData'
                                
                                % 4 5       6       7           8           9           10      11  12      13      14      15      16   17     18
                                % b trial   alpha   env_speed   drift_speed trajLength  initPos ecc fixOff  tFix    tBeg    tFixOff tSac tEnd   sacRT
                                
                                trial3 = str2double(char(la(5)));
                                block = str2double(char(la(4)));
                                alpha  = str2double(char(la(6)));
                                env_speed = str2double(char(la(7)));
                                drift_speed = str2double(char(la(8)));
                                trajLength = str2double(char(la(9)));
                                initPos = str2double(char(la(10)));
                                ecc = str2double(char(la(11)));
                                fixOff = str2double(char(la(12)));
                                
                                % presentation times are relative to the first frame of motion ([tFix tBeg tFixOff tSac tEnd]-tBeg)
                                tFix = str2double(char(la(13)));
                                tBeg = str2double(char(la(14)));
                                tFixOff = str2double(char(la(15)));
                                tSac = str2double(char(la(16)));
                                tEnd = str2double(char(la(17)));
                                
                                sacRT = str2double(char(la(18)));
                                other1 = NaN;
                                other2 = NaN;
                                other3 = NaN;
                                
                                stillTheSameTrial = 0;
         
                        end
                    end
                end
            end
            
            % check if trial ok and all messages available
            if trial==trial3 && sum(isnan([trial trial3 tedfFix tedfTOn tedfFixOff tedfSac]))==0
                everythingAvailable = 1;
            else
                everythingAvailable = 0;
            end

            if everythingAvailable
%                 tedfTOn = tedfTOn - tedfFix;
%                 tedfSac = tedfSac - tedfFix;
%                 tedfTOf = tedfTOf - tedfFix;
%                 tedfClr = tedfClr - tedfFix;
%                 tedfChange = tedfChange - tedfFix;
                
                % create uniform data matrix containing any potential
                % information concerning a trial
                tab = [tab; block trial3 alpha env_speed drift_speed trajLength initPos ecc tFixOff tFix tBeg tFixOff tSac tEnd sacRT other1 other2 other3 tedfFix tedfTOn tedfFixOff tedfSac];
                
            elseif trial~=trial2
                fprintf(1,'\nMissing Message between TRIALID %i and trialData %i (ignore if last trial)',trial,trial2);
            end
        end
        fclose(msgfid);

        outname = sprintf('%s%s.tab',tabpath,vpcode);
        fout = fopen(outname,'w');
        
        % block trial3 alpha env_speed drift_speed trajLength initPos ecc fixOff tFix tBeg tFixOff tSac tEnd sacRT other1 other2 other3 tedfFix tedfTOn tedfFixOff tedfSac
        % %i\t%i\t.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\n
        
        fprintf(fout,'%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%.4f\t%.4f\t%.4f\t%i\t%i\t%i\t%i\n',tab');

        
        fclose(fout);
    end
end
fclose(subfid);
fprintf(1,'\n\nOK!!\n');
