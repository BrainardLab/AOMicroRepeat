%% CombineTrials
%
% Combine all of the psychophysical data files into one big .mat
% file.  We do this in part to simplify the programs, and in part because
% the fitting program is crashing randomly and perhaps separating these
% two parts will help track it down.
%
% Also identify catch trials in QUEST runs, which are somewhat obscurely
% indicated.

%% Initialize
close all; clear all;

%% Get path to data
%
% Path to data tree on this machine
%
% This is set up by TbTb local hook file, but
% you can also use
%  setpath('AOMicroRepeat','dataDir',theDataDir)
% to do this, where theDataDir is the path to the
% files.
%
% In the code example I got, this was
%  path = 'W:\Data\11125\20231103\AO_Psychophysics\MOCS\8\group2';
% but that does not match the actuall data tree I received for 11002.
dataDir = getpref('AOMicroRepeat','dataDir');

%% Also analysis output dir, same idea
analysisDir = getpref('AOMicroRepeat','analysisDir');

%% Some parameters
log0Value = -3.5;

%% Define what data we are analyzing here
%
% Generally speaking, this routine loops over everything and does its thing
%
% Define subjects
theParticipants = {'11002' '11108' '11118' '11119' '11125'};

% Define sessions (1 or 2)
theSessions = [1 2];

% Define session splits
theSplits = {'All', 'FirstHalf', 'SecondHalf'};

% Define sizes (8 or 43)
theDiameters = [8 43];

% Define methods ('MOCS' or 'QUEST')
%
% When we create COMBINED data below, we count
% on this being as it is.  Don't change wihtout care.
theMethods = {'MOCS' 'QUEST', 'COMBINED'};

%% Get the AOM lookup table info.
%
% We use this below to masssage the nominal data, given
% the lookup table
%
% As far as I can guess
%   Column 1: nominal linear intensities
%   Column 2: LUT linear intensities
%   Column 3: nominal log10 intensities
%   Column 4: looks like a non-linear mapping from 10 to 8 bits
AOM = load('green_AOM_LUT_processing');

%% Loop over everything
tableRow = 1;
for pp = 1:length(theParticipants)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            for hh = 1:length(theSplits)
                checkSessionDate = [];
                MOCSFileTimes = [];
                QUESTFileTimes = [];
                for mm = 1:length(theMethods)

                    % Store info for what we are analyzing in this run
                    theMethod{tableRow,1} = theMethods{mm};
                    theSubject{tableRow,1} = theParticipants{pp};
                    theDiameter(tableRow,1) = theDiameters(dd);
                    theSession(tableRow,1) = theSessions(ss);
                    theSplit{tableRow,1} = theSplits{hh};

                    % Handle MOCS and QUEST from the data
                    if (strcmp(theMethods{mm},'MOCS') | strcmp(theMethods{mm},'QUEST'))

                        % Form path the the directory with the data sitting in it
                        pathToData = fullfile(dataDir,theSubject{tableRow},['Session' num2str(theSession(tableRow))],['Size' num2str(theDiameter(tableRow))],theMethod{tableRow});
                        pathToAnalysis = fullfile(analysisDir,theSubject{tableRow},['Session' num2str(theSession(tableRow))],['Size' num2str(theDiameter(tableRow))],theMethod{tableRow},theSplit{tableRow});
                        if (~exist(pathToAnalysis,'dir'))
                            mkdir(pathToAnalysis);
                        end

                        % Get list of data files
                        dirOffset = 1;
                        trial_videos = dir(fullfile(pathToData,'*.mat'));
                        num_trial_videos = length(trial_videos(dirOffset:end));
                        if (num_trial_videos ~= 4)
                            error('Expect 4 data files');
                        end

                        % Read and concatenate the data
                        all_trials{pp,dd,ss,hh,mm} = {};
                        all_trials_unpacked{pp,dd,ss,hh,mm} = [];                
                        fprintf('\tReadng videos\n');
                        for i = 1:num_trial_videos
                            % Get check date from filename.  This should be the
                            % same for all the MOCS and QUEST files from the
                            % same subject session.  The method is the inner
                            % loop variable, so we can do the check here.
                            if (isempty(checkSessionDate))
                                checkSessionDate = trial_videos(i+dirOffset-1).name(7:16);
                            else
                                if (~strcmp(checkSessionDate,trial_videos(i+dirOffset-1).name(7:16)))
                                    error('Not all data files from same session are on the same date');
                                end
                            end

                            % Check that QUEST file times are later than last
                            % MOCS time.  This code assumes that it runs to
                            % completion on the same calendar day on which it
                            % was started, which will almost always be true.
                            %
                            % Count underscores to find start of the time
                            % string.
                            nUnderscore = 0;
                            for cc = 1:length(trial_videos(i+dirOffset-1).name)
                                if (trial_videos(i+dirOffset-1).name(cc) == '_')
                                    nUnderscore = nUnderscore + 1;
                                end
                                if (nUnderscore == 4)
                                    timeStartIndex = cc + 1;
                                    break;
                                end
                            end
                            fileTimeStr = trial_videos(i+dirOffset-1).name(timeStartIndex:timeStartIndex+4);

                            % Handle way single digit time numbers get written
                            % in the filename.  Brute force this. Ugh.
                            if (fileTimeStr(2) == '_')
                                fileTimeStr = ['0' fileTimeStr];
                            end
                            if (fileTimeStr(5) == '_')
                                fileTimeStr(5) = fileTimeStr(4);
                                fileTimeStr(4) = '0';
                            end
                            fileTimeStr = fileTimeStr(1:5);

                            % Convert to format MATLAB can operate on
                            fileTime = datetime(fileTimeStr,'InputFormat','HH_mm');
                            if (strcmp(theMethod(tableRow),'MOCS'))
                                MOCSFileTimes = [MOCSFileTimes ; fileTime];
                            else
                                QUESTFileTimes = [QUESTFileTimes ; fileTime];
                            end

                            % Read the data file
                            % fprintf('Reading file %s\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name));
                            loadedData{pp,dd,ss,hh,mm,i} = load(fullfile(pathToData,trial_videos(i+dirOffset-1).name));

                            % Grab the trial data we need
                            if (isfield(loadedData{pp,dd,ss,hh,mm,i},'trial_vector'))
                                if (~strcmp(theMethod{tableRow},'MOCS'))
                                    fprintf('File %s\n\tUnder QUEST, appears to be MOCS\n\tTrials: %d\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name),length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector));
                                end
                                if (length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector) ~= 90 & length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector) ~= 100)
                                    fprintf('File %s: Wrong number of trials in MOCS data file\n');
                                end
                                all_trials{pp,dd,ss,hh,mm}{i,1} = loadedData{pp,dd,ss,hh,mm,i}.trial_vector;
                                all_trials{pp,dd,ss,hh,mm}{i,2} = loadedData{pp,dd,ss,hh,mm,i}.response_vector;
                            elseif (isfield(loadedData{pp,dd,ss,hh,mm,i},'trial_matrix'))
                                if (~strcmp(theMethod{tableRow},'QUEST'))
                                    fprintf('File %s\n\tUnder MOCS, appears to be QUEST\n\tTrials: %d\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name),length(loadedData{pp,dd,ss,hh,mm,i}.trial_matrix));
                                end
                                if (length(loadedData{pp,dd,ss,hh,mm,i}.trial_matrix) ~= 44)
                                    fprintf('File %s: Wrong number of trials in QUEST data file\n');
                                end
                                all_trials{pp,dd,ss,hh,mm}{i,1} = loadedData{pp,dd,ss,hh,mm,i}.trial_matrix;
                                questCatchTrialIndex = find(isnan(loadedData{pp,dd,ss,hh,mm,i}.theThreshold));
                                if (length(questCatchTrialIndex) ~= 4)
                                    fprintf('Wrong number of catch trials in QUEST data file');
                                end
                                all_trials{pp,dd,ss,hh,mm}{i,1}(questCatchTrialIndex) = log0Value;
                                all_trials{pp,dd,ss,hh,mm}{i,2} = loadedData{pp,dd,ss,hh,mm,i}.response_matrix;
                            else
                                error('Do not understand data format');
                            end

                            % Get/check number of trials
                            if (i == 1)
                                nTrialsPerSession = length(all_trials{pp,dd,ss,hh,mm}{i,1});
                            else
                                if (length(all_trials{pp,dd,ss,hh,mm}{i,1}) ~= nTrialsPerSession)
                                    error('Trial mismatch across runs');
                                end
                            end

                            % Concatenate
                            all_trials_unpacked{pp,dd,ss,hh,mm} = [all_trials_unpacked{pp,dd,ss,hh,mm} ; [all_trials{pp,dd,ss,hh,mm}{i,1} all_trials{pp,dd,ss,hh,mm}{i,2}] ];
                        end

                        if (strcmp(theMethod{tableRow},'MOCS') & size(MOCSFileTimes,1) ~= 4)
                            error('Wrong number of MOCS files somewhere in data tree');
                        end
                        if (strcmp(theMethod{tableRow},'QUEST')  & size(QUESTFileTimes,1) ~= 4)
                            error('Wrong number of QUEST files somewhere in data tree');
                        end
                        if (strcmp(theMethod{tableRow},'QUEST'))
                            QUESTFileTimesSorted = sort(QUESTFileTimes);
                            MOCSFileTimesSorted = sort(MOCSFileTimes);
                            for ff = 1:size(MOCSFileTimes,1)
                                if (QUESTFileTimesSorted(ff) < MOCSFileTimesSorted(ff))
                                    fprintf('File time order error\n');
                                end
                            end
                        end
                    else
                        % Concatenate MOCS and QUEST data into COMBINED
                        all_trials_unpacked{pp,dd,ss,hh,mm} = [all_trials_unpacked{pp,dd,ss,hh,1} ; all_trials_unpacked{pp,dd,ss,hh,2}];
                    end

                    % Bump table row
                    tableRow = tableRow + 1;
                end
            end
        end
    end
end

% Save out one nice big combined file
save(fullfile(analysisDir,'combinedData.mat'),'all_trials','all_trials_unpacked','log0Value','theParticipants','theDiameters','theSessions','theSplits','theMethods','AOM','-v7.3');


