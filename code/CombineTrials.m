%% CombineTrials
%
% Combine all of the psychophysical data files into one big .mat
% file.  We do this to simplify the other programs, since the logic
% here involves lots of checks that everything is as it should be.
%
% Also identify catch trials in QUEST runs, which are somewhat obscurely
% indicated, and fix them up so later programs don't have to worry about
% this.
%
% Although this might look like it is splitting the data into two halves, that
% actually happens elsewhere, and here all the data are stored under all the splits.
% This is a vestige from the time this code was part of the full data analysis loop, and
% if we were better people we'd skip this step here.  Maybe we will become better at some
% point.
%
% The output of this program is a mat file 'CombinedData.mat' that gets stored in the
% analysis output tree, and that contains one big cell array with all of the data.
%
% See also: FitTrials, FigureN

%% Initialize
close all; clear all;

%% Get path to data
%
% Path to data tree on this machine
%
% This is set up by TbTb local hook file, but you can also do it
% by setting up the relevant preferences one time on your machine,
% for example like this:
  setpref('AOMicroRepeat','dataDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\New_analysis_20250912\dataDir');
  setpref('AOMicroRepeat','analysisDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\figure-rerun');
%
% Do not hard code directories into this program.
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
%
% We just do 'All' in the loop, and then hand code in the 'Group1' and 'Group2' splits,
% because the way the code is structured makes that easier to write clearly. Be sure not
% to change this order because it is used in FitTrials.m
% theSplits = {'All','Group1','Group2'};
theSplits = {'All'};

% Define sizes (8 or 43)
theDiameters = [8 43];

% When we create COMBINED data below, we count
% on this being as it is.  Don't change wihtout care.
theMethods = {'MOCS' 'QUEST', 'COMBINED'};

%% Get the AOM lookup table info.
%
% This is loaded and saved here with the data, so that the programs that
% read the output of this program don't have to go and find it.
%
% As far as I can guess
%   Column 1: nominal linear intensities
%   Column 2: LUT linear intensities
%   Column 3: nominal log10 intensities
%   Column 4: looks like a non-linear mapping from 10 to 8 bits
AOM = load('green_AOM_LUT_processing');

%% Lock the rng
rng(0);

%% Loop over everything
tableRow = 1;
for pp = 1:length(theParticipants)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            for hh = 1:length(theSplits)
                checkSessionDate = [];
                checkSessionTime = [];
                MOCSFileTimes = [];
                QUESTFileTimes = [];
                runOrderSplitThisSessionForBothMethods = Shuffle([1 1 2 2]);
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

                            % Grab the trial data we need.  Also do some
                            % MOCS and QUEST specific reality checks
                            if (isfield(loadedData{pp,dd,ss,hh,mm,i},'trial_vector'))
                                % It's a MOCS data file by what's in it
                                if (~strcmp(theMethod{tableRow},'MOCS'))
                                    fprintf('File %s\n\tUnder QUEST, appears to be MOCS\n\tTrials: %d\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name),length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector));
                                end
                                if (length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector) ~= 90 & length(loadedData{pp,dd,ss,hh,mm,i}.trial_vector) ~= 100)
                                    fprintf('File %s: Wrong number of trials in MOCS data file\n');
                                end
                                all_trials{pp,dd,ss,hh,mm}{i,1} = loadedData{pp,dd,ss,hh,mm,i}.trial_vector;
                                all_trials{pp,dd,ss,hh,mm}{i,2} = loadedData{pp,dd,ss,hh,mm,i}.response_vector;

                                % Reality check on MOCS stimulus levels.
                                % This should be the same across runs
                                % within session, I think, but are not
                                % always. This doesn't error out, but does
                                % report the unexpected cases.
                                for j = i-1:-1:1
                                    if (any(unique(loadedData{pp,dd,ss,hh,mm,j}.trial_vector) ~= unique(loadedData{pp,dd,ss,hh,mm,i}.trial_vector)))
                                        fprintf('\t\tMOCS mismatch in trial levels across runs\n')
                                        fprintf('\t\t%s, session %d, size %d, run loaded #%d vs run loaded #%d\n',theParticipants{pp},theSessions(ss), theDiameters(dd),i,j);
                                        fprintf('\t\tFile loaded #%d: %s\n',i,fullfile(trial_videos(i+dirOffset-1).name));
                                        fprintf('\t\tFile loaded #%d: %s\n',j,fullfile(trial_videos(j+dirOffset-1).name));
                                    end
                                end

                            elseif (isfield(loadedData{pp,dd,ss,hh,mm,i},'trial_matrix'))
                                % It's a QUEST data file by what's in it
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

                            % Get/check number of trials and make sure they
                            % are the same for each run within
                            % method/session.
                            if (i == 1)
                                nTrialsPerSession = length(all_trials{pp,dd,ss,hh,mm}{i,1});
                            else
                                if (length(all_trials{pp,dd,ss,hh,mm}{i,1}) ~= nTrialsPerSession)
                                    error('Trial mismatch across runs');
                                end
                            end

                            % Concatenate across runs into one long pair of vectors.
                            all_trials_unpacked{pp,dd,ss,hh,mm} = [all_trials_unpacked{pp,dd,ss,hh,mm} ; [all_trials{pp,dd,ss,hh,mm}{i,1} all_trials{pp,dd,ss,hh,mm}{i,2}] ];

                            % Do the splits for this session.  Half the runs go into one
                            % split, half into the other.  Here's what the
                            % pre-registration says about how we'll do: What the
                            % pre-registration says we'll do is: "To measure intra-session
                            % variation in sensitivity threshold, the data from each
                            % session will be randomly split in half by experimental run,
                            % psychometric functions will be fit to each half of the data
                            % and thresholds will be estimated and compared for each
                            % condition." This is not completely clear about whether we
                            % should keep MOCS and QUEST paired runs together, but we
                            % decided that is most sensible.

                            % Make sure hh == 1 here, because logic now depends on it
                            % never being bigger than 1.
                            if (hh ~= 1)
                                error('Someone added splits, but the logic does not support that anymore here');
                            end
                            if (runOrderSplitThisSessionForBothMethods(i) == 1)
                                all_trials_unpacked{pp,dd,ss,2,mm} = [all_trials_unpacked{pp,dd,ss,hh,mm} ; [all_trials{pp,dd,ss,hh,mm}{i,1} all_trials{pp,dd,ss,hh,mm}{i,2}] ];
                            else
                                all_trials_unpacked{pp,dd,ss,3,mm} = [all_trials_unpacked{pp,dd,ss,hh,mm} ; [all_trials{pp,dd,ss,hh,mm}{i,1} all_trials{pp,dd,ss,hh,mm}{i,2}] ];
                            end
                        
                        end

                        % Checks
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
                        % Concatenate MOCS and QUEST data into COMBINED.  If we're here,
                        % mm indeices COMBINED and we know that indices 1 and 2 index MOCS
                        % and QUEST.
                        if (~strcmp(theMethods{mm},'COMBINED') | ~strcmp(theMethods{1},'MOCS') | ~strcmp(theMethods{2},'QUEST') )
                            error('Code counts on specific order of methods in theMethods, and someone has changed that');
                        end
                        if (~strcmp(theSplits{hh},'All'))
                            error('In this program, the only split explicitly specified should be ''All''');
                        end

                        % Check that hh never exceeds 1
                        if (hh ~= 1)
                            error('Someone added splits, but the logic does not support that anymore here');
                        end

                        % Combine across methods for 'All' and the two groups.
                        all_trials_unpacked{pp,dd,ss,hh,mm} = [all_trials_unpacked{pp,dd,ss,hh,1} ; all_trials_unpacked{pp,dd,ss,hh,2}];
                        all_trials_unpacked{pp,dd,ss,2,mm} = [all_trials_unpacked{pp,dd,ss,2,1} ; all_trials_unpacked{pp,dd,ss,hh,2}];
                        all_trials_unpacked{pp,dd,ss,3,mm} = [all_trials_unpacked{pp,dd,ss,3,1} ; all_trials_unpacked{pp,dd,ss,hh,2}];

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


