
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
theMethods = {'MOCS' 'QUEST'};

%% Some parameters
log0TrialThreshold = -3;
log0Value = -3.5;
thresholdCriterion = 0.78;
nBootstraps = 500;
convertToDb = true;
justCheckFiles = false;

%% Freeze rng seed for repeatability
rng(101);

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

                    % Say hello
                    if (~justCheckFiles)
                        fprintf('Analyzing condition %d\n',tableRow);
                    end

                    % Store info for what we are analyzing in this run
                    theMethod{tableRow,1} = theMethods{mm};
                    theSubject{tableRow,1} = theParticipants{pp};
                    theDiameter(tableRow,1) = theDiameters(dd);
                    theSession(tableRow,1) = theSessions(ss);
                    theSplit{tableRow,1} = theSplits{hh};

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
                    all_trials = {};
                    all_trials_unpacked = [];
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
                        tempData = load(fullfile(pathToData,trial_videos(i+dirOffset-1).name));

                        % Grab the trial data we need
                        if (isfield(tempData,'trial_vector'))
                            if (~strcmp(theMethod{tableRow},'MOCS'))
                                fprintf('File %s\n\tUnder QUEST, appears to be MOCS\n\tTrials: %d\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name),length(tempData.trial_vector));
                                if (~justCheckFiles)
                                    error('Inconsistency in our understanding of method');
                                end

                            end
                            if (length(tempData.trial_vector) ~= 90 & length(tempData.trial_vector) ~= 100)
                                fprintf('File %s: Wrong number of trials in MOCS data file\n');
                                if (~justCheckFiles)
                                    error('Wrong number of trials in MOCS data file');
                                end
                            end
                            all_trials{i,1} = tempData.trial_vector;
                            all_trials{i,2} = tempData.response_vector;
                        elseif (isfield(tempData,'trial_matrix'))
                            if (~strcmp(theMethod{tableRow},'QUEST'))
                                fprintf('File %s\n\tUnder MOCS, appears to be QUEST\n\tTrials: %d\n',fullfile(pathToData,trial_videos(i+dirOffset-1).name),length(tempData.trial_matrix));
                                if (~justCheckFiles)
                                    error('Inconsistency in our understanding of method');
                                end
                            end
                            if (length(tempData.trial_matrix) ~= 44)
                                fprintf('File %s: Wrong number of trials in QUEST data file\n');
                                if (~justCheckFiles)
                                    error('Wrong number of trials in QUEST data file');
                                end
                            end
                            all_trials{i,1} = tempData.trial_matrix;
                            questCatchTrialIndex = find(isnan(tempData.theThreshold));
                            if (length(questCatchTrialIndex) ~= 4)
                                fprintf('Wrong number of catch trials in QUEST data file');
                                if (~justCheckFiles)
                                    error('Wrong number of catch trials in QUEST data file');
                                end
                            end
                            all_trials{i,1}(questCatchTrialIndex) = log0Value;
                            all_trials{i,2} = tempData.response_matrix;
                        else
                            error('Do not understand data format');
                        end

                        % Get/check number of trials
                        if (~justCheckFiles)
                            if (i == 1)
                                nTrialsPerSession = length(all_trials{i,1});
                            else
                                if (length(all_trials{i,1}) ~= nTrialsPerSession)
                                    error('Trial mismatch across runs');
                                end
                            end

                            % Concatenate
                            all_trials_unpacked = [all_trials_unpacked ; [all_trials{i,1} all_trials{i,2}]];
                        end
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
                                if (~justCheckFiles)
                                    error('File time order error');
                                end
                            end
                        end
                    end

                    % If just checking, break and don't do anything else
                    if (~justCheckFiles)
                        % Clear out variables from previous time through
                        % the loop.
                        clear log_intensity_nominal lin_intensity_nominal lin_intensity_rounded log_intensity_nominal log_intensity_rounded log_intensity_rounded_chk
                        clear lut_row_index dac_intensity_lut lin_intensity_lut
                        % Do the lookup table conversion
                        for i = 1:size(all_trials_unpacked,1)
                            % Get the nominal log intensity from the data
                            log_intensity_nominal(i) = all_trials_unpacked(i,1);

                            % Convert to nominal lin intensity
                            lin_intensity_nominal(i) = 10.^log_intensity_nominal(i);

                            % Log10 trial values less than 3 get rounded
                            % down to 0 intensity.  Deal with this.
                            if (log_intensity_nominal(i) < log0TrialThreshold)
                                log_intensity_nominal(i) = log0Value;
                                lin_intensity_nominal(i) = 0;
                            end

                            % Round linear intensity nominal.  Should match a row first column of lut
                            lin_intensity_rounded(i) = round(lin_intensity_nominal(i)*1000)/1000;

                            % Get the rounded log intensity, which we think is what
                            % is in the third column of the lut.  Except that if
                            % the linear intensity is 0, then the entry of the
                            % third column is zero and not -Inf.  Handle
                            % this for the check, and map the -Inf to the
                            % log value we use for 0 for further analysis.
                            log_intensity_rounded(i) = log10(lin_intensity_rounded(i));
                            log_intensity_rounded_chk(i) = log_intensity_rounded(i);
                            if (lin_intensity_rounded(i) == 0)
                                log_intensity_rounded(i) = log0Value;
                                log_intensity_rounded_chk(i) = 0;
                            end

                            % Compute the row index into the lut based on the linear intensity
                            lut_row_index(i) = lin_intensity_rounded(i)*1000+1;

                            % Get the corrected linear intensity by
                            % averaging the nominal linear intensities that
                            % correspond to the DAC value actually used.
                            % Since we don't know precisely what the AOM
                            % does in this case, since the LUT was created
                            % by interpolation of measurements that were
                            % made long ago, that seems as clever as
                            % anything else I can think of.  One could also
                            % use the min or the max of the ambiguous set
                            % of values, or any weighted average.
                            dac_intensity_lut(i) = AOM.green_AOM_lut(lut_row_index(i),4);
                            questCatchTrialIndex = find(AOM.green_AOM_lut(:,4) == dac_intensity_lut(i) );
                            lin_intensity_lut(i) = mean(AOM.green_AOM_lut(questCatchTrialIndex,1));

                            % Convert to log, taking log10(0) to be -3.5
                            log_intensity_lut(i) = log10(lin_intensity_lut(i));
                            if (lin_intensity_lut(i) == 0)
                                log_intensity_lut(i) = -3.5;
                            end

                            % Sanity checks
                            if (abs(lin_intensity_rounded(i)-AOM.green_AOM_lut(lut_row_index(i),1)) > 1e-5)
                                error('Do not understand lut table and/or its indexing');
                            end
                            if (abs(log_intensity_rounded_chk(i)-AOM.green_AOM_lut(lut_row_index(i),3)) > 1e-5)
                                error('Do not understand lut table and/or its indexing');
                            end

                            % Tag corrected log and linear intensities into the
                            % all_trials_unpacked variable that we save
                            all_trials_unpacked(i,3) = lin_intensity_lut(i);
                            all_trials_unpacked(i,4) = log_intensity_lut(i);
                        end

                        % Figure to make sure LUT conversion is basically
                        % the identity.  Uncomment if you want to look at
                        % this.
                        % figure; subplot(1,2,1);
                        % plot(lin_intensity_nominal,lin_intensity_lut,'ro');
                        % axis('square'); axis([0 1 0 1]);
                        % subplot(1,2,2);
                        % plot(lin_intensity_nominal,lin_intensity_lut,'ro');
                        % axis('square'); axis([0 0.1 0 0.1]);
                        % drawnow;

                        % Here is the format of all_trials_unpacked
                        %
                        % These values are not adjusted for power measurements.
                        % they are scaled with the instrument maximum having a
                        % linear intensity of 1.
                        %
                        %  Column 1 - nominal log10 intensity
                        %  Column 2 - 1 = seen, 0 = not seen
                        %  Column 3 - lut corrected linear intensity
                        %  Column 4 - lut corrected log intensity
                        %
                        % Split the data
                        nTrials = size(all_trials_unpacked,1);
                        shuffleIndex = Shuffle(1:nTrials);
                        switch (theSplit{tableRow})
                            case 'All'
                                dataIndex = 1:nTrials;
                            case 'FirstHalf'
                                dataIndex = shuffleIndex(1:round(nTrials/2));
                            case 'SecondHalf'
                                dataIndex = shuffleIndex(round(nTrials/2)+1:nTrials);
                            otherwise
                                error('Unknown split specified');
                        end
                        all_trials_unpacked = all_trials_unpacked(dataIndex,:);
                        nTrials = size(all_trials_unpacked,1);

                        % Extract the core data to fit
                        trial_log_intensities = all_trials_unpacked(:,4);
                        trial_responses = all_trials_unpacked(:,2);

                        % Convert to dB?
                        if (convertToDb)
                            trial_log_intensities = 10*trial_log_intensities;
                            log0Value = 10*log0Value;
                        end

                        % Fit the psychometric function.  The fitting routine makes
                        % a plot and we adjust the title here for things the
                        % fitting routine doesn't know about.
                        [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
                            log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(trial_log_intensities,trial_responses,log0Value,thresholdCriterion,true,convertToDb);
                        if (convertToDb)
                            title({ sprintf('Subject %s, %s, session %d, split %s',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow}) ; ...
                                sprintf('Diameter %d arcmin, threshold (dB) %0.2f',theDiameter(tableRow),corrected_log_threshold) ; ...
                                sprintf('Slope %0.1f, guess %0.2f, lapse %0.2f',psiParamsValues(2),psiParamsValues(3),psiParamsValues(4)) ; ''},'FontSize',20);
                        else
                            title({ sprintf('Subject %s, %s, session %d, split %s',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow}) ; ...
                                sprintf('Diameter %d arcmin, log threshold %0.2f',theDiameter(tableRow),corrected_log_threshold) ; ...
                                sprintf('Slope %0.1f, guess %0.2f, lapse %0.2f',psiParamsValues(2),psiParamsValues(3),psiParamsValues(4)) ; ''},'FontSize',20);
                        end
                        drawnow;
                        saveas(h,fullfile(pathToAnalysis,'psychometricFcn.tif'),'tif');

                        % Bootstrap the data
                        bootstrap_log_intensities = zeros(nTrials,nBootstraps);
                        bootstrap_responses = zeros(nTrials,nBootstraps);
                        switch(theMethod{tableRow})
                            case 'QUEST'
                                % For QUEST, we draw with replacement from all
                                % the trials except the design catch
                                % trials, which we draw from separately.
                                questNonCatchTrialIndex = setdiff(1:nTrials,questCatchTrialIndex);
                                for bb = 1:nBootstraps
                                    temp_boot_intensities = [];
                                    temp_boot_responses = [];
                                    
                                    % Catch trials
                                    uindex = questCatchTrialIndex;
                                    utemp_intensities = trial_log_intensities(uindex);
                                    utemp_responses = trial_responses(uindex);
                                    bootIndex = randi(length(uindex),length(uindex),1);
                                    temp_boot_intensities = [temp_boot_intensities ; utemp_intensities(bootIndex)];
                                    temp_boot_responses = [temp_boot_responses ; utemp_responses(bootIndex)];

                                    % Non catch trials
                                    uindex = questNonCatchTrialIndex;
                                    utemp_intensities = trial_log_intensities(uindex);
                                    utemp_responses = trial_responses(uindex);
                                    bootIndex = randi(length(uindex),length(uindex),1);
                                    temp_boot_intensities = [temp_boot_intensities ; utemp_intensities(bootIndex)];
                                    temp_boot_responses = [temp_boot_responses ; utemp_responses(bootIndex)];

                                    bootstrap_log_intensities(:,bb) = temp_boot_intensities;
                                    bootstrap_responses(:,bb) = temp_boot_responses;
                                end
                                clear questCatchTrialIndex questNonCatchTrialIndex

                            case 'MOCS'
                                % For MOCS, we respect the experimental design
                                % and bootstrap each stimulus intensity
                                % separately.
                                unique_log_intensities = unique(trial_log_intensities);
                                for bb = 1:nBootstraps
                                    temp_boot_intensities = [];
                                    temp_boot_responses = [];
                                    for uu = 1:length(unique_log_intensities)
                                        uindex = find(trial_log_intensities == unique_log_intensities(uu));
                                        utemp_intensities = trial_log_intensities(uindex);
                                        utemp_responses = trial_responses(uindex);

                                        bootIndex = randi(length(uindex),length(uindex),1);
                                        temp_boot_intensities = [temp_boot_intensities ; utemp_intensities(bootIndex)];
                                        temp_boot_responses = [temp_boot_responses ; utemp_responses(bootIndex)];
                                    end
                                    bootstrap_log_intensities(:,bb) = temp_boot_intensities;
                                    bootstrap_responses(:,bb) = temp_boot_responses;
                                end
                            otherwise
                                error('Unknown method specified');
                        end

                        % Bootstrap the fits and extract CI.  Add CI to plot a
                        % a thick blue horizontal bar.
                        for bb = 1:nBootstraps
                            if (rem(bb,100) == 0)
                                fprintf('\tBootstrap fit %d of %d\n',bb,nBootstraps);
                            end
                            [~,~,~,~,boot_corrected_threshold(bb),~,~] = FitPsychometricFunction(bootstrap_log_intensities(:,bb),bootstrap_responses(:,bb),log0Value,thresholdCriterion,false,convertToDb);
                        end
                        boot_quantiles = quantile(boot_corrected_threshold,[0.025 0.5 0.975]);
                        figure(h);
                        plot([boot_quantiles(1) boot_quantiles(3)],[thresholdCriterion thresholdCriterion],'b','LineWidth',4);
                        saveas(h,fullfile(pathToAnalysis,'psychometricFcnCI.tif'),'tif');

                        % Save what we learned
                        save(fullfile(pathToAnalysis,'fitOutput.mat'),'all_trials_unpacked','plot_log_intensities','plot_psychometric', ...
                            'corrected_psychometric','log_threshold','corrected_log_threshold','psiParamsValues','boot_corrected_threshold','boot_quantiles','-v7.3');

                        % Accumulate data
                        theTableLogThreshold(tableRow,1) = log_threshold;
                        theTableCorrectedLogThreshold(tableRow,1) = corrected_log_threshold;
                        theTableCorrectedLogThresholdBootMedian(tableRow,1)= boot_quantiles(2);
                        theTableCorrectedBootCILow(tableRow,1) = boot_quantiles(1);
                        theTableCorrectedBootCIHigh(tableRow,1) = boot_quantiles(3);
                        theTablePsiParamsValues(tableRow,:) = psiParamsValues;
                        theTableThresholdCriterion(tableRow,1) = thresholdCriterion;

                        % Clear
                        clear all_trials all_trials_unpacked

                        % Bump table row
                        tableRow = tableRow + 1;
                    end
                end
            end
        end
    end

    % Close figures this subject
    close all;

    % Make any summary figures for this subject here
    if (~justCheckFiles)
    end
end

% Write out full analysis data table
if (~justCheckFiles)
    if (convertToDb)
        tableVariableNames = {'Subject','Diameter','Session','Method','Split','Threshold (dB)','Corrected Threshold (dB)','CI Low','CI High', 'Alpha','Beta','Guess','Lapse','Criterion'};
        thresholdPlaces = 1;
    else
        tableVariableNames = {'Subject','Diameter','Session','Method','Split','Log10 Threshold','Corrected Log10 Threshold','CI Low','CI High', 'Alpha','Beta','Guess','Lapse','Criterion'};
        thresholdPlaces = 2;
    end
    dataTable = table(theSubject,theDiameter,theSession,theMethod,theSplit, ...
        round(theTableLogThreshold,thresholdPlaces),round(theTableCorrectedLogThreshold,thresholdPlaces),round(theTableCorrectedBootCILow,thresholdPlaces), round(theTableCorrectedBootCIHigh,thresholdPlaces), ...
        round(theTablePsiParamsValues(:,1),2),round(theTablePsiParamsValues(:,2),2), round(theTablePsiParamsValues(:,3),2), round(theTablePsiParamsValues(:,4),2), ...
        theTableThresholdCriterion, ...
        'VariableNames',tableVariableNames);
    writetable(dataTable,fullfile(analysisDir,'fitTable.xlsx'),'FileType','spreadsheet');
end


% FitPsychometricFunction
%
% Massage data, fit and plot psychometric function
function [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(log_intensity,responses,log0Value,thresholdCriterion,make_plot,convert_to_db)

% Sort data and get what we need
sorted_data = sortrows([log_intensity responses],1);
sorted_log_intensities = sorted_data(:,1);
sorted_log_intensities(isinf(sorted_log_intensities)) = log0Value;
sorted_responses = sorted_data(:,2);

% Find raw guess rate
catch_trial_index = find(sorted_log_intensities == log0Value);
guess_rate = sum(sorted_responses(catch_trial_index))/length(catch_trial_index);

% Get unique stimulus levels and data for those
[unique_log_intensities] = unique(sorted_log_intensities,'legacy');
numpresented = zeros(size(unique_log_intensities));
for i = 1:length(numpresented)
    num_presented(i) = length(find(sorted_log_intensities == unique_log_intensities(i)));
    num_seen(i) = sum(sorted_responses(find(sorted_log_intensities == unique_log_intensities(i))));
end

% Fit using a wrapper into Palemedes
[plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold, psiParamsValues] = ...
    fitPsychometric(unique_log_intensities,num_seen',num_presented',thresholdCriterion,convert_to_db);

% Plot
if (make_plot)
    if (convert_to_db)
        minIntensityLog=-35;
    else
        minIntensityLog=-3.5;
    end
    maxIntensityLog= 0;
    fontsize = 18; fwidth = 600; fheight = 600;
    h = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
    set(gca,'FontName','Helvetical','FontSize',14);
    if (convert_to_db)
        xlabel('Trial intensity (dB)','FontSize',fontsize);
    else
        xlabel('Log trial intensity (au)','FontSize',fontsize);
    end
    ylabel('Fraction seen','FontSize',fontsize);
    xlim([minIntensityLog maxIntensityLog]);
    ylim([0 1]);
    set(a0,'FontSize',fontsize);
    set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(a0,'Color','none');
    set(h,'Color',[1 1 1]);
    set(h,'PaperPositionMode','auto');
    set(h, 'renderer', 'painters');
    baseMarkerSize = 10;
    sizeNormalizer = 40;
    for tt = 1:length(unique_log_intensities)
        markerSize = 4+round(baseMarkerSize*num_presented(tt)/sizeNormalizer);
        plot(unique_log_intensities(tt), num_seen(tt)' ./ num_presented(tt)','ro','MarkerSize',markerSize,'MarkerFaceColor','r');
    end
    plot(plot_log_intensities,corrected_psychometric,'k-','LineWidth',3);
    plot(plot_log_intensities,plot_psychometric,'r:','LineWidth',2);
    plot([minIntensityLog corrected_log_threshold],[thresholdCriterion thresholdCriterion],'k:','LineWidth',2);
    plot([corrected_log_threshold corrected_log_threshold],[0 thresholdCriterion],'k:','LineWidth',2);
else
    h  = [];
end

end

% fitPsychometric
%
% Wrapper into Palamedes that fits the logistic with parameters reasonable
% for our data.
function [fit_log_intensity,fit_psychometric,corrected_fit_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues] = fitPsychometric(log_intensity,num_pos,out_of_num,threshold_criterion,convert_to_db)
%
% Usage:
%   [fit_log_intensity,fit_psychometric,corrected_fit_psychonmetric, ...
%        log_threshold,corrected_log_threshold,psiParamsValues] = fitPsychometric(log_intensity,num_pos,out_of_num,threshold_criterion)
%
% Inputs:
%   log_intensity : (vector) log10 stimulus levels tested
%   num_pos : (vector) for each stimulus level tested, the number of trials a positive response was given
%   out_of_num : (vector) for each stimulus level tested, the total number of trials
%   threshold_criterion : [scalar] criterion for threshold
%
% History:
%   07/05/21  amn  Wrote it.
%   06/28/25  dhb  Adopted for this project

% Calculate x-axis values to plot
%
% Increase the number of stimulus level values to plot for a smooth fit.
if (convert_to_db)
    fit_log_intensity = linspace(-35,0,1000);
else
    fit_log_intensity = linspace(-3.5,0,1000);
end

% Psychometric function form (alternative: PAL_Gumbel).
PF = @PAL_Logistic;

% 'paramsFree' is a boolean vector that determines what parameters get
% searched over (1: free parameter, 0: fixed parameter).
paramsFree = [1 1 1 1];

% Range info
guessUpper = 0.10;
lapseUpper = 0.05;

% Set up starting points:
if (convert_to_db)
    searchGrid.alpha = -35:.1:0;
    searchGrid.beta = 1:0.1:4;
else
    searchGrid.alpha = -3.5:.01:0;
    searchGrid.beta = 10:1:40;
end
searchGrid.gamma = linspace(0,guessUpper,3);
searchGrid.lambda = linspace(0,lapseUpper,6);

% Set up lapse limits.
guessLimits = [0 guessUpper];
lapseLimits = [0 lapseUpper];

% Set up standard options for Palamedes search.
options = PAL_minimize('options');

% Fit with Palemedes Toolbox.
[psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
    'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
    'checkLimits',false);

% Check whether beta is bigger or smaller than we would like,
% constrain it if so.
if (psiParamsValues(2) > max(searchGrid.beta))
    searchGrid.beta = max(searchGrid.beta);
    paramsFree = [1 0 1 1];
    [psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
        'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
        'checkLimits',false);
    % disp('Overrun on beta')
    % psiParamsValues
elseif (psiParamsValues(2) < min(searchGrid.beta))
    searchGrid.beta = min(searchGrid.beta);
    paramsFree = [1 0 1 1];
    [psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
        'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
        'checkLimits',false);
    % disp('Underunrun on beta')
    % psiParamsValues
end

% Get fitted curve values.
fit_psychometric = PF(psiParamsValues,fit_log_intensity);

% Calculate threshold: the difference between the stimulus levels at
% performances equal to 0.7602 and 0.5.
log_threshold = PF(psiParamsValues,threshold_criterion,'inverse');

% Now correct for guessing and lapse rate
psiParamsValuesCorrect = psiParamsValues;
psiParamsValuesCorrect(3) = 0;
psiParamsValuesCorrect(4) = 0;
corrected_fit_psychometric = PF(psiParamsValuesCorrect,fit_log_intensity);
corrected_log_threshold = PF(psiParamsValuesCorrect,threshold_criterion,'inverse');

end