
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
theParticipants = {'11002'};

% Define sessions (1 or 2)
theSessions = [1 2];

% Define session splits
theSplits = {'All', 'FirstHalf', 'SecondHalf'};

% Define sizes (8 or 43)
theDiameters = [8 43];

% Define methods ('MOCS' or 'QUEST')
theMethods = {'MOCS' 'QUEST'};

%% Some parameters
log0Value = -3.5;
thresholdCriterion = 0.92;

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
    for mm = 1:length(theMethods)
        for dd = 1:length(theDiameters)
            for ss = 1:length(theSessions)
                for hh = 1:length(theSplits)
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
                    %
                    % The three is to skip over '.' and '..' which come first
                    dirOffset = 3;
                    trial_videos = dir(pathToData);
                    num_trial_videos = length(trial_videos(dirOffset:end));
                    if (num_trial_videos ~= 4)
                        error('Expect 4 data files');
                    end

                    % Read and concatenate the data
                    all_trials = {};
                    all_trials_unpacked = [];
                    for i = 1:num_trial_videos
                        % Read the data file
                        tempData = load(fullfile(pathToData,trial_videos(i+dirOffset-1).name));

                        % Grab the trial data we need
                        if (isfield(tempData,'trial_vector'))
                            if (~strcmp(theMethod{tableRow},'MOCS'))
                                error('Inconsistency in our understanding of method');
                            end
                            all_trials{i,1} = tempData.trial_vector;
                            all_trials{i,2} = tempData.response_vector;
                        elseif (isfield(tempData,'trial_matrix'))
                            if (~strcmp(theMethod{tableRow},'QUEST'))
                                error('Inconsistency in our understanding of method');
                            end
                            all_trials{i,1} = tempData.trial_matrix;
                            all_trials{i,2} = tempData.response_matrix;
                        else
                            error('Do not understand data format');
                        end

                        % Get/check number of trials
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

                    % Do the lookup table conversion
                    for i = 1:size(all_trials_unpacked,1)
                        % Get the nominal log intensity from the data
                        log_intensity_nominal(i) = all_trials_unpacked(i,1);

                        % Convert to nominal lin intensity
                        lin_intensity_nominal(i) = 10.^log_intensity_nominal(i);

                        % Round this.  Should match first column of lut
                        lin_intensity_rounded(i) = round(lin_intensity_nominal(i)*1000)/1000;

                        % Get the rounded log intensity, which we think is what
                        % is in the third column of the lut.  Except that if
                        % the linear intensity is 0, then the entry of the
                        % third column is zero and not -Inf.  Since we won't
                        % use the rounded log intensity for anything other than
                        % checking, we just tweak it here.
                        log_intensity_rounded(i) = log10(lin_intensity_rounded(i));
                        log_intensity_rounded_chk(i) = log_intensity_rounded(i);
                        if (lin_intensity_rounded(i) == 0)
                            log_intensity_rounded(i) = -3.5;
                            log_intensity_rounded_chk(i) = 0;
                        end

                        % Compute the row index into the lut based on the linear intensity
                        lut_row_index(i) = lin_intensity_rounded(i)*1000+1;

                        % Get the corrected linear intensity by indexing into
                        % the second column of the LUT
                        lin_intensity_lut(i) = AOM.green_AOM_lut(lut_row_index(i),2);

                        % Convert to log, taking log10(0) to be -3.5
                        log_intensity_lut(i) = log10(lin_intensity_lut(i));

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


                    % Correct for trials that were nominally not zero, but
                    % really were zero.
                    % NEED TO DO THIS

                    % Here is the format of all_trials_unpacked
                    % These values are not adjusted for power measurements.
                    % they are scaled with the instrument maximum having a
                    % linear intensity of 1.
                    %
                    %  Column 1 - nominal log10 intensity
                    %  Column 2 - 1 = seen, 0 = not seen
                    %  Column 3 - lut corrected linear intensity
                    %  Column 4 - lut corrected log intensity

                    % Split the data
                    nTrials = size(all_trials_unpacked,1);
                    shuffleIndex = Shuffle(1:nTrials);
                    switch (theSplit{tableRow})
                        case 'All'
                            dataIndex = 1:nTrials;
                        case 'FirstHalf'
                            dataIndex = shuffleIndex(1:round(nTrials/2));
                        case 'SecondHalf'
                            dataIndex = shuffleIndex(round(nTrials/2):nTrials);
                        otherwise
                            error('Unknown split specified');
                    end
                    all_trials_unpacked = all_trials_unpacked(dataIndex,:);

                    % Fit the psychometric function.  The fitting routine makes
                    % a plot and we adjust the title here for things the
                    % fitting routine doesn't know about.
                    [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
                        log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(all_trials_unpacked(:,4),all_trials_unpacked(:,2),log0Value,thresholdCriterion);
                    title({ sprintf('Subject %s, %s, session %d, split %s',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow}) ; ...
                        sprintf('Diameter %d arcmin, log threshold %0.2f',theDiameter(tableRow),corrected_log_threshold) ; ...
                        sprintf('Slope %0.1f, guess %0.2f, lapse %0.2f',psiParamsValues(2),psiParamsValues(3),psiParamsValues(4)) ; ''},'FontSize',20);
                    drawnow;
                    saveas(h,fullfile(pathToAnalysis,'psychometricFcn.tif'),'tif');

                    % Save what we learned
                    save(fullfile(pathToAnalysis,'fitOutput.mat'),'all_trials_unpacked','plot_log_intensities','plot_psychometric', ...
                        'corrected_psychometric','log_threshold','corrected_log_threshold','psiParamsValues','-v7.3');

                    % Accumulate data
                    theTableLogThreshold(tableRow,1) = log_threshold;
                    theTableCorrectedLogThreshold(tableRow,1) = corrected_log_threshold;
                    theTableCorrectedDbThreshold(tableRow,1) = 10*corrected_log_threshold;
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

    % Close figures this subject
    close all;

    % Make summary figures for this subject here
end

% Write out full analysis data table
tableVariableNames = {'Subject','Diameter','Session','Method','Split','Log10 Threshold','Corrected Log10 Threshold','Corrected Threshold (dB)', 'Alpha','Beta','Guess','Lapse','Criterion'};
dataTable = table(theSubject,theDiameter,theSession,theMethod,theSplit, ...
    round(theTableLogThreshold,2),round(theTableCorrectedLogThreshold,2), round(theTableCorrectedDbThreshold,1), ...
    round(theTablePsiParamsValues(:,1),2),round(theTablePsiParamsValues(:,2),1), round(theTablePsiParamsValues(:,3),2), round(theTablePsiParamsValues(:,4),2), ...
    theTableThresholdCriterion, ...
    'VariableNames',tableVariableNames);
writetable(dataTable,fullfile(analysisDir,'fitTable.xlsx'),'FileType','spreadsheet');


% Overall summary figures here


% FitPsychometricFunction
%
% Massage data, fit and plot psychometric function
function [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(log_intensity,responses,log0Value,thresholdCriterion)

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
    fitPsychometric(unique_log_intensities,num_seen',num_presented',thresholdCriterion);

% Plot
minIntensityLog=-3.5;
maxIntensityLog= 0;
fontsize = 18; fwidth = 600; fheight = 600;
h = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
set(gca,'FontName','Helvetical','FontSize',14);
xlabel('Log trial intensity (au)','FontSize',fontsize);
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
plot([-3.5 corrected_log_threshold],[thresholdCriterion thresholdCriterion],'k:','LineWidth',2);
plot([corrected_log_threshold corrected_log_threshold],[0 thresholdCriterion],'k:','LineWidth',2);

end

% fitPsychometric
%
% Wrapper into Palamedes that fits the logistic with parameters reasonable
% for our data.
function [fit_log_intensity,fit_psychometric,corrected_fit_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues] = fitPsychometric(log_intensity,num_pos,out_of_num,threshold_criterion)
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
fit_log_intensity = linspace(-3.5,1,1000);

% Psychometric function form (alternative: PAL_Gumbel).
PF = @PAL_Logistic;

% 'paramsFree' is a boolean vector that determines what parameters get
% searched over (1: free parameter, 0: fixed parameter).
paramsFree = [1 1 1 1];

% Range info
guessUpper = 0.10;
lapseUpper = 0.05;

% Set up starting points:
searchGrid.alpha = -3.5:.01:0;
searchGrid.beta = 10:1:40;
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

% Check whether beta is bigger than we would like,
% constrain it at max if so.
if (psiParamsValues(2) > max(searchGrid.beta))
    searchGrid.beta = max(searchGrid.beta);
    paramsFree = [1 0 1 1];
    [psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
        'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
        'checkLimits',false);
    % disp('Overrun on beta')
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