
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

% Define sizes (8 or 43)
theDiameters = [8 43];

% Define methods ('MOCS' or 'QUEST')
theMethods = {'MOCS' 'QUEST'};

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
for mm = 1:length(theMethods)
    theMethod = theMethods{mm};
    for pp = 1:length(theParticipants)
        theSubject = theParticipants{pp};
        for dd = 1:length(theDiameters)
            theDiameter = theDiameters(dd);
            for ss = 1:length(theSessions)
                theSession = theSessions(ss);

                % Form path the the directory with the data sitting in it
                pathToData = fullfile(dataDir,theSubject,['Session' num2str(theSession)],['Size' num2str(theDiameter)],theMethod);
                pathToAnalysis = fullfile(analysisDir,theSubject,['Session' num2str(theSession)],['Size' num2str(theDiameter)],theMethod);
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
                        if (~strcmp(theMethod,'MOCS'))
                            error('Inconsistency in our understanding of method');
                        end
                        all_trials{i,1} = tempData.trial_vector;
                        all_trials{i,2} = tempData.response_vector;
                    elseif (isfield(tempData,'trial_matrix'))
                        if (~strcmp(theMethod,'QUEST'))
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
                %  Column 1 - nominal log10 intensity
                %  Column 2 - 1 = seen, 0 = not seen
                %  Column 3 - lut corrected log10 intensity
                %  Column 4 - lut corrected linear intensity
                % These values are not adjusted for power measurements

                [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
                    log_threshold,corrected_log_threshold,h] = FitPsychometricFunction(all_trials_unpacked);
                title({ sprintf('Subject %s, %s, session %d, all data',theSubject,theMethod,theSession) ; ...
                    sprintf('Diameter %d arcim,og threshold %0.2f',theDiameter,corrected_log_threshold) ; ''},'FontSize',20);
                    
                % Save what we learned
                save(fullfile(pathToAnalysis,'all_trials_unpacked.mat'),'all_trials_unpacked','-v7.3');

                % Clear
                clear all_trials all_trials_unpacked

            end
        end
    end
end

% Fit the psychometric function
function [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold,h] = FitPsychometricFunction(all_trials_unpacked)

% Sort data and get what we need
sorted_data = sortrows(all_trials_unpacked,4);
sorted_log_intensities = sorted_data(:,4);
sorted_log_intensities(isinf(sorted_log_intensities)) = -3.5;
sorted_responses = sorted_data(:,2);

% Find raw guess rate
catch_trial_index = find(sorted_log_intensities == -3.5);
guess_rate = sum(sorted_responses(catch_trial_index))/length(catch_trial_index);

% Get unique stimulus levels and data for those
[unique_log_intensities] = unique(sorted_log_intensities,'legacy');
numpresented = zeros(size(unique_log_intensities));
for i = 1:length(numpresented)
    num_presented(i) = length(find(sorted_log_intensities == unique_log_intensities(i)));
    num_seen(i) = sum(sorted_responses(find(sorted_log_intensities == unique_log_intensities(i))));
end

% Fit using a wrapper into Palemedes
thresholdCriterion = 0.78;
[plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold] = ...
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

% 
% % Set up Palamedes logistic fit
% PF = @PAL_Logistic;
% paramsFree = [1 1 0 0];
% searchGrid.alpha = -3.5:.01:0;
% searchGrid.beta = 3.5;
% searchGrid.gamma = 0;   % Scalar here (since fixed) but may be vector
% searchGrid.lambda = 0;  % Ditto % threshold, slope, guess rate, lapse rate
% 
% %Parameters for bootstrapping
% [unique_log_intensities,IA,IC] = unique(sorted_log_intensities,'legacy');
% % logIntensityLevels=C';
% numpresented = zeros(size(unique_log_intensities));
% for i = 1:length(numpresented)
%     numpresented(i) = length(find(sorted_log_intensities == unique_log_intensities(i)));
%     numseen(i) = sum(sorted_responses(find(sorted_log_intensities == unique_log_intensities(i))));
% end
% 
% numbootstraps = 500;
% f1 = figure;
% stimSizes=43;
%     set(f1, 'Units', 'centimeters')
%     numPlotRows = 2; numPlotCols = ceil(length(stimSizes)./numPlotRows);
%     fheight = 5*numPlotRows;
%     fwidth = (fheight./numPlotRows).*numPlotCols;
%     set(f1, 'Position', [1 5 fwidth fheight])
%     set(gca,'Color','none');
%     set(f1,'Color',[1 1 1]);
%     set(f1,'PaperPositionMode','auto');
%     set(f1, 'renderer', 'painters');
% 
% % Start Bootstrapping
% for k = 1:numbootstraps + 1    
%    if k ~= numbootstraps + 1
%     [sampledInt, ind] = datasample(Sorted_Intensities(:,1), size(Sorted_Intensities,1), 'Replace', true);
%     resp = Sorted_Intensities(:,2);
%     sampledRes = resp(ind);
% 
% % numPlotRows=2;
% % numPlotCols=1;
% logIntensityLevels = sampledInt;
% numSeen = sampledRes;
% numPresented = ones(size(Sorted_responses));
% %Correction for guessing
% guess = find(logIntensityLevels == -3.5);
% PFa = sum(numSeen(guess))/length(guess);
% % PHit = sum(numSeen)/length(numSeen);
% % PCorr = (PHit-PFa)/(1-PFa);
% 
% logIntensityLevels(guess,:)=[];
% numSeen(guess,:)=[];
% 
% 
% if PFa > 0
%     searchGrid.gamma = PFa;
% else 
%     searchGrid.gamma = 0;
% end
% 
% n=1;
% [paramsValues, ~, exitFlag] = PAL_PFML_Fit(logIntensityLevels, numSeen, numPresented, ...
%                                 searchGrid, paramsFree, PF, 'lapseLimits', [0 0.05], 'guessLimits', [0 0.05]);
% 
% stimLevelsCoarse = linspace(-4, 1.5, 1000);
% propCorrectPlot = PF(paramsValues, stimLevelsCoarse);
% propCorrectPlot_corrected = (propCorrectPlot-PFa)/(1-PFa); %Correction for guessing
% figure(f1)
% hold on
% subplot(1, 1, 1), plot(stimLevelsCoarse,propCorrectPlot_corrected,'-','color',[0 .7 0],'linewidth',3); hold on
% threshPercent = 78;
% fitThresh(n,1) = PF(paramsValues, threshPercent./100, 'inverse'); % to compute the value of y at x use inverse, also can used derivative
% 
% if fitThresh(n,1) < -4 || fitThresh(n,1) > 1
%   fitThresh(n,1) = NaN;
% end           
% % bootstrapThresholds(i,1)=fitThresh();
% 
%  % Plot thresholds and error estimates on fits to real data
%           if i == numbootstraps + 1  
%                 subplot(1, 2, 1), plot([fitThresh(n)-1.96*nanstd(bootstrapThresholds(:,n)) fitThresh(n)+1.96*nanstd(bootstrapThresholds(:,n))], ...
%                     [fitThresh(n) fitThresh(n)],'-','color',[0 .7 0],'linewidth',1.5)
%           end
% if i ~= numbootstraps+1
%     bootstrapThresholds(k,:) = fitThresh;
% end         
% 
% 
% % if i == numbootstraps+1
% %       plot([fitThresh(n)-1.96*nanstd(bootstrapThresholds(:,n)) fitThresh(n)+1.96*nanstd(bootstrapThresholds(:,n))], ...
% %       [fitThresh(n) fitThresh(n)],'-','color',[0 .7 0],'linewidth',1.5)
% % end
% 
% % Plot stuff        
% if isfinite(fitThresh(n,1))
%         fontColor = [0 0 0];
% else
%         fontColor = [1 0 0];
% end
% title([{[num2str(stimSizes(n)) 'px']}, {['T_{fit} (' num2str(threshPercent) '%): ' num2str(fitThresh(n,1), '%.2f')]}], 'FontSize', 6, 'HorizontalAlignment', 'Center', 'Color', fontColor);
%   end
% end
% 
% Sorted_Intensities = sortrows(all_trials_unpacked,1);
% guess_full = find(Sorted_Intensities(:,1) == -3.5);
% PFa_full = sum(Sorted_responses(guess_full))/length(guess_full);
% PHit_full= sum(Sorted_responses)/length(Sorted_responses);
% PCorr_full = (PHit_full-PFa_full)/(1-PFa_full);
% paramsFree1 = [1 1 0 0];
% searchGrid1.alpha = -3.5:.01:0;
% searchGrid1.beta = 3.5;
% searchGrid1.gamma = 0;
% searchGrid1.lambda = 0;
% 
% Sorted_Intensities(guess_full,:)=[];
% Sorted_responses(guess_full,:)=[];
% 
% %Parameters for final plot
% [C1,IA1,IC1] = unique(Sorted_Intensities(:,1),'legacy');
% % logIntensityLevels=C';
% numpresented_new=zeros(size(C1));
% for i = 1:length(numpresented_new)-1
%     numpresented_new(1) = IA1(1);
%     numpresented_new(i+1) = IA1(i+1)-IA1(i);
% end
% 
% numseen_new = nan(size(numpresented_new));
% numseen_new(1)= sum(Sorted_responses(1:IA1(1)));
% 
% 
% for i = 2 :length(IA1)
%     numseen_new(i) = sum(Sorted_responses(IA1(i-1)+1:IA1(i)));
% end
% 
% if PFa_full > 0
%     searchGrid1.gamma= PFa_full;
% else
%     searchGrid1.gamma= 0;
% end
% 
% 
% [paramsValues1, negLL, exitFlag1] = PAL_PFML_Fit(Sorted_Intensities(:,1), Sorted_responses, ones(size(Sorted_responses)), ...
%                                 searchGrid1, paramsFree1, PF, 'lapseLimits', [0 0.05], 'guessLimits', [0 0.05]);
% xEval = linspace(minIntensityLog, maxIntensityLog, 100);
% 
% propCorrectPlot1 = PF(paramsValues1, xEval);
% propCorrectPlotfull_corrected = (propCorrectPlot1-PFa_full)/(1-PFa_full);
% figure;plot(xEval, 100.*propCorrectPlotfull_corrected, '-', 'Color', [0 0.7 0], 'LineWidth', 2)
% hold on;
% % plot(C, 100.*(numseen./numpresented), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
% threshPercent = 78;
% fitThresh1(n,1) = PF(paramsValues1, threshPercent./100, 'inverse');
% 
% for p = 1:length(numpresented_new)
%      q=1*numpresented_new(p)/6;
%      plot(C1(p), 100.*(numseen_new(p)./numpresented_new(p)), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', q)
%      text(C1(p), 100.*(numseen_new(p)./numpresented_new(p)), num2str(numpresented_new(p)), 'color', 'r', 'FontSize', 10, 'HorizontalAlignment','left','VerticalAlignment', 'bottom', 'FontWeight', 'bold')
%      hold on
% end
% x=fitThresh1; y=78;
% plot(x,y,'bo','linewidth',.5)
% hold on;
% plot([x, x], [y, 0], '--k');  % Dotted line to y-axis
% plot([x, -3.5], [y, y], '--k');  % Dotted line to x-axis
% text(x,y,[num2str(fitThresh1,'%.2f')],'color', 'k', 'FontSize', 10, 'HorizontalAlignment','right','VerticalAlignment', 'cap', 'FontWeight', 'bold');
% plot(x,y,'o-', 'LineWidth', 2, 'MarkerSize', 8)
% %text(xEval,100.*PF(paramsValues1, xEval), num2str(fitThresh1))
% % 
% % plot(bootstrapThresholds,78,'xb','linewidth',5)
% 
% confidenceInterval = prctile(bootstrapThresholds, [5, 95]);
% line([confidenceInterval(1), confidenceInterval(1)], ylim, 'Color', [0 0.7 0], 'LineWidth', 2);
% line([confidenceInterval(2), confidenceInterval(2)], ylim, 'Color', [0 0.7 0], 'LineWidth', 2); 

end


function [fit_log_intensity,fit_psychometric,corrected_fit_psychometric, ...
    log_threshold,corrected_log_threshold] = fitPsychometric(log_intensity,NumPos,OutOfNum,thresholdCriterion)
%fitPsychometric
%
% Usage:
%   [fit_log_intensity,fit_psychometric,corrected_fit_psychonmetric, ...
%        log_threshold,corrected_log_threshold] = fitPsychometric(x,NumPos,OutOfNum,thresholdCriterion)
%
% Inputs:
%   x        : (vector) stimulus levels tested
%   NumPos   : (vector) for each stimulus level tested, the number of trials a positive response was given
%   OutOfNum : (vector) for each stimulus level tested, the total number of trials
%   thresholdCriterion : [scalar] criterion for threshold
%
% History:
%   07/05/21  amn  Wrote it.
%   06/28/25  dhb  Adopted for this project

% Calculate x-axis values to plot
%
% Increase the number of stimulus level values to plot for a smooth fit.
fit_log_intensity = linspace(min(log_intensity),max(log_intensity),1000);

% Psychometric function form (alternative: PAL_Gumbel).
PF = @PAL_Logistic;

% 'paramsFree' is a boolean vector that determines what parameters get
% searched over (1: free parameter, 0: fixed parameter).
paramsFree = [1 1 1 1];  

% Set up starting points:
%   1st (mean of the cumulative normal): set to mean of stimulus levels tested
%   2nd (standard deviation): set to a fraction of the range of stimulus levels tested
searchGrid = [mean(log_intensity) 1/(max(log_intensity)-min(log_intensity)) 0 0];

% Set up lapse limits.
guessLimits = [0 0.05];
lapseLimits = [0 0.05];

% Set up standard options for Palamedes search.
options = PAL_minimize('options');

% Fit with Palemedes Toolbox.
[paramsValues] = PAL_PFML_Fit(log_intensity,NumPos,OutOfNum,searchGrid,paramsFree,PF, ...
    'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false);

% Get fitted curve values.
fit_psychometric = PF(paramsValues,fit_log_intensity);

% Calculate threshold: the difference between the stimulus levels at
% performances equal to 0.7602 and 0.5.
log_threshold = PF(paramsValues,thresholdCriterion,'inverse');

% Now correct for guessing and lapse rate
paramsValues(3) = 0;
paramsValues(4) = 0;
corrected_fit_psychometric = PF(paramsValues,fit_log_intensity);
corrected_log_threshold = PF(paramsValues,thresholdCriterion,'inverse');

end