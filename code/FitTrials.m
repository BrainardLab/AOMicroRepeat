%% FitTrials
%
% Read in the combined data file and fit it all.

%% Initialize
close all; clear all;

%% Name of output path for the fits
outputName = 'SlopeFree';
switch (outputName)
    case 'SlopeFree'
        if (convertToDb)
            slopeLower = 0.5;
            slopeUpper = 7;
        else
            slopeLower = 5;
            slopeUpper = 70;
        end
    case 'SlopeFixed'
        if (convertToDb)
            slopeLower = 3.5;
            slopeUpper = 3.5;
        else
            slopeLower = 35;
            slopeUpper = 35;
        end
end


%% Figure visibility
%
% Maybe making invisible figures will avoid crashes?
figureVis = 'off';

%% Read combined data produced by CombineData
analysisDir = getpref('AOMicroRepeat','analysisDir');
d = load(fullfile(analysisDir,'combinedData.mat'),'all_trials','all_trials_unpacked','log0Value','theParticipants','theDiameters','theSessions','theSplits','theMethods','AOM');
all_trials_unpacked = d.all_trials_unpacked;
log0Value = d.log0Value;
theParticipants = d.theParticipants;
theDiameters = d.theDiameters;
theSessions = d.theSessions;
theSplits = d.theSplits;
theMethods = d.theMethods;
AOM = d.AOM;

%% Some parameters
guessUpper = 0.05;
lapseUpper = 0.05;
log0TrialThreshold = -3;
thresholdCriterion = 0.78;
nBootstraps = 500;
convertToDb = true;
if (convertToDb)
    slopeLower = 3.5;
    slopeUpper = 3.5;
else
    slopeLower = 35;
    slopeUpper = 35;
end

%% Freeze rng seed for repeatability
rng(101);

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

                    % If just checking, break and don't do anything else
                    fprintf('\tLUT correcting trial data\n');

                    % Clear out variables from previous time through
                    % the loop.
                    clear log_intensity_nominal lin_intensity_nominal lin_intensity_rounded log_intensity_nominal log_intensity_rounded log_intensity_rounded_chk
                    clear lut_row_index dac_intensity_lut lin_intensity_lut
                    % Do the lookup table conversion
                    fprintf('\tEntering LUT correction loop over trials ... ')
                    for i = 1:size(all_trials_unpacked{pp,dd,ss,hh,mm},1)
                        % Dir for this output
                        pathToAnalysis = fullfile(analysisDir,outputName,theSubject{tableRow},['Session' num2str(theSession(tableRow))],['Size' num2str(theDiameter(tableRow))],theMethod{tableRow},theSplit{tableRow});
                        if (~exist(pathToAnalysis,"dir"))
                            mkdir(pathToAnalysis);
                        end

                        % Get the nominal log intensity from the data
                        log_intensity_nominal(i) = all_trials_unpacked{pp,dd,ss,hh,mm}(i,1);

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
                        lutRowIndex = find(AOM.green_AOM_lut(:,4) == dac_intensity_lut(i) );
                        if (isempty(lutRowIndex))
                            error('No LUT rows correspond to trial dac intensity');
                        end
                        lin_intensity_lut(i) = mean(AOM.green_AOM_lut(lutRowIndex,1));

                        % Convert to log, taking log10(0) to be the
                        % log0Value
                        if (lin_intensity_lut(i) == 0)
                            log_intensity_lut(i) = log0Value;
                        else
                            log_intensity_lut(i) = log10(lin_intensity_lut(i));
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
                        all_trials_unpacked{pp,dd,ss,hh,mm}(i,3) = lin_intensity_lut(i);
                        all_trials_unpacked{pp,dd,ss,hh,mm}(i,4) = log_intensity_lut(i);
                    end
                    fprintf('\tdone\n');

                    % Figure to make sure LUT conversion is basically
                    % the identity.  Uncomment if you want to look at
                    % this.
                    h1 = figure('visible',figureVis); subplot(1,2,1); hold on;
                    plot(lin_intensity_nominal,lin_intensity_lut,'ro');
                    axis('square'); axis([0 1 0 1]);
                    plot([0 1],[0 1],'k:');
                    title({ sprintf('%s, %s, session %d, split %s, d %d',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow},theDiameter(tableRow)) ; ...
                                       ''},'FontSize',14);
                    subplot(1,2,2); hold on;
                    plot(lin_intensity_nominal,lin_intensity_lut,'ro');
                    axis('square');
                    if (theDiameters(dd) == min(theDiameters))
                        axis([0 0.1 0 0.1]);
                        plot([0 0.1],[0 0.1],'k:');
                    else
                        axis([0 0.01 0 0.01]);
                        plot([0 0.01],[0 0.01],'k:');
                    end
                    title({ sprintf('%s, %s, session %d, split %s, d %d',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow},theDiameter(tableRow)) ; ...
                        ''},'FontSize',14);
                    drawnow;
                    saveas(h1,fullfile(pathToAnalysis,'quantizationFig.tif'),'tif');

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
                    fprintf('\tSplitting data as specified\n')
                    nTrials = size(all_trials_unpacked{pp,dd,ss,hh,mm},1);
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
                    all_trials_unpacked{pp,dd,ss,hh,mm} = all_trials_unpacked{pp,dd,ss,hh,mm}(dataIndex,:);
                    nTrials = size(all_trials_unpacked{pp,dd,ss,hh,mm},1);

                    % Extract the core data to fit
                    trial_log_intensities = all_trials_unpacked{pp,dd,ss,hh,mm}(:,4);
                    trial_responses = all_trials_unpacked{pp,dd,ss,hh,mm}(:,2);

                    % Convert to dB?
                    if (convertToDb)
                        trial_log_intensities = 10*trial_log_intensities;
                        log0ValueUse = 10*log0Value;
                    end

                    % Staircase plot.  We only do this for full sessions.
                    % We take advantage of the fact that we know the data
                    % were concatenated run-by-run.
                    switch (theSplit{tableRow})
                        case 'All'
                            if (strcmp(theMethod{tableRow},'MOCS') | strcmp(theMethod{tableRow},'QUEST'))
                                nRuns = 4;
                                nTrialsPerRun = length(trial_responses)/nRuns;
                                if (round(nTrialsPerRun) ~= nTrialsPerRun)
                                    error('Do not understand how data were concatenated');
                                end
                                for zz = 1:nRuns
                                    run_log_intensities = trial_log_intensities((zz-1)*nTrialsPerRun+1:zz*nTrialsPerRun);
                                    run_responses = trial_responses((zz-1)*nTrialsPerRun+1:zz*nTrialsPerRun);
                                    hs = figure('Visible',figureVis); hold on
                                    set(gca,'FontName','Helvetica','FontSize',14);
                                    runIndicator = 1:nTrialsPerRun;
                                    plot(runIndicator,run_log_intensities,'r','LineWidth',2);
                                    index = find(run_responses == 1);
                                    plot(runIndicator(index),run_log_intensities(index),'o','MarkerSize',8,'Color',[0.5 1 0.5],'MarkerFaceColor',[0.5 1 0.5]);
                                    index = find(run_responses == 0);
                                    plot(runIndicator(index),run_log_intensities(index),'^','MarkerSize',8,'Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1]);
                                    xlabel('Trial number','FontSize',18)
                                    if (convertToDb)
                                        ylabel('Trial intensity (dB)','FontSize',18);
                                        ylim([-36 1]);
                                    else
                                        ylabel('Log trial intensity (au)','FontSize',18);
                                        ylim([-3.6 0.1]);
                                    end
                                    xlim(([0 100]));
                                    title({ sprintf('Subject %s, %s, session %d, split %s',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow}) ; ...
                                        sprintf('Diameter %d arcmin, run %d',theDiameter(tableRow),zz) ; ''},'FontSize',20);
                                    drawnow;
                                    saveas(hs,fullfile(pathToAnalysis,sprintf('staircasePlot_run%d.tif',zz)),'tif');
                                    title('');
                                    drawnow;
                                    saveas(hs,fullfile(pathToAnalysis,sprintf('staircasePlotNoTitle_run%d.tif',zz)),'tif');
                                end
                            end
                    end              

                    % Fit the psychometric function.  The fitting routine makes
                    % a plot and we adjust the title here for things the
                    % fitting routine doesn't know about.
                    fprintf('\tCalling top level PF fit routine\n');
                    [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
                        log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(trial_log_intensities,trial_responses,log0ValueUse,thresholdCriterion, ...
                        convertToDb,true,figureVis,guessUpper,lapseUpper,slopeLower,slopeUpper);
                    fprintf('\tDone with top level PF fit\n');
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
                    fprintf('\tBootstrapping data\n')
                    nTrialsForBootstrap = length(trial_log_intensities);
                    bootstrap_log_intensities = zeros(nTrialsForBootstrap,nBootstraps);
                    bootstrap_responses = zeros(nTrialsForBootstrap,nBootstraps);
                    switch(theMethod{tableRow})
                        case {'QUEST', 'COMBINED'}
                            % For QUEST and COMBINED, we draw with replacement from all
                            % the trials except the catch trials, which we draw from separately.
                            for bb = 1:nBootstraps
                                temp_boot_intensities = [];
                                temp_boot_responses = [];

                                % Catch trials
                                uindex = find(trial_log_intensities == log0ValueUse);
                                utemp_intensities = trial_log_intensities(uindex);
                                utemp_responses = trial_responses(uindex);
                                bootIndex = randi(length(uindex),length(uindex),1);
                                temp_boot_intensities = [temp_boot_intensities ; utemp_intensities(bootIndex)];
                                temp_boot_responses = [temp_boot_responses ; utemp_responses(bootIndex)];

                                % Non catch trials
                                uindex = find(trial_log_intensities ~= log0ValueUse);
                                utemp_intensities = trial_log_intensities(uindex);
                                utemp_responses = trial_responses(uindex);
                                bootIndex = randi(length(uindex),length(uindex),1);
                                temp_boot_intensities = [temp_boot_intensities ; utemp_intensities(bootIndex)];
                                temp_boot_responses = [temp_boot_responses ; utemp_responses(bootIndex)];

                                bootstrap_log_intensities(:,bb) = temp_boot_intensities;
                                bootstrap_responses(:,bb) = temp_boot_responses;
                            end
                            clear lutRowIndex

                        case 'MOCS'
                            % For MOCS, we respect the experimental design
                            % and bootstrap each stimulus intensity separately.
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
                    fprintf('\tFitting bootstrapped data\n');
                    for bb = 1:nBootstraps
                        if (rem(bb,100) == 0)
                            fprintf('\t\tBootstrap fit %d of %d\n',bb,nBootstraps);
                        end
                        [~,~,~,~,boot_corrected_threshold(bb),~,~] = FitPsychometricFunction(bootstrap_log_intensities(:,bb),bootstrap_responses(:,bb),log0ValueUse,thresholdCriterion, ...
                            convertToDb,false,false,guessUpper,lapseUpper,slopeLower,slopeUpper);
                    end
                    fprintf('\tDone fitting bootstrapped data\n');
                    boot_quantiles = quantile(boot_corrected_threshold,[0.025 0.5 0.975]);
                    % figure(h);
                    plot([boot_quantiles(1) boot_quantiles(3)],[thresholdCriterion thresholdCriterion],'b','LineWidth',4);
                    saveas(h,fullfile(pathToAnalysis,'psychometricFcnCI.tif'),'tif');

                    % Save a version without our informative title, for
                    % paper figures
                    title('');
                    saveas(h,fullfile(pathToAnalysis,'psychometricFcnCINoTitle.tif'),'tif');

                    % Save what we learned
                    fprintf('\tSaving ... ');
                    this_all_trials_unpacked = all_trials_unpacked{pp,dd,ss,hh,mm};
                    save(fullfile(pathToAnalysis,'fitOutput.mat'),'this_all_trials_unpacked','plot_log_intensities','plot_psychometric', ...
                        'corrected_psychometric','log_threshold','corrected_log_threshold','psiParamsValues','boot_corrected_threshold','boot_quantiles','-v7.3');
                    fprintf('done\n');

                    % Accumulate data
                    theTableLogThreshold(tableRow,1) = log_threshold;
                    theTableCorrectedLogThreshold(tableRow,1) = corrected_log_threshold;
                    theTableCorrectedLogThresholdBootMedian(tableRow,1)= boot_quantiles(2);
                    theTableCorrectedBootCILow(tableRow,1) = boot_quantiles(1);
                    theTableCorrectedBootCIHigh(tableRow,1) = boot_quantiles(3);
                    theTablePsiParamsValues(tableRow,:) = psiParamsValues;
                    theTableThresholdCriterion(tableRow,1) = thresholdCriterion;

                    % Clear and close
                    clear bootstrap_log_intensities bootstrap_responses boot_corrected_threshold boot_quantiles
                    clear dac_intensity_lut lin_intensity_lut lin_intensity_nominal lin_intensity_rounded
                    clear log_intensity_lut log_intensity_nominal log_intensity_rounded log_intensity_rounded_chk lut_row_index
                    clear plot_log_intensities plot_psychometric shuffleIndex
                    clear temp_boot_responses temp_boot_responses trial_log_intensities trial_responses

                    % Bump table row
                    tableRow = tableRow + 1;
                end
            end
        end
    end

    % Close figures this subject
    close all hidden;

    % Make any summary figures for this subject here

end

% Write out full analysis data table
fprintf('\nWriting xlsx file ... ')
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
fprintf('done\n');

% FitPsychometricFunction
%
% Massage data, fit and plot psychometric function
function [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(log_intensity,responses, ...
    log0Value,thresholdCriterion,convert_to_db,make_plot,figure_vis, ...
    guessUpper,lapseUpper,slopeLower,slopeUpper)

% Sort data and get what we need
sorted_data = sortrows([log_intensity responses],1);
sorted_log_intensities = sorted_data(:,1);
index = find(isinf(sorted_log_intensities));
if (~isempty(index))
    sorted_log_intensities(isinf(sorted_log_intensities)) = log0Value;
end
sorted_responses = sorted_data(:,2);

% Get unique stimulus levels and data for those
[unique_log_intensities] = unique(sorted_log_intensities,'legacy');
numpresented = zeros(size(unique_log_intensities));
for i = 1:length(numpresented)
    num_presented(i) = length(find(sorted_log_intensities == unique_log_intensities(i)));
    num_seen(i) = sum(sorted_responses(find(sorted_log_intensities == unique_log_intensities(i))));
end

% Fit using a wrapper into Palemedes
if (make_plot)
    fprintf('\t\tCalling Palemedes wrapper ... ');
end
[plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold, psiParamsValues] = ...
    fitPsychometric(unique_log_intensities,num_seen',num_presented',thresholdCriterion,convert_to_db,guessUpper,lapseUpper,slopeLower,slopeUpper);
if (make_plot)
    fprintf('done\n');
end

% Plot
if (make_plot)
    if (convert_to_db)
        maxIntensityLog= 1;
        minIntensityLog=-36;
    else
        maxIntensityLog= 0.1;
        minIntensityLog=-3.6;
    end
    fontsize = 18; fwidth = 600; fheight = 600;
    h = figure('Position', [400 200 fwidth fheight],'visible',figure_vis); hold on;
    a0 = gca;
    set(a0,'FontName','Helvetica','FontSize',14);
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

