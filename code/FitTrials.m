%% FitTrials
%
% Read in the combined data file and fit it all.

%% Initialize
close all hidden; clear all;

%% Variant
%
% The variant is a string that determines a set of fitting parameters.  Output for each
% variant is stored separately.  As an example, maybe we want to allow the slope of the
% psychometric function to be a free parameter versus fix it to a particular value, or
% maybe we want to experiment with different limits on the guess or lapse rate.  The
% variant string allows you to go look at how paremeters for each variant were set and
% then go find the output for that variant.
outputVariant = 'SlopeFree1';

%% Some parameters
guessUpper = 0.05;
log0TrialThreshold = -3;
thresholdCriterion = 0.78;
nBootstraps = 500;
convertToDb = true;

% This parameter determins how much above threshold a trial intensity
% needs to be for us to use it to estimate lapse rate. Different numerical
% values depending on whether you are using log10 units or dB. Indeed, this sort of
% numerical adjustment happens in other places below.
if (convertToDb)
    lapseEstIntensityThreshold = 6.0;
else
    lapseEstIntensityThreshold = 0.6;
end

%% Determine other parameters based on variant
%
% Note that locking lapse rate to 0 for slope fixed leads
% to some pathalogical fits, for reasons I don't fully understand.
% Finesse by using 1%, which is small enough for practical
% purposes.
%
% To get fixed slope values, first run with 'SlopeFree', and then
% use SummarizePFSlopes to get useful info, including the
% mean slopes from the slope free fits.
switch (outputVariant)
    case 'SlopeFree'
        if (convertToDb)
            slopeLower8 = 0.5;
            slopeUpper8 = 7;
            slopeLower43 = 0.5;
            slopeUpper43 = 7;
        else
            slopeLower8 = 5;
            slopeUpper8 = 70;
            slopeLower43 = 5;
            slopeUpper43 = 70;
        end
        lapseUpper = 0.05;
    case 'SlopeFree1'
        if (convertToDb)
            slopeLower8 = 0.5;
            slopeUpper8 = 7;
            slopeLower43 = 0.5;
            slopeUpper43 = 7;
        else
            slopeLower8 = 5;
            slopeUpper8 = 70;
            slopeLower43 = 5;
            slopeUpper43 = 70;
        end
        lapseUpper = 0.01;
    case 'SlopeFixed'
        if (convertToDb)
            slopeLower8 = 1.08;
            slopeUpper8 = 1.08;
            slopeLower43 = 2.69;
            slopeUpper43 = 2.69;
        else
            slopeLower8 = 10.8;
            slopeUpper8 = 10.8;
            slopeLower43 = 26.9;
            slopeUpper43 = 26.9;
        end
        lapseUpper = 0.01;
end

%% Figure visibility
%
% Maybe making invisible figures will avoid crashes?
% They still get saved out so easy to look at later.
figureVis = 'off';

%% Read combined data produced by CombineTrials
%
% See comments in CombineTrials about how to set preferences for where data and analsis 
% directories live.  Do not hard code the paths.
% setpref('AOMicroRepeat','dataDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\New_analysis_20250912\dataDir');
% setpref('AOMicroRepeat','analysisDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\figure-rerun');
analysisDir = getpref('AOMicroRepeat','analysisDir');
d = load(fullfile(analysisDir,'combinedData.mat'),'all_trials','all_trials_unpacked','log0Value','theParticipants','theDiameters','theSessions','theSplits','theMethods','AOM');
all_trials_unpacked = d.all_trials_unpacked;
log0Value = d.log0Value;
theParticipants = d.theParticipants;
theDiameters = d.theDiameters;
theSessions = d.theSessions;
theMethods = d.theMethods;
AOM = d.AOM;

% Check and handle how theSplits comes in and goes
% out.  This is tricky handshaking between this routine
% and CombineTrials.  Despite d.theSplits saying there
% is just one split, there are 3, with the latter two being
% the session based splits done according to the pre-registration.
%
% The way things are coded is a bit kluged up, but we check like
% mad that our various assumptions hold.
if (length(d.theSplits) > 1 | ~strcmp(d.theSplits,'All'))
    error('Unexpected form of d.theSplits. Did you re-run CombineData?');
end
theSplits = {'All' 'Group1' 'Group2'};

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
                        pathToAnalysis = fullfile(analysisDir,outputVariant,theSubject{tableRow},['Session' num2str(theSession(tableRow))],['Size' num2str(theDiameter(tableRow))],theMethod{tableRow},theSplit{tableRow});
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
                    plot(lin_intensity_nominal,lin_intensity_lut,'ro','MarkerFaceColor','r');
                    axis('square'); axis([0 1 0 1]);
                    plot([0 1],[0 1],'k:');
                    xlabel('Linear intensity nominal');
                    ylabel('Linear intensity LUT');
                    title({ sprintf('%s, %s, session %d, split %s, size %d',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow},theDiameter(tableRow)) ; ...
                        ''},'FontSize',10);
                    subplot(1,2,2); hold on;
                    plot(lin_intensity_nominal,lin_intensity_lut,'ro','MarkerFaceColor','r');
                    axis('square');
                    if (theDiameters(dd) == min(theDiameters))
                        axis([0 0.1 0 0.1]);
                        plot([0 0.1],[0 0.1],'k:');
                    else
                        axis([0 0.01 0 0.01]);
                        plot([0 0.01],[0 0.01],'k:');
                    end
                    xlabel('Linear intensity nominal');
                    ylabel('Linear intensity LUT');
                    title({ sprintf('%s, %s, session %d, split %s, size %d',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow},theDiameter(tableRow)) ; ...
                        ''},'FontSize',10);
                    drawnow;
                    saveas(h1,fullfile(pathToAnalysis,'quantizationFig.tif'),'tif');

                    h2 = figure('visible',figureVis); hold on;
                    if (convertToDb)
                        plot(10*log_intensity_nominal,10*log_intensity_lut,'ro','MarkerFaceColor','r');
                        plot([-36 1],[-36 1],'k:');
                        xlabel('Intensity nominal (dB)');
                        ylabel('Intensity LUT (dB)');
                        axis([-36 1 -36 1]);
                    else
                        plot(log_intensity_nominal,log_intensity_lut,'ro','MarkerFaceColor','r');
                        plot([-3.6 0.1],[-3.6 0.1],'k:');
                        xlabel('Log10 intensity nominal');
                        ylabel('Log10 intensity LUT');
                        axis([-3.6 0.1 -3.6 0.1]);
                    end
                    title({ sprintf('%s, %s, session %d, split %s, size %d',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow},theDiameter(tableRow)) ; ...
                        ''},'FontSize',10);
                    axis('square');
                    drawnow;
                    print(gcf,fullfile(pathToAnalysis,'quantizationLogFig.png'), '-dpng', '-r600'); % saves current figure as PNG at 600 dpi

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
                    trial_log_intensities = all_trials_unpacked{pp,dd,ss,hh,mm}(:,4);
                    trial_responses = all_trials_unpacked{pp,dd,ss,hh,mm}(:,2);

                    % Convert to dB?
                    if (convertToDb)
                        trial_log_intensities = 10*trial_log_intensities;
                        log0ValueUse = 10*log0Value;
                    end

                    % Staircase plot.  We only do this for full sessions and only for MOCS
                    % and QUEST, not COMBINED.
                    %
                    % We take advantage of the fact that we know the data
                    % were concatenated run-by-run.  Maybe not the most secure coding, but
                    % OK at least for now.
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
                                    fontsize = 14; fwidth = 3; fheight = 3;
                                    hs = figure('Units', 'inches','Position', [400 200 fwidth fheight]); hold on;
                                    set(gca,'FontName','Arial','FontSize',14);
                                    ax.LineWidth = 2;
                                    ax.LineWidth = 2;
                                    runIndicator = 1:nTrialsPerRun;
                                    plot(runIndicator,run_log_intensities,'r','LineWidth',2);
                                    index = find(run_responses == 1);
                                    plot(runIndicator(index),run_log_intensities(index),'o','MarkerSize',4,'Color',[0.5 1 0.5],'MarkerFaceColor',[0.5 1 0.5]);
                                    index = find(run_responses == 0);
                                    plot(runIndicator(index),run_log_intensities(index),'^','MarkerSize',4,'Color',[0.5 0.5 1],'MarkerFaceColor',[0.5 0.5 1]);
                                    xlabel('Trial number','FontSize',14)
                                    if (convertToDb)
                                        ylabel('Trial intensity (dB)','FontSize',14);
                                        ylim([-35 0]);
                                    else
                                        ylabel('Log trial intensity (au)','FontSize',14);
                                        ylim([-3.6 0.1]);
                                    end
                                    xlim(([0 100]));
                                    title({ sprintf('Subject %s, %s, session %d, split %s',theSubject{tableRow},theMethod{tableRow},theSession(tableRow),theSplit{tableRow}) ; ...
                                        sprintf('Diameter %d arcmin, run %d',theDiameter(tableRow),zz) ; ''},'FontSize',20);
                                    drawnow;
                                    saveas(hs,fullfile(pathToAnalysis,sprintf('staircasePlot_run%d.tif',zz)),'tif');
                                    title('');
                                    drawnow;
                                    print(hs, fullfile(pathToAnalysis,sprintf('staircasePlotNoTitle_run%d.png',zz)), '-dpng', '-r600'); % saves current figure as PNG at 600 dpi

                                end
                            end
                    end

                    % Set beta for this size
                    switch (theDiameter(tableRow))
                        case 8
                            slopeLowerUse(tableRow,1) = slopeLower8;
                            slopeUpperUse(tableRow,1) = slopeUpper8;
                        case 43
                            slopeLowerUse(tableRow,1) = slopeLower43;
                            slopeUpperUse(tableRow,1) = slopeUpper43;
                        otherwise
                    end
                    guessUpperUse(tableRow,1) = guessUpper;
                    lapseUpperUse(tableRow,1) = lapseUpper;

                    % Fit the psychometric function.  The fitting routine makes
                    % a plot and we adjust the title here for things the
                    % fitting routine doesn't know about.
                    fprintf('\tCalling top level PF fit routine\n');
                    [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
                        log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(trial_log_intensities,trial_responses,log0ValueUse,thresholdCriterion, ...
                        convertToDb,true,figureVis,guessUpperUse(tableRow),lapseUpperUse(tableRow),slopeLowerUse(tableRow),slopeUpperUse(tableRow));
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
                    title('');
                    print(h, fullfile(pathToAnalysis,sprintf('psychometricFcn_NoTitle.png',zz)), '-dpng', '-r600'); 

                    % Estimate guess and lapse rates from the data
                    index = find(trial_log_intensities == log0ValueUse);
                    if (~isempty(index))
                        guessEstFromData(tableRow,1) = mean(trial_responses(index));
                    else
                        guessEstFromData(tableRow,1) = -1;
                    end
                    nTrialsGuessEstFromData(tableRow,1) = length(index);
                    index = find(trial_log_intensities > corrected_log_threshold + lapseEstIntensityThreshold);
                    if (~isempty(index))
                        lapseEstFromData(tableRow,1) = 1-mean(trial_responses(index));
                    else
                        lapseEstFromData(tableRow,1) = -1;
                    end
                    nTrialsLapseEstFromData(tableRow,1) = length(index);

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
                            convertToDb,false,false,guessUpperUse(tableRow),lapseUpperUse(tableRow),slopeLowerUse(tableRow),slopeUpperUse(tableRow));
                    end
                    fprintf('\tDone fitting bootstrapped data\n');
                    boot_quantiles = quantile(boot_corrected_threshold,[0.025 0.5 0.975]);
                    % figure(h);
                    plot([boot_quantiles(1) boot_quantiles(3)],[thresholdCriterion thresholdCriterion],'b','LineWidth',4);
                    saveas(h,fullfile(pathToAnalysis,'psychometricFcnCI.tif'),'tif');

                    % Save a version without our informative title, for
                    % paper figures
                    title('');
                    print(gcf,fullfile(pathToAnalysis,'psychometricFcnCINoTitle.png'), '-dpng', '-r600'); % saves current figure as PNG at 600 dpi

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

    % Could make summary figures for this subject here

end

% Write out full analysis data table
fprintf('\nWriting xlsx file ... ')
if (convertToDb)
    tableVariableNames = {'Subject','Diameter','Session','Method','Split','Threshold (dB)','Corrected Threshold (dB)','CI Low','CI High', 'Alpha','Beta','Guess','Lapse', ...
        'Criterion','SlopeLimitLow','SlopeLimitHigh','guessLimit','lapseLimit','guessEstFromData','nTrialsGuessEstFromData','lapseEstFromData','nTrialsLapseEstFromData'};
    thresholdPlaces = 1;
else
    tableVariableNames = {'Subject','Diameter','Session','Method','Split','Log10 Threshold','Corrected Log10 Threshold','CI Low','CI High', 'Alpha','Beta','Guess','Lapse', ...
        'Criterion','SlopeLimitLow','SlopeLimitHigh','guessLimit','lapseLimit','guessEstFromData','nTrialsGuessEstFromData','lapseEstFromData','nTrialsLapseEstFromData'};
    thresholdPlaces = 2;
end
dataTable = table(theSubject,theDiameter,theSession,theMethod,theSplit, ...
    round(theTableLogThreshold,thresholdPlaces),round(theTableCorrectedLogThreshold,thresholdPlaces),round(theTableCorrectedBootCILow,thresholdPlaces), round(theTableCorrectedBootCIHigh,thresholdPlaces), ...
    round(theTablePsiParamsValues(:,1),2),round(theTablePsiParamsValues(:,2),2), round(theTablePsiParamsValues(:,3),3), round(theTablePsiParamsValues(:,4),3), ...
    theTableThresholdCriterion, ...
    slopeLowerUse,slopeUpperUse,guessUpperUse,lapseUpperUse, ...
    round(guessEstFromData,3), nTrialsGuessEstFromData,round(lapseEstFromData,3),nTrialsLapseEstFromData, ...
    'VariableNames',tableVariableNames);
writetable(dataTable,fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet');
fprintf('done\n');

