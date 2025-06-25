% Bootstrapping analysis of spatial summation data. If desired, this script
% will remove mis-delivered trials (for data collected with eye tracking)
% from spatial summation data and fit psychometric functions for each
% stimulus size via a maximum likelihood approach (i.e. using the Palamedes
% Toolbox). This process is repeated some number of times to generate error
% estimates for the fitted parameters -- specifically, thresholds and
% Ricco's areas derived from a IxA vs A plot). Bootstrap samples of the
% data are drawn with replacement. Because QUEST was used to collect this
% data, the distribution of test intensities is biased heavily towards the
% threshold percent seeing (78%) on the psychometric function; to keep the
% fits sensible, we constrain the resampling procedure from
% equally-populated bins along the stimulus intensity dimension.

% 4/20/2017     wst modified an earlier version of this script 
%7/20/2017      wst added new Ricco's fitting routine
%               (@riccosFittingFunction) which should hopefully speed
%               things up a bit

%% Housekeeping
% clear all
% close all
% clc
% load all_trials_unpacked
Sorted_Data = sortrows(all_trials_unpacked,1); 
intensity  =  Sorted_Data(:,1);
response  =  Sorted_Data(:,2);
threshold = Sorted_Data(:,3);
number_nan = isnan(threshold);
catch_trials = find(number_nan==1);
guess = find(response(catch_trials))==1;
guess_rate = length(guess)/length(number_nan);

%% Psychometric function fitting parameters
PF = @PAL_Logistic; % Fit a logistic function
%Threshold and Slope are free parameters, guess and lapse rate are fixed
paramsFree = [1 1 0 0]; %1: free parameter, 0: fixed parameter
%searchGrid.alpha = -6:.01:6;
searchGrid.alpha = -4:.01:0;
searchGrid.beta = 3.5;
searchGrid.gamma = 0;   % Scalar here (since fixed) but may be vector
searchGrid.lambda = 0;  % Ditto


% if isempty(guess)== true
%     guess_rate = 0;
%     intensity(catch_trials,:)=[];
%     response(catch_trials,:)=[];
% elseif isempty(guess)== false
%   
%     searchGrid.gamma=  guess_rate; % Update guess rate if the response is 1 during a catch_trial
%     intensity(catch_trials,:)=[];
%     response(catch_trials,:)=[];
% end
%% The main bootstrapping loop
    % Number of simulations to run for this folder
    numBootstrap = 500;
    % Can save a video of the bootstrapping procedure here if you want
    % (only on first folder)
    % Convert stimulus sizes from pixels to arcmin
    stimSizes = 1; %(here just to represent the number of stimulus sizes used)
    % Pre-allocate
    fitThresh = zeros(length(stimSizes),1);
    
    % Plotting so that you can observe the bootstrapping in action
    % PF fits go in Figure 1
    f1 = figure;
    set(f1, 'Units', 'centimeters')
%     fontsize = 14; markersize = 10; fwidth = 550; fheight = 350;
%     f0 = figure('Position', [400 200 fwidth fheight]); 
    numPlotRows = 2; numPlotCols = 2;
    fheight = 10*numPlotRows;
    fwidth = 10.*numPlotCols;
    set(f1, 'Position', [1 5 fwidth fheight])
    set(gca,'Color','none');
    set(f1,'Color',[1 1 1]);
    set(f1,'PaperPositionMode','auto');
    set(f1, 'renderer', 'painters');
    fontsize = 14; markersize = 6; 
    xlim([-4 0]);
    ylim([0 1]);
    xlabel('Trial intensity (\DeltaI)', 'FontSize', fontsize+2)
    ylabel('Proportion seen', 'FontSize', fontsize+2)
    
%     A waitbar is nice
%     h = waitbar(0);
%     set(h, 'Units', 'centimeters');
%     set(h,'Position', [15 18 9.5 2])
    
    % Pre-allocate
    bootstrapThresholds = nan(numBootstrap, length(stimSizes));
    bootstrapData ={};
    
    for i = 1:numBootstrap+1 % Real data are fit and plotted at the end
        % Cycle through stimulus sizes
            i
            logStimLevels = intensity;
            numSeen = response==1;
            numPresented = ones(size(numSeen));
            catch_threshold = threshold;

            find_nan = isnan (catch_threshold);
            index_nan = find(find_nan==1);
            catch_threshold(index_nan)=0;
            

            if i ~=numBootstrap+1
                % Bootstrapping stuff here
                testMatrix = [logStimLevels numSeen numPresented,catch_threshold];
                testMatrix = sortrows(testMatrix,1); % Sort these
                numPerBin = 5; % Bin the data for slightly constrained resampling
                numBins = ceil(length(logStimLevels)/numPerBin);
                
                if numBins ~= length(logStimLevels)/numPerBin % If the length of the data is not a multiple of the population per bin
                    testMatrix(length(logStimLevels)+1:numPerBin*numBins,:) = NaN; % Pad out the matrix so you can reshape it without crashing
                end
                
                stimLevelsReshape = reshape(testMatrix(:,1), [], numBins); % Each column is a bin from which you can bootstrap
                numSeenReshape = reshape(testMatrix(:,2), [], numBins); % Each column is a bin from which you can bootstrap
                numPresentedReshape = reshape(testMatrix(:,3), [], numBins); % Each column is a bin from which you can bootstrap
                numthreshReshape = reshape(testMatrix(:,4), [], numBins);

                samplingMatrix = randi([1 numPerBin], numPerBin, numBins);% resample the test Matrix
                stimLevelsSampled = nan(size(samplingMatrix)); % Pre-allocate
                numSeenSampled = stimLevelsSampled; % Pre-allocate
                numPresentedSampled = stimLevelsSampled; % Pre-allocate
                numthreshSampled = stimLevelsSampled;

                for j = 1:size(samplingMatrix,2)
                    stimLevelsSampled(:,j) = stimLevelsReshape(samplingMatrix(:,j),j);
                    numSeenSampled(:,j) = numSeenReshape(samplingMatrix(:,j),j);
                    numPresentedSampled(:,j) = numPresentedReshape(samplingMatrix(:,j),j);
                    numthreshSampled(:,j) = numthreshReshape(samplingMatrix(:,j),j);
                end
                 
                % Reassemble into vectors
                StimLevels = reshape(stimLevelsSampled,[],1);
                numSeen = reshape(numSeenSampled,[],1);
                numPresented = reshape(numPresentedSampled,[],1);
                numthresh = reshape(numthreshSampled,[],1);
                % Remove NaNs
                StimLevels = StimLevels(isfinite(StimLevels));
                numSeen = numSeen(isfinite(numSeen));
                numPresented = numPresented(isfinite(numPresented));
                numthresh = numthresh(isfinite(numthresh));
                bootstrapData{i,1} = StimLevels; 

                %Correction for guessing
                catch_trials_bootstrapping = find(numthresh==0);
                %guess_bootstrapping = find(response(catch_trials))==1;
                PFa = sum(numSeen(catch_trials_bootstrapping))/length(catch_trials_bootstrapping);
                % PHit = sum(numSeen)/length(numSeen);
                % PCorr = (PHit-PFa)/(1-PFa);

                StimLevels(catch_trials_bootstrapping,:)=[];
                numSeen(catch_trials_bootstrapping,:)=[];
                numPresented(catch_trials_bootstrapping,:)=[];

                if PFa > 0
                  searchGrid.gamma = PFa;
                else 
                  searchGrid.gamma = 0;
                end

                [paramsValues, ~, ~] = PAL_PFML_Fit((StimLevels),numSeen, ...
                numPresented,searchGrid,paramsFree,PF);
                
            end            
        
            % Now do the fitting using Palamedes toolbox
      
            n=1;
            % Coarse evaluation of the fit for plotting purposes
            stimLevelsCoarse = linspace(-4, 0, 1000);
            propCorrectPlot = PF(paramsValues, stimLevelsCoarse);
            propCorrectPlot_bootstrapping = (propCorrectPlot-PFa)/(1-PFa); %Correction for guessing
               % Find point in the fit that gives the QUEST threshold percentage;
            threshPercent = 78;
            fitThresh(n,1) = PF(paramsValues, threshPercent./100, 'inverse');
            if fitThresh(n,1) < -4 || fitThresh(n,1) > 1
                fitThresh(n,1) = NaN;
            end
            if i ~= numBootstrap+1
                 bootstrapThresholds(i,:) = fitThresh;
            end  


            figure(f1)
            hold on
            plot(stimLevelsCoarse,propCorrectPlot_bootstrapping,'-','color',[0 .7 0],'linewidth',1); 
            hold off
            fontsize = 14; markersize = 6; fwidth = 850; fheight = 650;
            minIntensityLog=-4;
            maxIntensityLog= 0;
           
            
            if isfinite(fitThresh(n,1))
                fontColor = [0 0 0];
            else
                fontColor = [1 0 0];
            end
            title([{[num2str(stimSizes(n)) 'px']}, {['T_{fit} (' num2str(threshPercent) '%): ' num2str(fitThresh(n,1), '%.2f')]}], 'FontSize', 10, 'HorizontalAlignment', 'Center', 'Color', fontColor);
            hold off
     end
            f0 = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
            xlabel('Log trial intensity (au)','FontSize',fontsize);
            ylabel('Percent seen (%)','FontSize',fontsize);
            xlim([minIntensityLog maxIntensityLog]);
            ylim([0 100]);
            set(a0,'FontSize',fontsize);
            set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
            set(a0,'Color','none');
            set(f0,'Color',[1 1 1]);
            set(f0,'PaperPositionMode','auto');
            set(f0, 'renderer', 'painters');
            figure(f0)
            hold on

            searchGrid.alpha = -4:.01:0;
            searchGrid.beta = 3.5;
            searchGrid.gamma = 0;   % Scalar here (since fixed) but may be vector
            searchGrid.lambda = 0;
             %Correction for guessing
             guess_full = catch_trials;
             PFa_full = guess_rate;
                % PHit = sum(numSeen)/length(numSeen);
                % PCorr = (PHit-PFa)/(1-PFa);
               
                int_full = intensity;
                resp_full = response;
                int_full(guess_full,:)=[];
                resp_full(guess_full,:)=[];


                if PFa_full > 0
                  searchGrid.gamma = PFa_full;
                else 
                  searchGrid.gamma = 0;
                end
           
            [paramsValues1, ~, exitFlag] = PAL_PFML_Fit(int_full, resp_full, ones(size(resp_full)), ...
            searchGrid, paramsFree, PF, 'lapseLimits', [0 0.05], 'guessLimits', [0 0.05]);
            threshPercent = 78;fitThresh1(n,1) = PF(paramsValues1, threshPercent./100, 'inverse');
            xEval = linspace(minIntensityLog, maxIntensityLog, 100);
            propCorrectPlot1 = PF(paramsValues1, xEval);
            propCorrectPlotfull_corrected = (propCorrectPlot1-PFa_full)/(1-PFa_full);
            plot(xEval,100.*propCorrectPlotfull_corrected, '-', 'Color', [0 0.7 0], 'LineWidth', 2); hold on

%             subplot(1, 1, 1), 
            % Plot thresholds and error estimates on fits to real data
%             if i == numBootstrap+1
%                 subplot(2, 1, 1), 
%                 plot([fitThresh(n)-1.96*nanstd(bootstrapThresholds(:,n)) fitThresh(n)+1.96*nanstd(bootstrapThresholds(:,n))], ...
%                     [fitThresh(n) fitThresh(n)],'-','color',[0 .7 0],'linewidth',1.5)
%             end
    
            % Plot stuff
%            stimLevels = logStimLevels;

            intensity(catch_trials,:)=[];
            response(catch_trials,:)=[];
            [N, edges] = histcounts(intensity, 'BinMethod', 'auto');
            markersize = 10;
            for j = 1:length(N)
                if N(j)~=0
                    %subplot(numPlotRows, numPlotCols, n), 
                    
                    hold on, plot(mean(intensity(intensity>=edges(j) & intensity<=edges(j+1))), sum(response(intensity>=edges(j) & intensity<=edges(j+1)))*100./N(j), ...
                        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', markersize)
                    %subplot(numPlotRows, numPlotCols, n), 
                    hold on, text(mean(intensity(intensity>=edges(j) & intensity<=edges(j+1))), (sum(response(intensity>=edges(j) & intensity<=edges(j+1)))*100./N(j)), num2str(N(j)), 'Color', [1 1 1], ...
                        'VerticalAlignment', 'middle','HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
                end
            end
   
            % Plot area shading to represent specified quantiles of the
%             % bootstrapping
%             qPercentLow = 5;
%             qPercentHigh = 95;
%             bootstrapQuantiles = log10(quantile(bootstrapThresholds, [qPercentLow/100 qPercentHigh/100]));
          
            % Plot the data +/- bootstrap simulations here
          
        
%         if i == numBootstrap+1
%             errorbar(xData, yData, 1.96.*nanstd(bootstrapThresholds), 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k', 'CapSize', 0), hold on
%         end
        
       
%         
%      hold on;
%      num_presented_new=ones(size(response));
%      numseen = response;
% for p = 1:length(ones(size(response)))
%      q=1*num_presented_new(p)/4
%      plot(logintensityLevels(p), 100.*(numseen(p)./numpresented(p)), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', q)
%      text(logintensityLevels(p), 100.*(numseen(p)./numpresented(p)), num2str(num_presented_new(p)), 'color', 'r', 'FontSize', 10, 'HorizontalAlignment','right','VerticalAlignment', 'cap', 'FontWeight', 'bold')
%      hold on
% end
% 
x=fitThresh1; y=78;
fontsize = 14; markersize = 6; fwidth = 750; fheight = 550;
plot(x,y,'bo','linewidth',.5)
hold on;
plot([x, x], [y, 0], '--k');  % Dotted line to y-axis
plot([x, -4], [y, y], '--k');  % Dotted line to x-axis
text(x,y,[num2str(fitThresh1,'%.2f')],'color', 'k', 'FontSize', 10, 'HorizontalAlignment','left','VerticalAlignment', 'cap', 'FontWeight', 'bold');
% plot(x,y,'o-', 'LineWidth', 2, 'MarkerSize', 8)
%text(xEval,100.*PF(paramsValues1, xEval), num2str(fitThresh1))
% 
% plot(bootstrapThresholds,78,'xb','linewidth',5)
 
confidenceInterval = prctile(bootstrapThresholds, [5, 95]);
line([confidenceInterval(1), confidenceInterval(1)], ylim, 'Color', [0 0.7 0], 'LineWidth', 1);
line([confidenceInterval(2), confidenceInterval(2)], ylim, 'Color', [0 0.7 0], 'LineWidth', 1);

          

