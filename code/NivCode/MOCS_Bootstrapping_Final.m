%MOCS

%FIT Logistic function
%Remove Catch Trials
Sorted_Intensities = sortrows(all_trials_unpacked,1);
Sorted_responses = Sorted_Intensities(:,2);
guess = find(Sorted_Intensities == -3.5);
guess_rate = sum(Sorted_responses(guess))/length(guess);

% Sorted_Intensities(guess,:)=[];
% Sorted_responses(guess,:)=[];
%Plot figure
minIntensityLog=-3.5;
maxIntensityLog= 0;
fontsize = 14; markersize = 6; fwidth = 350; fheight = 350;
% f0 = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
% xlabel('Log trial intensity (au)','FontSize',fontsize);
% ylabel('Percent seen (%)','FontSize',fontsize);
% xlim([minIntensityLog maxIntensityLog]);
% ylim([0 100]);
% set(a0,'FontSize',fontsize);
% set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
% set(a0,'Color','none');
% set(f0,'Color',[1 1 1]);
% set(f0,'PaperPositionMode','auto');
% set(f0, 'renderer', 'painters');
PF = @PAL_Logistic;
paramsFree = [1 1 0 0];
searchGrid.alpha = -3.5:.01:0;
searchGrid.beta = 3.5;
searchGrid.gamma = 0;   % Scalar here (since fixed) but may be vector
searchGrid.lambda = 0;  % Ditto % threshold, slope, guess rate, lapse rate


%Parameters for bootstrapping
[C,IA,IC] = unique(Sorted_Intensities(:,1),'legacy');
% logIntensityLevels=C';
numpresented=zeros(size(C));
for i = 1:length(numpresented)-1
    numpresented(1) = IA(1);
    numpresented(i+1) = IA(i+1)-IA(i);
end

numseen = nan(size(numpresented));
numseen(1)= sum(Sorted_responses(1:IA(1)));
for i = 2 :length(IA)
    numseen(i) = sum(Sorted_responses(IA(i-1)+1:IA(i)));
end


numbootstraps = 500;
f1 = figure;
stimSizes=43;
    set(f1, 'Units', 'centimeters')
    numPlotRows = 2; numPlotCols = ceil(length(stimSizes)./numPlotRows);
    fheight = 5*numPlotRows;
    fwidth = (fheight./numPlotRows).*numPlotCols;
    set(f1, 'Position', [1 5 fwidth fheight])
    set(gca,'Color','none');
    set(f1,'Color',[1 1 1]);
    set(f1,'PaperPositionMode','auto');
    set(f1, 'renderer', 'painters');
    
% Start Bootstrapping
for k = 1:numbootstraps + 1    
   if k ~= numbootstraps + 1
    [sampledInt, ind] = datasample(Sorted_Intensities(:,1), size(Sorted_Intensities,1), 'Replace', true);
    resp = Sorted_Intensities(:,2);
    sampledRes = resp(ind);
    
% numPlotRows=2;
% numPlotCols=1;
logIntensityLevels = sampledInt;
numSeen = sampledRes;
numPresented = ones(size(Sorted_responses));
%Correction for guessing
guess = find(logIntensityLevels == -3.5);
PFa = sum(numSeen(guess))/length(guess);
% PHit = sum(numSeen)/length(numSeen);
% PCorr = (PHit-PFa)/(1-PFa);

logIntensityLevels(guess,:)=[];
numSeen(guess,:)=[];


if PFa > 0
    searchGrid.gamma = PFa;
else 
    searchGrid.gamma = 0;
end

n=1;
[paramsValues, ~, exitFlag] = PAL_PFML_Fit(logIntensityLevels, numSeen, numPresented, ...
                                searchGrid, paramsFree, PF, 'lapseLimits', [0 0.05], 'guessLimits', [0 0.05]);

stimLevelsCoarse = linspace(-4, 1.5, 1000);
propCorrectPlot = PF(paramsValues, stimLevelsCoarse);
propCorrectPlot_corrected = (propCorrectPlot-PFa)/(1-PFa); %Correction for guessing
figure(f1)
hold on
subplot(1, 1, 1), plot(stimLevelsCoarse,propCorrectPlot_corrected,'-','color',[0 .7 0],'linewidth',3); hold on
threshPercent = 78;
fitThresh(n,1) = PF(paramsValues, threshPercent./100, 'inverse'); % to compute the value of y at x use inverse, also can used derivative

if fitThresh(n,1) < -4 || fitThresh(n,1) > 1
  fitThresh(n,1) = NaN;
end           
% bootstrapThresholds(i,1)=fitThresh();

 % Plot thresholds and error estimates on fits to real data
          if i == numbootstraps + 1  
                subplot(1, 2, 1), plot([fitThresh(n)-1.96*nanstd(bootstrapThresholds(:,n)) fitThresh(n)+1.96*nanstd(bootstrapThresholds(:,n))], ...
                    [fitThresh(n) fitThresh(n)],'-','color',[0 .7 0],'linewidth',1.5)
          end
if i ~= numbootstraps+1
    bootstrapThresholds(k,:) = fitThresh;
end         


% if i == numbootstraps+1
%       plot([fitThresh(n)-1.96*nanstd(bootstrapThresholds(:,n)) fitThresh(n)+1.96*nanstd(bootstrapThresholds(:,n))], ...
%       [fitThresh(n) fitThresh(n)],'-','color',[0 .7 0],'linewidth',1.5)
% end

% Plot stuff        
if isfinite(fitThresh(n,1))
        fontColor = [0 0 0];
else
        fontColor = [1 0 0];
end
title([{[num2str(stimSizes(n)) 'px']}, {['T_{fit} (' num2str(threshPercent) '%): ' num2str(fitThresh(n,1), '%.2f')]}], 'FontSize', 6, 'HorizontalAlignment', 'Center', 'Color', fontColor);
  end
end

Sorted_Intensities = sortrows(all_trials_unpacked,1);
guess_full = find(Sorted_Intensities(:,1) == -3.5);
PFa_full = sum(Sorted_responses(guess_full))/length(guess_full);
PHit_full= sum(Sorted_responses)/length(Sorted_responses);
PCorr_full = (PHit_full-PFa_full)/(1-PFa_full);
paramsFree1 = [1 1 0 0];
searchGrid1.alpha = -3.5:.01:0;
searchGrid1.beta = 3.5;
searchGrid1.gamma = 0;
searchGrid1.lambda = 0;

Sorted_Intensities(guess_full,:)=[];
Sorted_responses(guess_full,:)=[];

%Parameters for final plot
[C1,IA1,IC1] = unique(Sorted_Intensities(:,1),'legacy');
% logIntensityLevels=C';
numpresented_new=zeros(size(C1));
for i = 1:length(numpresented_new)-1
    numpresented_new(1) = IA1(1);
    numpresented_new(i+1) = IA1(i+1)-IA1(i);
end

numseen_new = nan(size(numpresented_new));
numseen_new(1)= sum(Sorted_responses(1:IA1(1)));


for i = 2 :length(IA1)
    numseen_new(i) = sum(Sorted_responses(IA1(i-1)+1:IA1(i)));
end

if PFa_full > 0
    searchGrid1.gamma= PFa_full;
else
    searchGrid1.gamma= 0;
end


[paramsValues1, negLL, exitFlag1] = PAL_PFML_Fit(Sorted_Intensities(:,1), Sorted_responses, ones(size(Sorted_responses)), ...
                                searchGrid1, paramsFree1, PF, 'lapseLimits', [0 0.05], 'guessLimits', [0 0.05]);
xEval = linspace(minIntensityLog, maxIntensityLog, 100);

propCorrectPlot1 = PF(paramsValues1, xEval);
propCorrectPlotfull_corrected = (propCorrectPlot1-PFa_full)/(1-PFa_full);
figure;plot(xEval, 100.*propCorrectPlotfull_corrected, '-', 'Color', [0 0.7 0], 'LineWidth', 2)
hold on;
% plot(C, 100.*(numseen./numpresented), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
threshPercent = 78;
fitThresh1(n,1) = PF(paramsValues1, threshPercent./100, 'inverse');

for p = 1:length(numpresented_new)
     q=1*numpresented_new(p)/6;
     plot(C1(p), 100.*(numseen_new(p)./numpresented_new(p)), 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', q)
     text(C1(p), 100.*(numseen_new(p)./numpresented_new(p)), num2str(numpresented_new(p)), 'color', 'r', 'FontSize', 10, 'HorizontalAlignment','left','VerticalAlignment', 'bottom', 'FontWeight', 'bold')
     hold on
end
x=fitThresh1; y=78;
plot(x,y,'bo','linewidth',.5)
hold on;
plot([x, x], [y, 0], '--k');  % Dotted line to y-axis
plot([x, -3.5], [y, y], '--k');  % Dotted line to x-axis
text(x,y,[num2str(fitThresh1,'%.2f')],'color', 'k', 'FontSize', 10, 'HorizontalAlignment','right','VerticalAlignment', 'cap', 'FontWeight', 'bold');
plot(x,y,'o-', 'LineWidth', 2, 'MarkerSize', 8)
%text(xEval,100.*PF(paramsValues1, xEval), num2str(fitThresh1))
% 
% plot(bootstrapThresholds,78,'xb','linewidth',5)
 
confidenceInterval = prctile(bootstrapThresholds, [5, 95]);
line([confidenceInterval(1), confidenceInterval(1)], ylim, 'Color', [0 0.7 0], 'LineWidth', 2);
line([confidenceInterval(2), confidenceInterval(2)], ylim, 'Color', [0 0.7 0], 'LineWidth', 2);