%% Make Figure 3 for the paper

%% Clear
clear; close all;

%% Output variant
outputVariant = 'SlopeFree1';

%% Set up for the figures and read data
FigureSetup;

%% Get all MOCS data for full sessions
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'MOCS') & strcmp(dataTable.Split,'All'));
            if (sum(index) ~= 1)
                error('Have not properly set up condition to pick out just one sensitivity');
            end
            sensitivityMOCS(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);

            % Get sensitivity confidence limits. The table gives threshold, so the high CI value
            % is the lower limit and vice-versa when we are thinking about sensitivity.
            CILower_MOCS(pp,dd,ss) = -dataTable.CIHigh(index); %negative sensitivity, so high value corresonds to lower CI)
            CIUpper_MOCS(pp,dd,ss) = -dataTable.CILow(index);  %negative sensitivity, so low value corresonds to higher CI)

            % Negative error bar.  We want the magnitude of the downgoing error bar for
            % the MOCS sensitivities.  I think this should be 
            %   xneg(pp,dd,ss) = sensitivityMOCS(pp,dd,ss) - CILower_MOCS(pp,dd,ss);   
            % But here this value is the one for the positive error bar below.
            xneg(pp,dd,ss) = sensitivityMOCS(pp,dd,ss) - CIUpper_MOCS(pp,dd,ss);

            % Positive error bar.  We want the magnitude of the upgoing error bar.  I
            % think this should be 
            %   xpos(pp,dd,ss) = CIUpper_MOCS(pp,dd,ss) - sensitivityMOCS(pp,dd,ss);
            % But here this is the negative of the value for the negative error bar.
            %
            % Note. The errorbar() plotting function appears to take the absolute value of
            % its error bar magnitude arguments, so it probably doesn't hurt that the
            % values computed by the line below are negative, but cleaner to make them
            % positive as I have done in my suggested code above.
            xpos(pp,dd,ss) = sensitivityMOCS(pp,dd,ss) - CILower_MOCS(pp,dd,ss);  

            % David's fix (as per comments above).
            xneg(pp,dd,ss) = sensitivityMOCS(pp,dd,ss) - CILower_MOCS(pp,dd,ss);  
            xpos(pp,dd,ss) = CIUpper_MOCS(pp,dd,ss) - sensitivityMOCS(pp,dd,ss);
            if (any(xneg < 0) | any(xpos < 0))
                error('Check error bar magnitude calculation');
            end
        end
    end
end

%% Get all QUEST data for full sessions
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'QUEST') & strcmp(dataTable.Split,'All'));
            if (sum(index) ~= 1)
                error('Have not properly set up condition to pick out just one sensitivity');
            end
            sensitivityQUEST(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);

            % See comments where error bars are computed above
            CILower_QUEST(pp,dd,ss) = -dataTable.CIHigh(index); 
            CIUpper_QUEST(pp,dd,ss) = -dataTable.CILow(index);
            yneg(pp,dd,ss) = sensitivityQUEST(pp,dd,ss) - CIUpper_QUEST(pp,dd,ss);
            ypos(pp,dd,ss) = sensitivityQUEST(pp,dd,ss)- CILower_QUEST(pp,dd,ss);

            % David's fix 
            yneg(pp,dd,ss) = sensitivityQUEST(pp,dd,ss) - CILower_QUEST(pp,dd,ss);  
            ypos(pp,dd,ss) = CIUpper_QUEST(pp,dd,ss) - sensitivityQUEST(pp,dd,ss);
            if (any(yneg < 0) | any(ypos < 0))
                error('Check error bar magnitude calculation');
            end
        end
    end
end

%% Make Figure 3a

%% For Paper figure convert the diameters to arcmin^2 in label because we used arcmin^2 for stimulus size throughout the paper
theDiameters = [1.05,30.4]; % 8 pixels corrrespond to 1.05arcmin^2 and 43 pixels correspond to 30.4arcmin^2
plotSize = [100 100 400 400];
f = figure('Position',plotSize); clf; hold on

errorbar(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1),yneg(:,1,1),ypos(:,1,1),xneg(:,1,1),xpos(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2),yneg(:,1,2),ypos(:,1,2),xneg(:,1,2),xpos(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1),yneg(:,2,1),ypos(:,2,1),xneg(:,2,1),xpos(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2),yneg(:,2,2),ypos(:,2,2),xneg(:,2,2),xpos(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);

plot([limMin limMax],[limMin limMax],'k:','LineWidth',2);
xlabel('MOCS Sensitivity (dB)', 'FontWeight','bold', 'FontSize', 10, 'FontName', 'Times New Roman');
ylabel('QUEST Sensitivity (dB)','FontWeight','bold', 'FontSize', 10, 'FontName', 'Times New Roman');
set(gca, 'FontName', 'Arial')
legend( ...
    {sprintf('Session %d, %.2f arcmin^2',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %.2f arcmin^2',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %.1f arcmin^2',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %.1f arcmin^2',theSessions(2),theDiameters(2)) ; ...
    sprintf('Line of Equality')
    }, ...
    'Location','SouthEast');
axis('square');

lgd = legend('show');
legend boxoff
lgd.FontSize = 8; 
ax=gca;
ax.XAxis.FontWeight = 'bold'; % Make x-axis tick labels bold
ax.YAxis.FontWeight = 'bold'; % Make y-axis tick labels bold
% set(gca, 'FontName', 'Arial','FontWeight','bold')
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
axis([limMin limMax limMin limMax]);
ticks = 16:2:32; 
xticks(ticks);
yticks(ticks);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure3a.png'),'-dpng','-r600');

%% t-tests
[~,p(1,1)] = ttest(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1));
[~,p(1,2)] = ttest(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2));
[~,p(2,1)] = ttest(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1));
[~,p(2,2)] = ttest(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2));
fprintf('MOCS vs QUEST t-test p values\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.2f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

%% Wilcoxon signed-rank test
[p(1,1),h,~] = signrank(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1));
[p(1,2),h,~] = signrank(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2));
[p(2,1),h,~] = signrank(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1));
[p(2,2),h,~] = signrank(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2));
fprintf('MOCS vs QUEST t-test p values-Wilcoxon Test\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.2f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

%% Figure 3b (Bland Altman plot),  MOCS VS QUEST
%
% Recall that the second and third indices are subject, size (8 and 43), session (1 and 2)
Session1_8pixels_M = sensitivityMOCS(:,1,1);
Session2_8pixels_M = sensitivityMOCS(:,1,2);
Session1_43pixels_M = sensitivityMOCS(:,2,1);
Session2_43pixels_M = sensitivityMOCS(:,2,2);
Session1_8pixels_Q = sensitivityQUEST(:,1,1);
Session2_8pixels_Q = sensitivityQUEST(:,1,2);
Session1_43pixels_Q = sensitivityQUEST(:,2,1);
Session2_43pixels_Q = sensitivityQUEST(:,2,2);

% Print out mean and standard errors by session/method/size for reporting in the paper
fprintf('\n');
fprintf('MOCS : Session1, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session1_8pixels_M),std(Session1_8pixels_M)/sqrt(length(Session1_8pixels_M)));
fprintf('QUEST : Session1, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session1_8pixels_Q),std(Session1_8pixels_Q)/sqrt(length(Session1_8pixels_Q)));
fprintf('MOCS : Session2, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session2_8pixels_M),std(Session2_8pixels_M)/sqrt(length(Session2_8pixels_M)));
fprintf('QUEST : Session2, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session2_8pixels_Q),std(Session2_8pixels_Q)/sqrt(length(Session2_8pixels_Q)));
fprintf('MOCS : Session1, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session1_43pixels_M),std(Session1_43pixels_M)/sqrt(length(Session1_43pixels_M)));
fprintf('QUEST : Session1, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session1_43pixels_Q),std(Session1_43pixels_Q)/sqrt(length(Session1_43pixels_Q)));
fprintf('MOCS : Session2, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session2_43pixels_M),std(Session2_43pixels_M)/sqrt(length(Session2_43pixels_M)));
fprintf('QUEST : Session2, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(Session2_43pixels_Q),std(Session2_43pixels_Q)/sqrt(length(Session2_43pixels_Q)));

% Calcs for Bland-Altmann 
[mean_Session1_8pixels, diff_Session1_8pixels] = calculate_bland_altman(Session1_8pixels_M, Session1_8pixels_Q);
[mean_Session2_8pixels, diff_Session2_8pixels] = calculate_bland_altman(Session2_8pixels_M, Session2_8pixels_Q);
[mean_Session1_43pixels, diff_Session1_43pixels] = calculate_bland_altman(Session1_43pixels_M, Session1_43pixels_Q);
[mean_Session2_43pixels, diff_Session2_43pixels] = calculate_bland_altman(Session2_43pixels_M, Session2_43pixels_Q);

mean_diff_Session1_8pixels = mean(diff_Session1_8pixels);
std_diff_Session1_8pixels = std(diff_Session1_8pixels);
mean_diff_Session2_8pixels = mean(diff_Session2_8pixels);
std_diff_Session2_8pixels = std(diff_Session2_8pixels);
mean_diff_Session1_43pixels = mean(diff_Session1_43pixels);
std_diff_Session1_43pixels = std(diff_Session1_43pixels);
mean_diff_Session2_43pixels = mean(diff_Session2_43pixels);
std_diff_Session2_43pixels = std(diff_Session2_43pixels);

% (3b) Session1 8 Pixels 
LoA_Session1_8pixels = [mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels];
h = figure; 
set(h, 'Units', 'inches');
set(h, 'Position', [2, 2, 1.8, 3]);
scatter(mean_Session1_8pixels, diff_Session1_8pixels, 25, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S1');% Plot Bland-Altman data for each comparison
line([18, 23], [mean_diff_Session1_8pixels, mean_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS8, S1)');% Plot mean differences and limits of agreement
line([18, 23], [LoA_Session1_8pixels(2), LoA_Session1_8pixels(2)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8, S1)');
line([18, 23], [LoA_Session1_8pixels(1), LoA_Session1_8pixels(1)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 10)

ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
daspect([4 1 1]);
xlabel({'Mean Sensitivity (dB)'}'');
ylabel('MOCS - QUEST Sensitivitity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure3b.png'),'-dpng','-r600');

% (3c) Session 2 8 pixels
LoA_Session2_8pixels = [mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels,mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels];
h = figure; 
set(h, 'Units', 'inches');
set(h, 'Position', [2, 2, 1.8, 3]);
scatter(mean_Session2_8pixels, diff_Session2_8pixels, 25, 'blue', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S2');
line([18, 23], [mean_diff_Session2_8pixels, mean_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS8, S2)');
line([18, 23], [LoA_Session2_8pixels(2), LoA_Session2_8pixels(2)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS8, S2)');
line([18, 23], [LoA_Session2_8pixels(1), LoA_Session2_8pixels(1)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS8, S2)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 10)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
daspect([4 1 1]);
xlabel({'Mean Sensitivity (dB)'}'');
ylabel('MOCS - QUEST Sensitivitity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure3c.png'),'-dpng','-r600');

% (3d) Session 1 43 pixels
LoA_Session1_43pixels = [mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels,mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels];
h = figure; 
set(h, 'Units', 'inches');
set(h, 'Position', [2, 2, 1.8, 3]);
scatter(mean_Session1_43pixels, diff_Session1_43pixels, 25, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S1');
line([26, 30], [mean_diff_Session1_43pixels, mean_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS43, S1)');
line([26, 30], [LoA_Session1_43pixels(2), LoA_Session1_43pixels(2)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43, S1)');
line([26, 30], [LoA_Session1_43pixels(1), LoA_Session1_43pixels(1)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 10)
ax.XLim = [24, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
daspect([4 1 1]);
ticks = 24:2:32; 
xticks(ticks);
xlabel({'Mean Sensitivity (dB)'}'');
ylabel('MOCS - QUEST Sensitivitity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure3d.png'), '-dpng', '-r600');

% (3e) Session 2 43 pixels
LoA_Session2_43pixels = [mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels,mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels];
h = figure; 
set(h, 'Units', 'inches');
set(h, 'Position', [2, 2, 1.8, 3]);
scatter(mean_Session2_43pixels, diff_Session2_43pixels, 25, 'red', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S2');
line([26, 30], [mean_diff_Session2_43pixels, mean_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS43, S2)');
line([26, 30], [LoA_Session2_43pixels(2), LoA_Session2_43pixels(2)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS43, S2)');
line([26, 30], [LoA_Session2_43pixels(1), LoA_Session2_43pixels(1)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS43, S2)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 10)
ax.XLim = [24, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 10;
ax.YAxis.FontSize = 10;
daspect([4 1 1]);
ticks = 24:2:32; 
xticks(ticks);
xlabel({'Mean Sensitivity (dB)'}'');
ylabel('MOCS - QUEST Sensitivitity (dB)');
set(gcf,'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure3e.png'),'-dpng','-r600');

%% t-test for 8 Vs 43 pixel stimulus
[~,p(1,1)] = ttest(Session1_8pixels_M,Session1_43pixels_M);
[~,p(1,2)] = ttest(Session1_8pixels_Q,Session1_43pixels_Q);
[~,p(2,1)] = ttest(Session2_8pixels_M,Session2_43pixels_M);
[~,p(2,2)] = ttest(Session2_8pixels_Q,Session2_43pixels_Q);
theSplits= {'MOCS','QUEST'}; % to compare 43 Vs 8 stimulus sizes for two testing paradigms for two sessions
fprintf('\nWithin-Session (8 Vs 43 pixels) t-test p values\n');
for dd = 1:length(theSplits)
    for ss = 1:length(theSessions)
        fprintf('\t%s , session %d, p = %0.4f\n',theSplits{dd},theSessions(ss),p(dd,ss));
    end
end

%% Print limits of agreement
fprintf('\n');
fprintf('MOCS Vs QUEST : Session1, 8 pixels, LoA: %.2f, %.2f\n', LoA_Session1_8pixels);
fprintf('MOCS Vs QUEST : Session1, 43 pixels, LoA: %.2f, %.2f\n', LoA_Session1_43pixels);
fprintf('MOCS Vs QUEST : Session2, 8 pixels, LoA: %.2f, %.2f\n', LoA_Session2_8pixels);
fprintf('MOCS Vs QUEST : Session2, 43 pixels, LoA: %.2f, %.2f\n', LoA_Session2_43pixels);
