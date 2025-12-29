% Figure 4
clear; close all;

%% Output variant
outputVariant = 'SlopeFree1';

%% Set up for the figures and read data
FigureSetup;

%%  Get Group1 data for combined trial (all participants, 2 sessions, and 2 stimulus sizes)
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'COMBINED') & strcmp(dataTable.Split,'Group1'));
            if (sum(index) ~= 1)
                error('Have not properly set up condition to pick out just one sensitivity');
            end

            % See comments in Figure3.m about this.  David's fix implemented here.
            sensitivityGroup1(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
            CILower_Group1(pp,dd,ss)= -dataTable.CIHigh(index); 
            CIUpper_Group1(pp,dd,ss)= -dataTable.CILow(index);  
            xneg(pp,dd,ss) = sensitivityGroup1(pp,dd,ss) - CILower_Group1(pp,dd,ss);
            xpos(pp,dd,ss) = CIUpper_Group1(pp,dd,ss) - sensitivityGroup1(pp,dd,ss);
            if (any(xneg < 0) | any(xpos < 0))
                error('Check error bar magnitude calculation');
            end
        end
    end
end

%%  Get Group2 data for combined trial (all participants, 2 sessions, and 2 stimulus sizes)
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'COMBINED') & strcmp(dataTable.Split,'Group2'));
            if (sum(index) ~= 1)
                error('Have not properly set up condition to pick out just one sensitivity');
            end

            % See comments in Figure3.m about this.  David's fix implemented here.
            sensitivityGroup2(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
            CILower_Group2(pp,dd,ss)= -dataTable.CIHigh(index); 
            CIUpper_Group2(pp,dd,ss)= -dataTable.CILow(index);
            yneg(pp,dd,ss) = sensitivityGroup2(pp,dd,ss) - CILower_Group2(pp,dd,ss);
            ypos(pp,dd,ss) = CIUpper_Group2(pp,dd,ss) - sensitivityGroup2(pp,dd,ss);
            if (any(yneg < 0) | any(ypos < 0))
                error('Check error bar magnitude calculation');
            end
        end
    end
end

%% Make figure 4a
plotSize = [100, 100,400,400];
f = figure('Position',plotSize); clf; hold on

errorbar(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1),yneg(:,1,1),ypos(:,1,1),xneg(:,1,1),xpos(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2),yneg(:,1,2),ypos(:,1,2),xneg(:,1,2),xpos(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1),yneg(:,2,1),ypos(:,2,1),xneg(:,2,1),xpos(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
errorbar(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2),yneg(:,2,2),ypos(:,2,2),xneg(:,2,2),xpos(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);

plot([limMin limMax],[limMin limMax],'k:','LineWidth',2);
set(gca, 'FontName', 'Arial')
legend( ...
    {sprintf('Session %d, %d pixels',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(2)) ; ...
    sprintf('Line of Equality') ...
    }, ...
    'Location','SouthEast');
axis('square');
lgd = legend('show');
lgd.FontSize = 8; 
ax=gca;
xlabel('Group 1 sensitivity (dB)','FontName','Arial','FontWeight','Bold');
ylabel('Group 2 sensitivity (dB)','FontName','Arial','FontWeight','Bold');
ax = gca; % Get current axes handle
ax.XAxis.FontWeight = 'bold'; % Make x-axis tick labels bold
ax.YAxis.FontWeight = 'bold'; % Make y-axis tick labels bold
ticks = 16:2:32; 
xticks(ticks);
yticks(ticks);
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
axis([limMin limMax limMin limMax]);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure4a.png'),'-dpng','-r600');

%% t-tests
[~,p(1,1)] = ttest(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1));
[~,p(1,2)] = ttest(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2));
[~,p(2,1)] = ttest(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1));
[~,p(2,2)] = ttest(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2));
fprintf('Within-Session t-test p values\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.2f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

%% Wilcoxon signed-rank test
[p(1,1),h,~] = signrank(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1));
[p(1,2),h,~] = signrank(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2));
[p(2,1),h,~] = signrank(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1));
[p(2,2),h,~] = signrank(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2));
fprintf('Within-Session p values-Wilcoxon Test\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.2f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

fprintf('\n');
fprintf('Session1, Group 1, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup1(:,1,1)),std(sensitivityGroup1(:,1,1))/sqrt(length(sensitivityGroup1(:,1,1))));
fprintf('Session1, Group2, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup2(:,1,1)),std(sensitivityGroup2(:,1,1))/sqrt(length(sensitivityGroup2(:,1,1))));
fprintf('Session1, Group 1, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup1(:,2,1)),std(sensitivityGroup1(:,2,1))/sqrt(length(sensitivityGroup1(:,2,1))));
fprintf('Session1, Group2, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup2(:,2,1)),std(sensitivityGroup2(:,2,1))/sqrt(length(sensitivityGroup2(:,2,1))));
fprintf('Session2, Group 1, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup1(:,1,2)),std(sensitivityGroup1(:,1,2))/sqrt(length(sensitivityGroup1(:,1,2))));
fprintf('Session2, Group2, 8 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup2(:,1,2)),std(sensitivityGroup2(:,1,2))/sqrt(length(sensitivityGroup2(:,1,2))));
fprintf('Session2, Group 1, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup1(:,2,2)),std(sensitivityGroup1(:,2,2))/sqrt(length(sensitivityGroup1(:,2,2))));
fprintf('Session2, Group2, 43 pixels, mean +/- stderr %0.1f +/- %0.1f\n', mean(sensitivityGroup2(:,2,2)),std(sensitivityGroup2(:,2,2))/sqrt(length(sensitivityGroup2(:,2,2))));

%% Figure 4b (Bland Altman plot), Group1 Vs Group2
% Recall that indices are subject, size (8 and 43), session (1 and 2)
Session1_8pixels_G1 = sensitivityGroup1(:,1,1);
Session2_8pixels_G1 = sensitivityGroup1(:,1,2);
Session1_43pixels_G1 = sensitivityGroup1(:,2,1);
Session2_43pixels_G1 = sensitivityGroup1(:,2,2);
Session1_8pixels_G2 = sensitivityGroup2(:,1,1);
Session2_8pixels_G2 = sensitivityGroup2(:,1,2);
Session1_43pixels_G2 = sensitivityGroup2(:,2,1);
Session2_43pixels_G2 = sensitivityGroup2(:,2,2);
    
[mean_Session1_8pixels, diff_Session1_8pixels] = calculate_bland_altman(Session1_8pixels_G1, Session1_8pixels_G2);
[mean_Session2_8pixels, diff_Session2_8pixels] = calculate_bland_altman(Session2_8pixels_G1, Session2_8pixels_G2);
[mean_Session1_43pixels, diff_Session1_43pixels] = calculate_bland_altman(Session1_43pixels_G1, Session1_43pixels_G2);
[mean_Session2_43pixels, diff_Session2_43pixels] = calculate_bland_altman(Session2_43pixels_G1, Session2_43pixels_G2);

mean_diff_Session1_8pixels = mean(diff_Session1_8pixels);
std_diff_Session1_8pixels = std(diff_Session1_8pixels);
mean_diff_Session2_8pixels = mean(diff_Session2_8pixels);
std_diff_Session2_8pixels = std(diff_Session2_8pixels);
mean_diff_Session1_43pixels = mean(diff_Session1_43pixels);
std_diff_Session1_43pixels = std(diff_Session1_43pixels);
mean_diff_Session2_43pixels = mean(diff_Session2_43pixels);
std_diff_Session2_43pixels = std(diff_Session2_43pixels);

% Plot Bland-Altman data for each comparison
LoA_Session1_8pixels = [mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels];
plotSize = [100 100 200 400];
figure('Position', plotSize);
scatter(mean_Session1_8pixels, diff_Session1_8pixels, 25, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S1');
hold on;
line([18, 23], [mean_diff_Session1_8pixels, mean_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS8, S1)');
line([18, 23], [LoA_Session1_8pixels(2), LoA_Session1_8pixels(2)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8, S1)');
line([18, 23], [LoA_Session1_8pixels(1), LoA_Session1_8pixels(1)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 12)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Mean Sensitivity (dB)'});
ylabel('Group 1 - Group 2 Sensitivity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure4b.png'),'-dpng','-r600');

% Figure 4c
LoA_Session2_8pixels = [mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels,mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels];
figure('Position', plotSize);
scatter(mean_Session2_8pixels, diff_Session2_8pixels, 25, 'blue', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S2');
hold on;
line([18, 23], [mean_diff_Session2_8pixels, mean_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS8, S2)');
line([18, 23], [LoA_Session2_8pixels(2), LoA_Session2_8pixels(2)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS8, S2)');
line([18, 23], [LoA_Session2_8pixels(1), LoA_Session2_8pixels(1)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS8, S2)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 12)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Mean Sensitivity (dB)'});
ylabel('Group 1 - Group 2 Sensitivity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure4c.png'),'-dpng','-r600');

% Figure 4d
LoA_Session1_43pixels = [mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels,mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels];
figure('Position', plotSize);
scatter(mean_Session1_43pixels, diff_Session1_43pixels, 25, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S1');
hold on;
line([24, 32], [mean_diff_Session1_43pixels, mean_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS43, S1)');
line([24, 32], [LoA_Session1_43pixels(2), LoA_Session1_43pixels(2)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43, S1)');
line([24, 32], [LoA_Session1_43pixels(1), LoA_Session1_43pixels(1)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [22, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Mean Sensitivity (dB)'});
ylabel('Group 1 - Group 2 Sensitivity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure4d.png'),'-dpng','-r600');

% Figure 4e
LoA_Session2_43pixels = [mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels,mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels];
figure('Position', plotSize);
scatter(mean_Session2_43pixels, diff_Session2_43pixels, 25, 'red', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S2');
hold on;
line([24, 32], [mean_diff_Session2_43pixels, mean_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS43, S2)');
line([24, 32], [LoA_Session2_43pixels(2), LoA_Session2_43pixels(2)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS43, S2)');
line([24, 32], [LoA_Session2_43pixels(1), LoA_Session2_43pixels(1)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS43, S2)');
LoA_Session2_43pixels = [mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels,mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels];
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [22, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Mean Sensitivity (dB)'});
ylabel('Group 1 - Group 2 Sensitivity (dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure4e.png'),'-dpng','-r600');

%% t-test for 8 Vs 43 pixel stimulus
[~,p(1,1)] = ttest(Session1_8pixels_G1,Session1_43pixels_G1);
[~,p(1,2)] = ttest(Session1_8pixels_G2,Session1_43pixels_G2);
[~,p(2,1)] = ttest(Session2_8pixels_G1,Session2_43pixels_G1);
[~,p(2,2)] = ttest(Session2_8pixels_G2,Session2_43pixels_G2);
theSplits= [1,2]; 
fprintf('\nWithin-Session (8 Vs 43 pixels) t-test p values\n');
for dd = 1:length(theSplits)
    for ss = 1:length(theSessions)
        fprintf('\t Group %d, session %d, p = %0.2f\n',theSplits(dd),theSessions(ss),p(dd,ss));
    end
end

%% Print limits of agreement
fprintf('\n');
fprintf('Group1 Vs Group2 : Session1, 8 pixels, LoA: %.2f, %.2f\n', LoA_Session1_8pixels);
fprintf('Group1 Vs Group2 : Session1, 43 pixels, LoA: %.2f, %.2f\n', LoA_Session1_43pixels);
fprintf('Group1 Vs Group2 : Session2, 8 pixels, LoA: %.2f, %.2f\n', LoA_Session2_8pixels);
fprintf('Group1 Vs Group2 : Session2, 43 pixels, LoA: %.2f, %.2f\n', LoA_Session2_43pixels);

