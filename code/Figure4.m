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
            sensitivityGroup1(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
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
            sensitivityGroup2(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
        end
    end
end
%% Make figure 4a
f = figure('Position',plotSize); clf; hold on

plot(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
plot(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);
plot([limMin limMax],[limMin limMax],'k:','LineWidth',2);
xlabel('Group1 sensitivity (dB)', 'FontWeight','bold', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('Group2 sensitivity (dB)','FontWeight','bold', 'FontSize', 18, 'FontName', 'Times New Roman');
legend( ...
    {sprintf('Session %d, %d pixels',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(2)) ; ...
    sprintf('Line of Equality')
    ''}, ...
    'Location','SouthEast');
axis('square');

lgd = legend('show');
lgd.FontSize = 12; 
ax=gca;
set(gca, 'FontName', 'Arial','FontWeight','bold')
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
axis([limMin limMax limMin limMax]);
saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3a.pdf'),'pdf');


%% %% t-tests
[~,p(1,1)] = ttest(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1));
[~,p(1,2)] = ttest(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2));
[~,p(2,1)] = ttest(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1));
[~,p(2,2)] = ttest(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2));
fprintf('MOCS vs QUEST t-test p values\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end
%% %% Wilcoxon test% Wilcoxon signed-rank test
[p(1,1),h,~] = signrank(sensitivityGroup1(:,1,1),sensitivityGroup2(:,1,1));
[p(1,2),h,~] = signrank(sensitivityGroup1(:,1,2),sensitivityGroup2(:,1,2));
[p(2,1),h,~] = signrank(sensitivityGroup1(:,2,1),sensitivityGroup2(:,2,1));
[p(2,2),h,~] = signrank(sensitivityGroup1(:,2,2),sensitivityGroup2(:,2,2));
% Results
fprintf('MOCS vs QUEST t-test p values-Wilcoxon Test\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end
%% Figure 4b (Bland Altman plot) %% Group1 Vs Group2
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

figure('Position', plotSize);

% Plot Bland-Altman data for each comparison
scatter(mean_Session1_8pixels, diff_Session1_8pixels, 150, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S1');
hold on;
scatter(mean_Session2_8pixels, diff_Session2_8pixels, 150, 'blue', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S2');
scatter(mean_Session1_43pixels, diff_Session1_43pixels, 150, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S1');
scatter(mean_Session2_43pixels, diff_Session2_43pixels, 150, 'red', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S2');

% Plot mean differences and limits of agreement
line([18, 23], [mean_diff_Session1_8pixels, mean_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS8, S1)');
line([18, 23], [mean_diff_Session2_8pixels, mean_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS8, S2)');
line([26, 30], [mean_diff_Session1_43pixels, mean_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS43, S1)');
line([26, 30], [mean_diff_Session2_43pixels, mean_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS43, S2)');

line([18, 23], [mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8, S1)');
line([18, 23], [mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8, S1)');
line([18, 23], [mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels, mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS8, S2)');
line([18, 23], [mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels, mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS8, S2)');

line([26, 30], [mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels, mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43, S1)');
line([26, 30], [mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels, mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43, S1)');
line([26, 30], [mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels, mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS43, S2)');
line([26, 30], [mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels, mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS43, S2)');

ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [16, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean of Group1 and Group2 Sensitivity (dB)');
ylabel('(Group1-Group2)Sensitivity (dB)');


