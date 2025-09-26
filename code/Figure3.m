%% Make Figure 3 for the paper
%
% TODO
%    a) Make 3a and 3b match in terms of size, fonts, etc.

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
        end
    end
end

%% Make Figure 3a
limMin = 18; limMax = 30;
f = figure('Position',plotSize); clf; hold on
set(gca,'FontName', 'Helvetica','FontSize',axisFontSize);
plot(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);
plot([limMin limMax],[limMin limMax],'k:','LineWidth',1);
xlabel('MOCS sensitivity (dB)','FontSize',labelFontSize);
ylabel('QUEST sensitivity (dB)','FontSize',labelFontSize);
legend( ...
    {sprintf('Session %d, %d pixels',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(2)) ; ...
    ''}, ...
    'Location','SouthEast');
axis('square');
axis([limMin limMax limMin limMax]);
saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3a.pdf'),'pdf');

%% t-tests
[~,p(1,1)] = ttest(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1));
[~,p(1,2)] = ttest(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2));
[~,p(2,1)] = ttest(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1));
[~,p(2,2)] = ttest(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2));
fprintf('MOCS vs QUEST t-test p values\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

%% Get all Session 1 and Session 2 data for combined trials
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'COMBINED') & strcmp(dataTable.Split,'All'));
            if (sum(index) ~= 1)
                error('Have not properly set up condition to pick out just one sensitivity');
            end
            sensitivitySessionwise(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
        end
    end
end

%% Figure 3b (Bland Altman plot) %% MOCS VS QUEST
%
% Recall that indices are subject, size (8 and 43), session (1 and 2)
Session1_8pixels_M = sensitivityMOCS(:,1,1);
Session2_8pixels_M = sensitivityMOCS(:,1,2);
Session1_43pixels_M = sensitivityMOCS(:,2,1);
Session2_43pixels_M = sensitivityMOCS(:,2,2);
Session1_8pixels_Q = sensitivityQUEST(:,1,1);
Session2_8pixels_Q = sensitivityQUEST(:,1,2);
Session1_43pixels_Q = sensitivityQUEST(:,2,1);
Session2_43pixels_Q = sensitivityQUEST(:,2,2);
    
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
ax.XLim = [16, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean (MOCS,QUEST) Sensitivity in dB', 'FontWeight', 'bold', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('(MOCS-QUEST) Sensitivity', 'FontWeight', 'bold', 'FontSize', 18, 'FontName', 'Times New Roman');




