
%% Clear
clear; close all;

%% Output variant
outputVariant = 'SlopeFree1';

%% Parameters
axisFontSize = 14;
labelFontSize = 18;
titleFontSize = 20;
legendFontSize = 14;
markerSize = 12;

%% Path to analysis output dir
analysisDir = getpref('AOMicroRepeat','analysisDir');

%% Read table of fit information
wState = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet'); %,'VariableNamingRule','preserve');
warning(wState);

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
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
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
f = figure; clf; hold on
set(gca,'FontName', 'Helvetica','FontSize',axisFontSize);
plot(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);
plot([limMin limMax],[limMin limMax],'k:','LineWidth',1);
xlabel('MOCS sensitivity (dB)','FontSize',labelFontSize);
ylabel('QUEST sensitivity (dB)','FontSize',labelFontSize);
legend( ...
    {sprintf('Session %d, %d arcmin',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(2),theDiameters(2)) ; ...
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
        fprintf('\t%d arcmin, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end
%% Get all Session 1 and Session 2 data for combined trials
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
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

%% Figure 5a
f = figure; clf; hold on
set(gca,'FontName', 'Helvetica','FontSize',axisFontSize);
plot(sensitivitySessionwise(:,1,1),sensitivitySessionwise(:,1,2),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivitySessionwise(:,2,1),sensitivitySessionwise(:,2,2),'ro','MarkerFaceColor','r','MarkerSize',markerSize);

plot([limMin limMax],[limMin limMax],'k:','LineWidth',1);
xlabel('Session 1 Sensitivity (dB)','FontSize',labelFontSize);
ylabel('Session 2 Sensitivity (dB)','FontSize',labelFontSize);
legend( ...
    {sprintf('%d arcmin',theDiameters(1)) ; ...
    sprintf('%d arcmin',theDiameters(2)) ; ...
    sprintf('Line of Equality')
    ''}, ...
    'Location','SouthEast');
axis('square');
lgd = legend('show');
lgd.FontSize = 10; 
axis([limMin limMax limMin limMax]);
saveas(gcf,fullfile(analysisDir,outputVariant,'Figure4a.pdf'),'pdf');
 %% t-tests for between session comparision
[~,p(1,1)] = ttest(sensitivitySessionwise(:,1,1),sensitivitySessionwise(:,1,2));
[~,p(1,2)] = ttest(sensitivitySessionwise(:,2,1),sensitivitySessionwise(:,2,2));

fprintf('Session 1 vs Session 2 t-test p values\n');
for dd = 1:length(theDiameters)
        fprintf('\t%d arcmin, p = %0.3f\n',theDiameters(dd),p(dd,ss));
end
%% Figure 3B (Bland Altman plot) %% MOCS VS QUEST
Session1_8arcmin_M = sensitivityMOCS(:,1,1);
Session2_8arcmin_M = sensitivityMOCS(:,1,2);
Session1_43arcmin_M = sensitivityMOCS(:,2,1);
Session2_43arcmin_M = sensitivityMOCS(:,2,2);
Session1_8arcmin_Q = sensitivityQUEST(:,1,1);
Session2_8arcmin_Q = sensitivityQUEST(:,1,2);
Session1_43arcmin_Q = sensitivityQUEST(:,2,1);
Session2_43arcmin_Q = sensitivityQUEST(:,2,2);
    
    
[mean_a, diff_a] = calculate_bland_altman(Session1_8arcmin_M, Session1_8arcmin_Q);
[mean_b, diff_b] = calculate_bland_altman(Session2_8arcmin_M, Session2_8arcmin_Q);
[mean_c, diff_c] = calculate_bland_altman(Session1_43arcmin_M, Session1_43arcmin_Q);
[mean_d, diff_d] = calculate_bland_altman(Session2_43arcmin_M, Session2_43arcmin_Q);

mean_diff_S1_SS8 = mean(diff_a);
std_diff_S1_SS8 = std(diff_a);
mean_diff_S2_SS8 = mean(diff_b);
std_diff_S2_SS8 = std(diff_b);
mean_diff_S1_SS43 = mean(diff_c);
std_diff_S1_SS43 = std(diff_c);
mean_diff_S2_SS43 = mean(diff_d);
std_diff_S2_SS43 = std(diff_d);

figure('Position', [100, 100,800,800]);

% Plot Bland-Altman data for each comparison
scatter(mean_a, diff_a, 150, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S1');
hold on;
scatter(mean_b, diff_b, 150, 'blue', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S2');
scatter(mean_c, diff_c, 150, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S1');
scatter(mean_d, diff_d, 150, 'red', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S2');

% Plot mean differences and limits of agreement
line([18, 23], [mean_diff_S1_SS8, mean_diff_S1_SS8], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS8, S1)');
line([18, 23], [mean_diff_S2_SS8, mean_diff_S2_SS8], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS8, S2)');
line([26, 30], [mean_diff_S1_SS43, mean_diff_S1_SS43], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS43, S1)');
line([26, 30], [mean_diff_S2_SS43, mean_diff_S2_SS43], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS43, S2)');

line([18, 23], [mean_diff_S1_SS8 + 1.96 * std_diff_S1_SS8, mean_diff_S1_SS8 + 1.96 * std_diff_S1_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8, S1)');
line([18, 23], [mean_diff_S1_SS8 - 1.96 * std_diff_S1_SS8, mean_diff_S1_SS8 - 1.96 * std_diff_S1_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8, S1)');
line([18, 23], [mean_diff_S2_SS8 + 1.96 * std_diff_S2_SS8, mean_diff_S2_SS8 + 1.96 * std_diff_S2_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS8, S2)');
line([18, 23], [mean_diff_S2_SS8 - 1.96 * std_diff_S2_SS8, mean_diff_S2_SS8 - 1.96 * std_diff_S2_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS8, S2)');

line([26, 30], [mean_diff_S1_SS43 + 1.96 * std_diff_S1_SS43, mean_diff_S1_SS43 + 1.96 * std_diff_S1_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43, S1)');
line([26, 30], [mean_diff_S1_SS43 - 1.96 * std_diff_S1_SS43, mean_diff_S1_SS43 - 1.96 * std_diff_S1_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43, S1)');
line([26, 30], [mean_diff_S2_SS43 + 1.96 * std_diff_S2_SS43, mean_diff_S2_SS43 + 1.96 * std_diff_S2_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS43, S2)');
line([26, 30], [mean_diff_S2_SS43 - 1.96 * std_diff_S2_SS43, mean_diff_S2_SS43 - 1.96 * std_diff_S2_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS43, S2)');

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

%% Figure 5B (Bland Altman plot) %% Session 1 VS Session 2

Session1_8arcmin = sensitivitySessionwise(:,1,1);

Session2_8arcmin = sensitivitySessionwise(:,1,2);

Session1_43arcmin = sensitivitySessionwise(:,2,1);

Session2_43arcmin = sensitivitySessionwise(:,2,2);


[mean_a, diff_a] = calculate_bland_altman(Session1_8arcmin, Session2_8arcmin);
[mean_b, diff_b] = calculate_bland_altman(Session1_43arcmin, Session2_43arcmin);


mean_diff_SS8 = mean(diff_a);
std_diff_SS8 = std(diff_a);

mean_diff_SS43 = mean(diff_b);
std_diff_SS43 = std(diff_b);

figure('Position', [100, 100,800,800]);

% Plot Bland-Altman data for each comparison
scatter(mean_a, diff_a, 150, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8');
hold on;
scatter(mean_b, diff_b, 150, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43');

% Plot mean differences and limits of agreement
line([18, 23], [mean_diff_SS8, mean_diff_SS8], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference SS8');
line([26, 30], [mean_diff_SS43, mean_diff_SS43], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference SS43');


line([18, 23], [mean_diff_SS8 + 1.96 * std_diff_SS8, mean_diff_SS8 + 1.96 * std_diff_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8)');
line([18, 23], [mean_diff_SS8 - 1.96 * std_diff_SS8, mean_diff_SS8 - 1.96 * std_diff_SS8], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8)');

line([26, 30], [mean_diff_SS43 + 1.96 * std_diff_SS43, mean_diff_SS43 + 1.96 * std_diff_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43)');
line([26, 30], [mean_diff_SS43 - 1.96 * std_diff_SS43, mean_diff_SS43 - 1.96 * std_diff_SS43], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43)');


ax = gca;
ax.XLim = [16, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.PlotBoxAspectRatio = [1 1 1]; 

xlabel('Mean (Session 1 & 2) Sensitivity in dB', 'FontWeight', 'bold', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('(Session 1 - Session 2) Sensitivity', 'FontWeight', 'bold', 'FontSize', 18, 'FontName', 'Times New Roman');

%% Function to calculate mean and differences
function [mean_val, diff_val] = calculate_bland_altman(x1, x2)
    mean_val = (x1 + x2) / 2;
    diff_val = x1 - x2;
end
