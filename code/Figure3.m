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
            CILower_MOCS(pp,dd,ss)= -dataTable.CIHigh(index); %negative sensitivity, so high value corresonds to lower CI)
            CIUpper_MOCS(pp,dd,ss)= -dataTable.CILow(index);%negative sensitivity, so low value corresonds to higher CI)
            xneg(pp,dd,ss) = sensitivityMOCS(pp,dd,ss) - CIUpper_MOCS(pp,dd,ss);
            xpos(pp,dd,ss) = sensitivityMOCS(pp,dd,ss)- CILower_MOCS(pp,dd,ss);
            
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
            CILower_QUEST(pp,dd,ss)= -dataTable.CIHigh(index); %negative sensitivity, so high value corresonds to lower CI
            CIUpper_QUEST(pp,dd,ss)= -dataTable.CILow(index);%negative sensitivity, so low value corresonds to higher CI
            yneg(pp,dd,ss) = sensitivityQUEST(pp,dd,ss) - CIUpper_QUEST(pp,dd,ss);
            ypos(pp,dd,ss) = sensitivityQUEST(pp,dd,ss)- CILower_QUEST(pp,dd,ss);
          
        end
    end
end

%% Make Figure 3a

f = figure('Position',plotSize); clf; hold on

errorbar(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1),yneg(:,1,1),ypos(:,1,1),xneg(:,1,1),xpos(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2),yneg(:,1,2),ypos(:,1,2),xneg(:,1,2),xpos(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1),yneg(:,2,1),ypos(:,2,1),xneg(:,2,1),xpos(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
errorbar(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2),yneg(:,2,2),ypos(:,2,2),xneg(:,2,2),xpos(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);

plot([limMin limMax],[limMin limMax],'k:','LineWidth',2);
xlabel('MOCS sensitivity (dB)', 'FontWeight','bold', 'FontSize', 18, 'FontName', 'Times New Roman');
ylabel('QUEST sensitivity (dB)','FontWeight','bold', 'FontSize', 18, 'FontName', 'Times New Roman');
legend( ...
    {sprintf('Session %d, %d pixels',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d pixels',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d pixels',theSessions(2),theDiameters(2)) ; ...
    sprintf('Line of Equality')
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
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3a.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figure3a.png', '-dpng', '-r600');
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
%% Wilcoxon test% Wilcoxon signed-rank test
[p(1,1),h,~] = signrank(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1));
[p(1,2),h,~] = signrank(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2));
[p(2,1),h,~] = signrank(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1));
[p(2,2),h,~] = signrank(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2));
% Results
fprintf('MOCS vs QUEST t-test p values-Wilcoxon Test\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d pixels, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end

%% Figure 3b (Bland Altman plot) %% MOCS VS QUEST
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


%% Plotting seperately 
%% 

%%(3b) Session1 8 Pixels 
figure('Position', plotSize);
scatter(mean_Session1_8pixels, diff_Session1_8pixels, 150, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S1');% Plot Bland-Altman data for each comparison
line([18, 23], [mean_diff_Session1_8pixels, mean_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS8, S1)');% Plot mean differences and limits of agreement
line([18, 23], [mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels + 1.96 * std_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8, S1)');
line([18, 23], [mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels, mean_diff_Session1_8pixels - 1.96 * std_diff_Session1_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean of MOCS and QUEST Sensitivity (dB)');
ylabel('Difference of MOCS and QUEST Sensitivities (dB)');
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3b.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figure3b.png', '-dpng', '-r600');
%% (3c) Session 2 8 pixels
figure('Position', plotSize);
scatter(mean_Session2_8pixels, diff_Session2_8pixels, 150, 'blue', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8S2');
line([18, 23], [mean_diff_Session2_8pixels, mean_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS8, S2)');
line([18, 23], [mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels, mean_diff_Session2_8pixels + 1.96 * std_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS8, S2)');
line([18, 23], [mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels, mean_diff_Session2_8pixels - 1.96 * std_diff_Session2_8pixels], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS8, S2)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean of MOCS and QUEST Sensitivity (dB)');
ylabel('Difference of MOCS and QUEST Sensitivities(dB)');
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3c.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figure3c.png', '-dpng', '-r600');
%% (3d) Session 1 43 pixels
figure('Position', plotSize);
scatter(mean_Session1_43pixels, diff_Session1_43pixels, 150, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S1');
line([26, 30], [mean_diff_Session1_43pixels, mean_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference (SS43, S1)');
line([26, 30], [mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels, mean_diff_Session1_43pixels - 1.96 * std_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43, S1)');
line([26, 30], [mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels, mean_diff_Session1_43pixels + 1.96 * std_diff_Session1_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43, S1)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [24, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean of MOCS and QUEST Sensitivities (dB)');
ylabel('Difference of MOCS and QUEST Sensitivities(dB)');
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3d.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figure3d.png', '-dpng', '-r600');
%% (3e) Session 2 43 pixels
figure('Position', plotSize);
scatter(mean_Session2_43pixels, diff_Session2_43pixels, 150, 'red', 'filled', '^', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43S2');
line([26, 30], [mean_diff_Session2_43pixels, mean_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '--', 'DisplayName', 'Mean Difference (SS43, S2)');
line([26, 30], [mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels, mean_diff_Session2_43pixels + 1.96 * std_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Upper Limit (SS43, S2)');
line([26, 30], [mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels, mean_diff_Session2_43pixels - 1.96 * std_diff_Session2_43pixels], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--',  'DisplayName', 'Lower Limit (SS43, S2)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 18)
ax.XLim = [24, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 18;
ax.PlotBoxAspectRatio = [1 1 1]; % Maintain aspect ratio

xlabel('Mean of MOCS and QUEST Sensitivity (dB)');
ylabel('Difference of MOCS and QUEST Sensitivities (dB)');
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3e.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf, 'figure3e.png', '-dpng', '-r600');



