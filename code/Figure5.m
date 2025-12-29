%% Session 1 VS Session 2 plots
%% Clear
clear; close all;

%% Output variant
outputVariant = 'SlopeFree1';

%% Set up for the figures and read data
FigureSetup;

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

            % See comments in Figure3.m about this.  David's fix implemented here.
            CILower_SessionWise(pp,dd,ss)= -dataTable.CIHigh(index);
            CIUpper_SessionWise(pp,dd,ss)= -dataTable.CILow(index);
            neg(pp,dd,ss) = sensitivitySessionwise(pp,dd,ss) - CILower_SessionWise(pp,dd,ss);
            pos(pp,dd,ss) = CIUpper_SessionWise(pp,dd,ss) - sensitivitySessionwise(pp,dd,ss);
        end
    end
end

%% Figure 5a
plotSize = [100 100 400 400];
f = figure('Position',plotSize); clf; hold on
errorbar(sensitivitySessionwise(:,1,1),sensitivitySessionwise(:,1,2),neg(:,1,1),pos(:,1,1),neg(:,1,2),pos(:,1,2),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
errorbar(sensitivitySessionwise(:,2,1),sensitivitySessionwise(:,2,2),neg(:,2,1),pos(:,2,1),neg(:,2,2),pos(:,2,2),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
plot([limMin limMax],[limMin limMax],'k:','LineWidth',2);
xlabel('Session 1 Sensitivity (dB)');
ylabel('Session 2 Sensitivity (dB)');
legend( ...
    {sprintf('%d pixels',theDiameters(1)) ; ...
    sprintf('%d pixels',theDiameters(2)) ; ...
    sprintf('Line of Equality')}, ...
    'Location','SouthEast');
axis('square');
ax=gca;
set(gca, 'FontName', 'Arial','FontWeight','bold')
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ticks = 16:2:32; 
xticks(ticks);
yticks(ticks);
lgd = legend('show');
lgd.FontSize = 12; 
axis([limMin limMax limMin limMax]);
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure5a.png'),'-dpng','-r600');

 %% t-tests for between session comparision
[~,p(1,1)] = ttest(sensitivitySessionwise(:,1,1),sensitivitySessionwise(:,1,2));
[~,p(1,2)] = ttest(sensitivitySessionwise(:,2,1),sensitivitySessionwise(:,2,2));
fprintf('Session 1 vs Session 2 t-test p values\n');
for dd = 1:length(theDiameters)
    fprintf('\t%d pixels, p = %0.2f\n',theDiameters(dd),p(dd));
end

%% Wilcoxon signed-rank test
[p(1,1),h,~] = signrank(sensitivitySessionwise(:,1,1),sensitivitySessionwise(:,1,2));
[p(1,2),h,~] = signrank(sensitivitySessionwise(:,2,1),sensitivitySessionwise(:,2,2));
fprintf('Session 1 vs Session 2 t-test p values (Wilcoxon-test)\n');
for dd = 1:length(theDiameters)
    fprintf('\t%d pixels, p = %0.2f\n',theDiameters(dd),p(dd));
end

%% Figure 5b (Bland Altman plot) 
Session1_8pixels = sensitivitySessionwise(:,1,1);
Session2_8pixels = sensitivitySessionwise(:,1,2);
Session1_43pixels = sensitivitySessionwise(:,2,1);
Session2_43pixels = sensitivitySessionwise(:,2,2);

[mean_SS8, diff_SS8] = calculate_bland_altman(Session1_8pixels, Session2_8pixels);
[mean_SS43, diff_SS43] = calculate_bland_altman(Session1_43pixels, Session2_43pixels);

mean_diff_SS8 = mean(diff_SS8);
std_diff_SS8 = std(diff_SS8);
mean_diff_SS43 = mean(diff_SS43);
std_diff_SS43 = std(diff_SS43);

% Figure 5b
plotSize = [100 100 200 400];
figure('Position', plotSize);

% Plot Bland-Altman data for each comparison
LOA_8pixels = [mean_diff_SS8 - 1.96 * std_diff_SS8,mean_diff_SS8 + 1.96 * std_diff_SS8];
scatter(mean_SS8, diff_SS8, 25, 'blue', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS8');
hold on;
line([18, 23], [mean_diff_SS8, mean_diff_SS8], 'Color', 'blue', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference SS8');
line([18, 23], [LOA_8pixels(2), LOA_8pixels(2)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Upper Limit (SS8)');
line([18, 23], [LOA_8pixels(1), LOA_8pixels(1)], 'Color', 'blue', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS8)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 12)
ax.XLim = [16, 24];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel({'Mean Sensitivity (dB)'});
ylabel('Session 1 - Session 2 Sensitivity(dB)');
% saveas(gcf,fullfile(analysisDir,outputVariant,'Figure5b.pdf'),'pdf');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure5b.png'),'-dpng','-r600');

% Figure 5c
LOA_43pixels = [mean_diff_SS43 - 1.96 * std_diff_SS43,mean_diff_SS43 + 1.96 * std_diff_SS43];
figure('Position', plotSize);
scatter(mean_SS43, diff_SS43, 25, 'red', 'filled', 'o', 'MarkerFaceAlpha', 0.6, 'DisplayName', 'SS43');
hold on;
line([26, 30], [mean_diff_SS43, mean_diff_SS43], 'Color', 'red', 'LineWidth', 4, 'LineStyle', '-', 'DisplayName', 'Mean Difference SS43');
line([26, 30], [LOA_43pixels(2) LOA_43pixels(2)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Lower Limit (SS43)');
line([26, 30], [LOA_43pixels(1) LOA_43pixels(1)], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '-',  'DisplayName', 'Upper Limit (SS43)');
ax = gca;
set(gca, 'FontName', 'Arial','FontWeight','bold','FontSize', 12)
ax.XLim = [24, 32];
ax.YLim = [-2.5, 2.5];
ax.XAxis.LineWidth = 2;
ax.YAxis.LineWidth = 2;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ticks = 24:2:32; 
xticks(ticks);
xlabel({'Mean Sensitivity (dB)'});
ylabel('Session 1 - Session 2 Sensitivity(dB)');
set(gcf, 'PaperPositionMode', 'auto');
print(gcf,fullfile(analysisDir,outputVariant,'figure5c.png'),'-dpng','-r600');

%% t-test for 8 Vs 43 pixel stimulus
[~,p(1,1)] = ttest(Session1_8pixels,Session1_43pixels);
[~,p(1,2)] = ttest(Session2_8pixels,Session2_43pixels);
fprintf('Between-Session (8 Vs 43 pixels) t-test p values\n');
for ss = 1:length(theSessions)
    fprintf('\t Session %d, p = %0.4f\n',theSessions(ss),p(ss));
end

%% %% Print limits of agreement
fprintf('\n');
fprintf('Session 1 Vs Session 2 : 8 pixels, LoA: %.1f, %.1f\n', LOA_8pixels);
fprintf('Session 1 Vs Session 2 : 43 pixels, LoA: %.1f, %.1f\n', LOA_43pixels);
fprintf('Session 1 Vs Session 2 : 8 pixels, CoR: %.1f\n', 1.96 * std_diff_SS8);
fprintf('Session 1 Vs Session 2 : 43 pixels, CoR: %.1f\n', 1.96 * std_diff_SS43);