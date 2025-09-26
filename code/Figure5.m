%% Session 1 VS Session 2 plots

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

%% Figure 5b (Bland Altman plot) 
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

% Open the figures
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