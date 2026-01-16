% FitPsychometricFunction
%
% Massage data, fit and plot psychometric function
function [plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues,h] = FitPsychometricFunction(log_intensity,responses, ...
    log0Value,thresholdCriterion,convert_to_db,make_plot,figure_vis, ...
    guessUpper,lapseUpper,slopeLower,slopeUpper)

% Sort data and get what we need
sorted_data = sortrows([log_intensity responses],1);
sorted_log_intensities = sorted_data(:,1);
index = find(isinf(sorted_log_intensities));
if (~isempty(index))
    sorted_log_intensities(isinf(sorted_log_intensities)) = log0Value;
end
sorted_responses = sorted_data(:,2);

% Get unique stimulus levels and data for those
[unique_log_intensities] = unique(sorted_log_intensities,'legacy');
num_presented = zeros(size(unique_log_intensities));
for i = 1:length(num_presented)
    num_presented(i) = length(find(sorted_log_intensities == unique_log_intensities(i)));
    num_seen(i) = sum(sorted_responses(find(sorted_log_intensities == unique_log_intensities(i))));
end

% Fit using a wrapper into Palemedes
if (make_plot)
    fprintf('\t\tCalling Palemedes wrapper ... ');
end
[plot_log_intensities,plot_psychometric,corrected_psychometric, ...
    log_threshold,corrected_log_threshold, psiParamsValues] = ...
    fitPsychometric(unique_log_intensities,num_seen',num_presented,thresholdCriterion,convert_to_db,guessUpper,lapseUpper,slopeLower,slopeUpper);
if (make_plot)
    fprintf('done\n');
end

% Plot
if (make_plot)
    if (convert_to_db)
        maxIntensityLog= 1;
        minIntensityLog=-36;
    else
        maxIntensityLog= 0.1;
        minIntensityLog=-3.6;
    end
    
    fontsize = 14; fwidth = 3.5; fheight = 3.5;
    h = figure('Units', 'inches','Position', [400 200 fwidth fheight],'visible',figure_vis); hold on;
    a0 = gca;
    ax.LineWidth = 2;ax.LineWidth = 2;
    set(a0,'FontName','Arial','FontSize',14);
    if (convert_to_db)
        xlabel('Trial intensity (dB)','FontSize',fontsize);
    else
        xlabel('Log trial intensity (au)','FontSize',fontsize);
    end
    ylabel('Fraction seen','FontSize',fontsize);
    xlim([minIntensityLog maxIntensityLog]);
    ylim([0 1]);
    set(a0,'FontSize',fontsize);
    set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(a0,'Color','none');
    set(h,'Color',[1 1 1]);
    set(h,'PaperPositionMode','auto');
    set(h, 'renderer', 'painters');
    baseMarkerSize = 5;
    sizeNormalizer = 40;
    for tt = 1:length(unique_log_intensities)
        markerSize = 4+round(baseMarkerSize*num_presented(tt)/sizeNormalizer);
        plot(unique_log_intensities(tt), num_seen(tt)' ./ num_presented(tt)','ro','MarkerSize',markerSize,'MarkerFaceColor','r');
    end

    plot(plot_log_intensities,corrected_psychometric,'k-','LineWidth',3,'HandleVisibility', 'off');
    plot(plot_log_intensities,plot_psychometric,'r:','LineWidth',2,'HandleVisibility', 'off');
    plot([minIntensityLog corrected_log_threshold],[thresholdCriterion thresholdCriterion],'k:','LineWidth',2,'HandleVisibility', 'off');
    plot([corrected_log_threshold corrected_log_threshold],[0 thresholdCriterion],'k:','LineWidth',2,'HandleVisibility', 'off');
    
    
[minVal, ~] = min(num_presented);
[maxVal, ~] = max(num_presented);

% Calculate the marker sizes 
szMin = 4+round(baseMarkerSize*minVal/sizeNormalizer);
szMax = 4+round(baseMarkerSize*maxVal/sizeNormalizer);

% Create dummy points for legend 
% These are invisible on the plot but appear in the legend
hMin = plot(NaN, NaN, 'ro', 'MarkerSize', szMin, 'MarkerFaceColor', 'r', ...
    'DisplayName', [num2str(minVal) ' Trials']);

hMax = plot(NaN, NaN, 'ro', 'MarkerSize', szMax, 'MarkerFaceColor', 'r', ...
    'DisplayName', [num2str(maxVal) ' Trials']);

% Define pluralization for the legend 
if minVal == 1
    minTrialSuffix = 'Trial';
else
    minTrialSuffix = 'Trials';
end

% Build the legend only with trial markers

if minVal == maxVal
    % Only one handle (hMin) and one label
    hArray = [hMin];
    labels = {[num2str(minVal) ' ' minTrialSuffix]};
else
    % Two handles (hMin, hMax) and two labels
    hArray = [hMin, hMax];
    labels = {[num2str(minVal) ' ' minTrialSuffix], [num2str(maxVal) ' Trials']};
end

% Legend 
% This will now only show the red circle(s) with trial counts
lgd = legend(hArray, labels, 'Location', 'southeast', 'FontSize',12);
legend boxoff;
else
    h  = [];
end

end

