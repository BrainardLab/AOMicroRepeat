% Make some plots of the computational observer results
%
% Description:
%   Run out the computational observer for the two spot experiment with a
%   variety of defocus and pupil sizes.

% Clear
clear; close all;

% Parameters
defocusDioptersList = [0 0.05 0.10 0.15];
pupilDiameterMmList = [2 3 4 5 6 7];
testing = false;

% Directory stuff
baseProject = 'AOCompObserver';
analysisBaseDir = getpref(baseProject,'analysisDir');

% Do them
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        pupilDiameterMm = pupilDiameterMmList(pp);
        defocusDiopters = defocusDioptersList(dd);

        analysisOutDir = fullfile(analysisBaseDir,sprintf('IncrDecr1_%s_%d',num2str(round(1000*defocusDiopters)),pupilDiameterMm));
        theData{dd,pp} = load(fullfile(analysisOutDir,'CompObserver'));
    end
end

%% Plot as a function of pupil size
%
% Choose one defocus
plotDefocusDiopters = 0;

% Set up and plot
theColors = ['r' 'g' 'k' 'b' 'c' 'm'];
colorIndex = 1;
thePupilSizeFig = figure; clf; hold on
defocusIndex = find(defocusDioptersList == defocusDiopters);
legendText = {};
for pp = 1:length(pupilDiameterMmList)
    theEllipse = theData{defocusIndex,pp}.fitEllipse;
    theEllipseNorm = vecnorm(theEllipse);
    if (pp == 1)
        theBaseEllipse = theEllipse;
        theBaseEllipseNorm = vecnorm(theBaseEllipse);
    end
    subplot(1,2,1); hold on
    plot(theEllipse(1,:),theEllipse(2,:),theColors(colorIndex),'LineWidth',3);
    
    % Scaled ellipse
    scaleEllipse = (theEllipseNorm'\theBaseEllipseNorm')*theEllipse;
    subplot(1,2,2); hold on
    plot(scaleEllipse(1,:),scaleEllipse(2,:),theColors(colorIndex),'LineWidth',3);
    
    % Bump color index
    colorIndex = colorIndex+1;
    if (colorIndex > length(theColors))
        colorIndex = 1;
    end
    
    % Legend text
    legendText{pp} = sprintf('%d mm',pupilDiameterMmList(pp));
    
end

% Plot finish up up
set(thePupilSizeFig, 'Position',  [800, 0, 1000, 600]);
contrastLim = 0.10;
subplot(1,2,1);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('Contrast 1');
ylabel('Contrsast 2');
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');
title(sprintf('Defocus %d D, raw scale',plotDefocusDiopters));
legend(legendText,'Location','NorthEast');

subplot(1,2,2);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('Contrast 1');
ylabel('Contrsast 2');
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');
title(sprintf('Defocus %d D, same scale',plotDefocusDiopters));
legend(legendText,'Location','NorthEast');

%% Plot as a function of diopters
%
% Choose one defocus
plotPupilSizeMm = 7;

% Set up and plot
theColors = ['r' 'g' 'k' 'b' 'c' 'm'];
colorIndex = 1;
legendIndex = 1;
theDioptersFig = figure; clf; hold on
pupilDiameterIndex = find(pupilDiameterMmList == plotPupilSizeMm );
legendText = {};
for dd = length(defocusDioptersList):-1:1
    theEllipse = theData{dd,pupilDiameterIndex}.fitEllipse;
    theEllipseNorm = vecnorm(theEllipse);
    if (dd == length(defocusDioptersList))
        theBaseEllipse = theEllipse;
        theBaseEllipseNorm = vecnorm(theBaseEllipse);
    end
    subplot(1,2,1); hold on
    plot(theEllipse(1,:),theEllipse(2,:),theColors(colorIndex),'LineWidth',3);
    
    % Scaled ellipse
    scaleEllipse = (theEllipseNorm'\theBaseEllipseNorm')*theEllipse;
    subplot(1,2,2); hold on
    plot(scaleEllipse(1,:),scaleEllipse(2,:),theColors(colorIndex),'LineWidth',3);
    
    % Bump color index
    colorIndex = colorIndex+1;
    if (colorIndex > length(theColors))
        colorIndex = 1;
    end
    
    % Add to legend text
    legendText{legendIndex} = sprintf('%0.2f D',defocusDioptersList(dd));
    legendIndex = legendIndex+1;
    
end

% Plot finish up up
set(theDioptersFig, 'Position',  [800, 0, 1000, 600]);
contrastLim = 0.02;
subplot(1,2,1);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('Contrast 1');
ylabel('Contrsast 2');
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');
title(sprintf('Pupil diameter %d mm, raw scale',plotPupilSizeMm));
legend(legendText,'Location','NorthEast');

subplot(1,2,2);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('Contrast 1');
ylabel('Contrsast 2');
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');
title(sprintf('Pupil diameter %d mm, normalized',plotPupilSizeMm));
legend(legendText,'Location','NorthEast');

% if (options.write)
%     print(theContourFig, fullfile(analysisOutDir,sprintf('CompObsEllipse.tiff')), '-dtiff');
% end

