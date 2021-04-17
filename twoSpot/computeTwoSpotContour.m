function computeTwoSpotContour(options)
% Compute isothreshold contour for mixtures of incremental and decremental spots.
%
% Description:
%    Use ISETBioCSFGenerator to run out an isothreshold contour in the
%    incr-decr contrast plane. This example uses an ideal Poisson 
%    TAFC observer.
%

% History:
%   03/16/21  dhb   Wrote it from color threshold example in ISETBioCSFGenerator.

%% Pick up optional arguments
arguments
    options.defocusDiopters (1,1) double = 0.05;
    options.pupilDiameterMm (1,1) double = 6;
    options.visualizeMosaicResponses (1,1) logical = false;
    options.testing (1,1) logical = false;
    options.write (1,1) logical = true;
    options.verbose (1,1) logical = false;
end

%% Clear and close
close all; ieInit;

%% Directory stuff
baseProject = 'AOCompObserver';
analysisBaseDir = getpref(baseProject,'analysisDir');

%% Testing or running out full computations?
if (options.testing)
    nPixels = 128;
    nTrialsTest = 64;
    angleList = [0 45 135 180 225];
else
    nPixels = 256;
    nTrialsTest = 512;
	angleList = [0 22.5 45 67.5 90 112.5 135 157.5 180 202.5 225 247.5 270 292.5 315 337.5];
end

% Ancillary stimulus parameters
%
% Define basic parameters of the AO stimulus
%
% Measured power was 76 nW for 1.5 by 1 degree field.
% Need to adjust for what it would have been for same
% radiance with a different size field.
wls = 400:10:750;
spotWavelengthNm = 550;
fieldSizeDegs = 0.25;
pupilDiameterMm = options.pupilDiameterMm;
defocusDiopters = options.defocusDiopters;
analysisOutDir = fullfile(analysisBaseDir,sprintf('IncrDecr1_%s_%d',num2str(round(1000*defocusDiopters)),pupilDiameterMm));
if (~exist(analysisOutDir,'dir'))
    mkdir(analysisOutDir);
end

% Set up two spot params
twoSpotParams = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimAngle', 0, ...                             % stimulus angle in incr/decr plane
    'spotWl', spotWavelengthNm, ...                 % spot wavelength: in nm
    'spotFWHM', 20, ...                             % spot full width at half max: in nm
    'spotWidthDegs', 1.3/60, ...                    % spot width: in degrees
    'spotHeightDegs', 1/60, ...                     % spot height: in degrees
    'spotVerticalSepDegs', 1/60, ...                % spot center to center vertical separation: in degrees
    'spotHorizontalSepDegs', 0, ...                 % spot center to center horizontal separation: in degrees
    'spotBgDegs', fieldSizeDegs, ...                % spot background: in degrees
    'spotBgPowerUW', (fieldSizeDegs^2)*76/1000*(2/3), ... % spot background power: in uW
    'imagingWl', 750, ...                           % imaging beam wavelength: in nm
    'imagingFWHM', 20, ...                          % imaging beam full width at half max: in nm
    'imagingBgPowerUW', 0, ...                      % imaging beam power entering eye: in uW
    'fovDegs', fieldSizeDegs, ...                   % spatial: scene field of view, in degrees
    'pixelsNum', nPixels, ...                       % spatial: desired size of stimulus in pixels
    'temporalModulation', 'flashed', ...            % temporal modulation mode: choose between {'flashed'}
    'temporalModulationParams', struct(...          % temporal: modulation params struct
    'stimOnFrameIndices', [1], ...                  %   params relevant to the temporalModulationMode
    'stimDurationFramesNum', 1), ...                %   params relevant to the temporalModulationMode
    'frameDurationSeconds', 3/16, ...               % temporal: frame duration, in seconds
    'pupilDiameterMm', pupilDiameterMm ...          % pupil diameter mm
    );

%% Create neural response engine
%
% This calculations isomerizations in a patch of cone mosaic with Poisson
% noise, and includes optical blur.
%neuralParams = nreAOPhotopigmentExcitationsWithNoEyeMovements;
neuralParams = nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic;

% Set optics params
neuralParams.opticsParams.wls = wls;
neuralParams.opticsParams.pupilDiameterMM = pupilDiameterMm;
neuralParams.opticsParams.defocusAmount = defocusDiopters;
neuralParams.opticsParams.accommodatedWl = spotWavelengthNm;
neuralParams.opticsParams.zCoeffs = zeros(66,1);
neuralParams.opticsParams.defeatLCA = false;
neuralParams.verbose = options.verbose;

% Cone params
neuralParams.coneMosaicParams.fovDegs = fieldSizeDegs;
theNeuralEngine = neuralResponseEngine(@nreAOPhotopigmentExcitationsWithNoEyeMovementsCMosaic, neuralParams);

%% Instantiate the PoissonTAFC responseClassifierEngine
%
% PoissonTAFC makes decision by performing the Poisson likelihood ratio test
% Also set up parameters associated with use of this classifier.
classifierEngine = responseClassifierEngine(@rcePoissonTAFC);
classifierPara = struct('trainFlag', 'none', ...
                        'testFlag', 'random', ...
                        'nTrain', 1, 'nTest', nTrialsTest);

%% Parameters for threshold estimation/quest engine
% The actual threshold varies enough with the different engines that we
% need to adjust the contrast range that Quest+ searches over, as well as
% the range of psychometric function slopes.  Threshold limits are computed
% as 10^-logThreshLimitVal.
thresholdPara = struct('logThreshLimitLow', 3, ...
                       'logThreshLimitHigh', 1, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.5);

% Parameter for running the QUEST+
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 1280, 'maxTrial', 1280, ...
                         'numEstimator', 1, 'stopCriterion', 0.05);
                     
% Create a static two-spot AO scene with a particular incr-decr direction,
% and other relevant parameters
% Compute function handle for two-spot AO stimuli
sceneComputeFunction = @sceTwoSpot;

%% Compute threshold for each spatial direction
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
nDirs = length(angleList);
dataFig = figure();
logThreshold = zeros(1, nDirs);
whichFrameToVisualize = 3;
for ii = 1:nDirs

    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle
    % and the custom params.
    twoSpotParams.stimAngle = angleList(ii);
    twoSpotScene = sceneEngine(sceneComputeFunction, twoSpotParams);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(ii), questObj] = ...
        computeThresholdTAFC(twoSpotScene, theNeuralEngine, classifierEngine, classifierPara, ...
        thresholdPara, questEnginePara, 'visualizeAllComponents', options.visualizeMosaicResponses, ...
        'extraVerbose',options.verbose);
    
    % Plot stimulus
    figure(dataFig);
    subplot(nDirs, 3, ii * 3 - 2);
    visualizationContrast = 1.0;
    [theSceneSequence{ii}] = twoSpotScene.compute(visualizationContrast);
    twoSpotScene.visualizeStaticFrame(theSceneSequence{ii}, ...
        'skipOutOfGamutCheck', true);
    
    % Plot optical image
    subplot(nDirs, 3, ii * 3 - 1);
    visualizationContrast = 1.0;
    [theSceneSequence{ii}] = twoSpotScene.compute(visualizationContrast);
    twoSpotScene.visualizeStaticFrame(theSceneSequence{ii}, ...
        'skipOutOfGamutCheck', true, ...
        'opticalImageInsteadOfScene', theNeuralEngine.neuralPipeline.optics);
    
    % Plot data and psychometric curve 
    % with a marker size of 2.5
    subplot(nDirs, 3, ii * 3 );
    questObj.plotMLE(2.5);
    drawnow;
end
set(dataFig, 'Position',  [0, 0, 800, 800]);
if (options.write)
    print(dataFig,fullfile(analysisOutDir,sprintf('CompObsPsychoFig.tiff')), '-dtiff');
end

% Convert returned log threshold to linear threshold
threshold = 10 .^ logThreshold;

% Threshold cone contrasts
thresholdContrasts = [threshold.*cosd(angleList) ; threshold.*sind(angleList)];

%% Fit an ellipse to the data.  See EllipseTest and FitEllipseQ.
%
% The use of scaleFactor to scale up the data and scale down the fit by the
% same amount is fmincon black magic.  Doing this puts the objective
% function into a better range for the default size of search steps.
%
% We constrain the ellipse to line up with the x and y axes.  Change flag
% below to relax this.  Doesn't make very much difference inthis case.
scaleFactor = 1/(max(abs(thresholdContrasts(:))));
[compEllParams,compFitA,compFitAinv,compFitQ] = FitEllipseQ(scaleFactor*thresholdContrasts(1:2,:),'lockAngleAt0',false, ...
    'initialParams',[0.9 1.1 45]);
nThetaEllipse = 200;
circleIn2D = UnitCircleGenerate(nThetaEllipse);
fitEllipse = PointsOnEllipseQ(compFitQ,circleIn2D)/scaleFactor;

%% Plot
contrastLim = 0.05;
theContourFig = figure; clf; hold on
plot(fitEllipse(1,:),fitEllipse(2,:),'r','LineWidth',3);
plot([-contrastLim contrastLim],[0 0],'k:','LineWidth',1);
plot([0 0],[-contrastLim contrastLim],'k:','LineWidth',1);
xlabel('Contrast 1');
ylabel('Contrsast 2');
set(theContourFig, 'Position',  [800, 0, 600, 800]);
xlim([-contrastLim contrastLim]); ylim([-contrastLim contrastLim]);
axis('square');
if (options.write)
    print(theContourFig, fullfile(analysisOutDir,sprintf('CompObsEllipse.tiff')), '-dtiff');
end

% Add comp data to plot
plot(thresholdContrasts(1,:), thresholdContrasts(2,:), 'ok', 'MarkerFaceColor','k', 'MarkerSize',12);
if (options.write)
    print(theContourFig, fullfile(analysisOutDir,sprintf('CompObsEllipseWithData.tiff')), '-dtiff');
end

%% Save
close(dataFig);
close(theContourFig);
if (options.write)
    save(fullfile(analysisOutDir,'CompObserver'));
end
