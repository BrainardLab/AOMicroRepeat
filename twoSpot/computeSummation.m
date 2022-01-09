function computeSummation(options)
% Compute summation curve for circular spots.
%
% Description:
%    Use ISETBioCSFGenerator to run out a summation curve for circular
%    spots.
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
    options.smoothingParam (1,1) double = 0.995;
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
    minTrials = 1024;
    diameterList = [0.25 0.5 1 2]/60;
else
    nPixels = 256;
    nTrialsTest = 32;
    minTrials = 5000;
	diameterList = logspace(log10(0.15),log10(2),20)/60;
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
analysisOutDir = fullfile(analysisBaseDir,sprintf('Summation_%s_%d',num2str(round(1000*defocusDiopters)),pupilDiameterMm));
if (~exist(analysisOutDir,'dir'))
    mkdir(analysisOutDir);
end

% Set up two spot params
circularSpotParams = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimDiameter', 0, ...                          % stimulus spot diameter in degrees.
    'spotWl', spotWavelengthNm, ...                 % spot wavelength: in nm
    'spotFWHM', 20, ...                             % spot full width at half max: in nm
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
                       'logThreshLimitHigh', 0, ...
                       'logThreshLimitDelta', 0.02, ...
                       'slopeRangeLow', 1, ...
                       'slopeRangeHigh', 50, ...
                       'slopeDelta', 2.5);

% Parameter for running the QUEST+
%
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive). See t_thresholdEngine and
% questThresholdEngine for more info.
%
% Here we set maxTrial = minTrial, and don't use an additinal stop
% criterion, so fundamentally we just run minTrials trials. This is total
% trials evaluated, not number of separate constrasts.  So if you increase 
% nTrialsTest, you'll get fewer contrasts tested unless you also increase
% minTrials in proportion.
questEnginePara = struct('minTrial', minTrials, 'maxTrial', minTrials, ...
                         'numEstimator', 1, 'stopCriterion', []);
                     
% Create a static two-spot AO scene with a particular incr-decr direction,
% and other relevant parameters
% Compute function handle for two-spot AO stimuli
sceneComputeFunction = @sceAOCircularSpot;

%% Compute threshold for each spatial direction
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
nDirs = length(diameterList);
dataFig = figure();
logThreshold = zeros(1, nDirs);

for ii = 1:nDirs

    % Instantiate a sceneEngine with the above sceneComputeFunctionHandle
    % and the custom params.
    circularSpotParams.stimDiameter = diameterList(ii);
    circularSpotScene = sceneEngine(sceneComputeFunction, circularSpotParams);
    
    % Compute the threshold for our grating scene with the previously
    % defined neural and classifier engine.  This function does a lot of
    % work, see t_tresholdEngine and the function itself, as well as
    % function computePerformanceTAFC.
    [logThreshold(ii), questObj] = ...
        computeThresholdTAFC(circularSpotScene, theNeuralEngine, classifierEngine, classifierPara, ...
        thresholdPara, questEnginePara, 'visualizeAllComponents', options.visualizeMosaicResponses, ...
        'extraVerbose',options.verbose);
    
    % Plot stimulus
    figure(dataFig);
    subplot(nDirs, 3, ii * 3 - 2);
    visualizationContrast = 1.0;
    [theSceneSequence{ii}] = circularSpotScene.compute(visualizationContrast);
    circularSpotScene.visualizeStaticFrame(theSceneSequence{ii}, ...
        'skipOutOfGamutCheck', true);
    
    % Plot optical image
    subplot(nDirs, 3, ii * 3 - 1);
    visualizationContrast = 1.0;
    [theSceneSequence{ii}] = circularSpotScene.compute(visualizationContrast);
    circularSpotScene.visualizeStaticFrame(theSceneSequence{ii}, ...
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

% Threshold energy
thresholdEnergy = pi*((diameterList/2).^2).*threshold;

%% Smooth curve through the data.  MATLAB's smoothing spline.
% Vary smoothness parameter to control.  Bigger is less smooth.
% Want this pretty big.
diameterPlot = linspace(diameterList(1),diameterList(end),100);
fObj = fit(log10(60*diameterList)',log10(thresholdEnergy)','smoothingspline','smoothingparam',options.smoothingParam);
log10ThresholdEnergySmooth = feval(fObj,log10(60*diameterPlot));

%% Plot
theSummationFig = figure; clf; hold on
plot(log10(60*diameterList),log10(thresholdEnergy),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(log10(60*diameterPlot),log10ThresholdEnergySmooth,'b','LineWidth',3);
xlabel('Log Spot Diameter (minutes)');
ylabel('Log Threshold Energy');
set(theSummationFig, 'Position',  [800, 0, 600, 800]);
if (options.write)
    print(theSummationFig, fullfile(analysisOutDir,sprintf('CompSummation.tiff')), '-dtiff');
end

%% Save
close(dataFig);
close(theSummationFig);
if (options.write)
    save(fullfile(analysisOutDir,'CompObserver'));
end

end

