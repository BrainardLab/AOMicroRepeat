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
    options.fFactor (1,1) double = 1e5;
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
    diameterList = [0.25 0.5 1 2]/60;
else
    nPixels = 256;
    nTrialsTest = 512;
	diameterList = linspace(0.15,2,20)/60;
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
% See t_thresholdEngine.m for more on options of the two different mode of
% operation (fixed numer of trials vs. adaptive)
questEnginePara = struct('minTrial', 1280, 'maxTrial', 1280, ...
                         'numEstimator', 1, 'stopCriterion', 0.05);
                     
% Create a static two-spot AO scene with a particular incr-decr direction,
% and other relevant parameters
% Compute function handle for two-spot AO stimuli
sceneComputeFunction = @sceCircularSpot;

%% Compute threshold for each spatial direction
% 
% See toolbox/helpers for functions createGratingScene computeThresholdTAFC
nDirs = length(diameterList);
dataFig = figure();
logThreshold = zeros(1, nDirs);
whichFrameToVisualize = 3;
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

%% Smooth curve through the data
criterionDPrime = 1;
sigmaExternal = diameterList(end);

sigmaInternal0 = thresholdEnergy(1);
externalIntrusionFactor0 = sigmaInternal0/diameterList(2);
signalExponent0 = 1;
x0Out = [sigmaInternal0 externalIntrusionFactor0 signalExponent0];
xFactor0 = 1 ./ x0Out;
x0 = x0Out .* xFactor0;

% Set reasonable bounds on parameters
vlb0 = [1e-8 1e-8 1];
vub0 = [1e6 1e6 1];

vlb = vlb0 .* xFactor0;
vub = vub0 .* xFactor0;
minoptions = optimset('fmincon');
minoptions = optimset(options,'Diagnostics','off','Display','iter','LargeScale','off','Algorithm','active-set');

% Fit with exponent locked
x1Temp = fmincon(@(x)FitTvNFunction(x,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor0,options.fFactor),x0,[],[],[],[],vlb,vub,[],minoptions);
f1 = FitTvNFunction(x1Temp,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor0,options.fFactor);
x1Out = x1Temp ./ xFactor0;
xFactor1 = 1 ./ x1Out;
x1 = x1Out .* xFactor1;
f1Out = FitTvNFunction(x1,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor1,options.fFactor);

% Exponent free
vlb0(3) = 0.33;
vub0(3) = 3;
vlb = vlb0 .* xFactor1;
vub = vub0 .* xFactor1;
x2Temp = fmincon(@(x)FitTvNFunction(x,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor1,options.fFactor),x1,[],[],[],[],vlb,vub,[],minoptions);
x2Out = x2Temp ./ xFactor1;
x2Out(1) = thresholdEnergy(1);
xFactor2 = 1 ./ x2Out;
x2 = x2Out .* xFactor2;
f2Out = FitTvNFunction(x2,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor2,options.fFactor);
xOut = x2 ./ xFactor2;

% Once more with exponent locked
%vlb0(3) = x2(3) ./ xFactor2(3);
%vub0(3) = x2(3) ./ xFactor2(3);
vlb = vlb0 .* xFactor2;
vub = vub0 .* xFactor2;
x3Temp = fmincon(@(x)FitTvNFunction(x2,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor2,options.fFactor),x2,[],[],[],[],vlb,vub,[],minoptions);
x3Out = x3Temp ./ xFactor2;
xFactor3 = 1 ./ x3Out;
x3 = x3Out .* xFactor3;
f3Out = FitTvNFunction(x3,thresholdEnergy,diameterList,sigmaExternal,criterionDPrime,xFactor3,options.fFactor);
xOut = x3 ./ xFactor3;

% Get parameters back out again
sigmaInternal = xOut(1);
externalIntrusionFactor = xOut(2);
signalExponent = xOut(3);

% Predict
diameterPlot = linspace(diameterList(1),diameterList(end),100);
smoothThresholdEnergy = ComputeTvNThreshold(diameterPlot,criterionDPrime,sigmaInternal,externalIntrusionFactor,signalExponent);

%% Plot
theSummationFig = figure; clf; hold on
plot(log10(60*diameterList),log10(thresholdEnergy),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(log10(60*diameterPlot),log10(smoothThresholdEnergy),'b','LineWidth',3);
xlabel('Log Spot Diameter (minutes)');
ylabel('Log Threshold Energy');
set(theSummationFig, 'Position',  [800, 0, 600, 800]);
if (options.write)
    print(theSummationFig, fullfile(analysisOutDir,sprintf('CompObsEllipse.tiff')), '-dtiff');
end

%% Save
close(dataFig);
close(theSummationFig);
if (options.write)
    save(fullfile(analysisOutDir,'CompObserver'));
end

end

%% Error function for fit
function f = FitTvNFunction(x,thresholds,sigmaExternals,sigmaExternal,criterionDPrime,xFactor,fFactor)

    x = x ./ xFactor;
    sigmaInternal = x(1);
    externalIntrusionFactor = x(2);
    signalExponent = x(3);
    
    predictedThresholds = ComputeTvNThreshold(sigmaExternals,criterionDPrime,sigmaInternal,externalIntrusionFactor,signalExponent);
    f = sqrt(mean(log10(thresholds)-log10(predictedThresholds)).^2);
    f = f*fFactor;

end

% Simple parametric form for TvN
function thresholds = ComputeTvNThreshold(sigmaExternals,criterionDPrime,sigmaInternal,externalIntrusionFactor,signalExponent)
    sigma = sqrt(sigmaInternal^2 + sigmaExternals.^2*externalIntrusionFactor);
    exponentiatedThreshold = criterionDPrime*sigma;
    thresholds = exponentiatedThreshold.^(1/signalExponent);
end
