% Run out a bunch of computational observers
%
% Description:
%   Run out the computational observer for the two spot experiment with a
%   variety of defocus and pupil sizes.
%
% NOTE: When updating, need to change directory strings as
%       well as numbers used in the computations.  It would
%       be better to generate the strings from the numbers.

%% Set parameters depending on mode
testMode = true;
if (testMode)
    defocusDioptersList = [0.05];
    pupilDiameterMmList = [7];
    testing = true;
    write = false;
    verbose = true;
else
    defocusDioptersList = [0 0.05 0.10 0.15];
    pupilDiameterMmList = [7];
    testing = false;
    write = true;
    verbose = false;
end

%% Degrees per pixel
degsPerPixel = 1/415;

% 7x9, various separations, optics
spotHeightDegs = 7*degsPerPixel;
spotWidthDegs = 9*degsPerPixel;
spotVertSepDegs = spotHeightDegs + 0*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_0','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose,'visualizeMosaicResponses',true,'visualizeStimulus',false,'angleList',330);
    end
end

