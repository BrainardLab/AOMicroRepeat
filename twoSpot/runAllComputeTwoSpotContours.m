% Run out a bunch of computational observers
%
% Description:
%   Run out the computational observer for the two spot experiment with a
%   variety of defocus and pupil sizes.

%% Set parameters depending on mode
testMode = false;
if (testMode)
    defocusDioptersList = [0];
    pupilDiameterMmList = [7];
    testing = true;
    write = false;
    verbose = true;
else
    defocusDioptersList = [0 0.05 0.10];
    pupilDiameterMmList = [7];
    testing = false;
    write = true;
    verbose = false;
end

%% 2020-01-31
%
% 7x9
degsPerPixel = (1/7)*(1/60);
spotWidthDegs = 9*degsPerPixel;
spotHeightDegs = 7*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','IncrDecr1','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotHeightDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

%% 2021-09-14
%
% 7x9
degsPerPixel = (1/7)*(1/60);
spotWidthDegs = 9*degsPerPixel;
spotHeightDegs = 7*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','IncrDecr2','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotHeightDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

%% 2021-10-18
%
% 5x7
degsPerPixel = (1/7)*(1/60);
spotWidthDegs = 7*degsPerPixel;
spotHeightDegs = 5*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','IncrDecr3','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotHeightDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end