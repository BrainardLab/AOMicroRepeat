% Run out a bunch of computational observers
%
% Description:
%   Run out the computational observer for the circular spot summation experiment with a
%   variety of defocus and pupil sizes.

% Set parameters depending on mode
testMode = false;
if (testMode)
    defocusDioptersList = [0];
    pupilDiameterMmList = [7];
    testing = true;
    write = false;
    verbose = true;
else
    defocusDioptersList = [0 0.05 0.10 0.15];
    pupilDiameterMmList = [2 3 4 5 6 7];
    defocusDioptersList = [0];
    pupilDiameterMmList = [7];
    testing = false;
    write = true;
    verbose = false;
end

% Do them
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeSummation('defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end