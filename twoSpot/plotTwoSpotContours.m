% Make some plots of the computational observer results
%
% Description:
%   Run out the computational observer for the two spot experiment with a
%   variety of defocus and pupil sizes.

% Paramters
defocusDioptersList = [0 0.05 0.10 0.15];
pupilDiameterMmList = [2 3 4 5 6 7];
testing = false;

% Do them
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        theData = computeTwoSpotContour('defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'testing',testing);
    end
end