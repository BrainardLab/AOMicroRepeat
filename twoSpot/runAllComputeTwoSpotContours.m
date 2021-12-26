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
    defocusDioptersList = [0 0.05 0.10 0.15];
    pupilDiameterMmList = [7];
    testing = false;
    write = true;
    verbose = false;
end

%% Degrees per pixel
degsPerPixel = 1/415;

%% 7x9, various separations, optics
spotHeightDegs = 7*degsPerPixel;
spotWidthDegs = 9*degsPerPixel;
spotVertSepDegs = spotHeightDegs + 0*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_0','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 1*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_1','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 2*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_2','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 3*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_3','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 4*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_4','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 5*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_5','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 6*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_6','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 8*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_8','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 12*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_12','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 16*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_16','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 20*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_20','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 24*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','7_9_24','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

%% 6x8
spotHeightDegs = 6*degsPerPixel;
spotWidthDegs = 8*degsPerPixel;
spotVertSepDegs = spotHeightDegs + 0*degsPerPixel;
% for dd = 1:length(defocusDioptersList)
%     for pp = 1:length(pupilDiameterMmList)
%         computeTwoSpotContour('conditionName','6_8_0','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
%             'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
%             'testing',testing,'write', write,'verbose',verbose);
%     end
% end

spotVertSepDegs = spotHeightDegs + 1*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_0','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 2*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_2','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 3*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_3','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

% spotVertSepDegs = spotHeightDegs + 4*degsPerPixel;
% for dd = 1:length(defocusDioptersList)
%     for pp = 1:length(pupilDiameterMmList)
%         computeTwoSpotContour('conditionName','6_8_4','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
%             'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
%             'testing',testing,'write', write,'verbose',verbose);
%     end
% end

spotVertSepDegs = spotHeightDegs + 6*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_6','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 8*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_8','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 12*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_12','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 16*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_16','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 20*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_20','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

spotVertSepDegs = spotHeightDegs + 24*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','6_8_24','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end

%% 5x7
spotHeightDegs = 5*degsPerPixel;
spotWidthDegs = 7*degsPerPixel;
spotVertSepDegs = spotHeightDegs + 0*degsPerPixel;
for dd = 1:length(defocusDioptersList)
    for pp = 1:length(pupilDiameterMmList)
        computeTwoSpotContour('conditionName','5_7_0','defocusDiopters',defocusDioptersList(dd),'pupilDiameterMm',pupilDiameterMmList(pp), ...
            'spotWidthDegs',spotWidthDegs ,'spotHeightDegs',spotHeightDegs,'spotVerticalSepDegs',spotVertSepDegs, ...
            'testing',testing,'write', write,'verbose',verbose);
    end
end