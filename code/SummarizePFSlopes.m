
%% Clear
clear; close all;

%% Output variant
outputVariant = 'SlopeFree1';

%% Parameters
axisFontSize = 14;
labelFontSize = 18;
titleFontSize = 20;
legendFontSize = 14;
markerSize = 12;

%% Add parameters so we can check limits.
%
% Copied over from FitTrials.
convertToDb = true;
switch (outputVariant)
    case 'SlopeFree'
        if (convertToDb)
            slopeLower8 = 0.5;
            slopeUpper8 = 7;
            slopeLower43 = 0.5;
            slopeUpper43 = 7;
        else
            slopeLower8 = 5;
            slopeUpper8 = 70;
            slopeLower43 = 5;
            slopeUpper43 = 70;
        end
        lapseUpper = 0.05;
    case 'SlopeFree1'
        if (convertToDb)
            slopeLower8 = 0.5;
            slopeUpper8 = 7;
            slopeLower43 = 0.5;
            slopeUpper43 = 7;
        else
            slopeLower8 = 5;
            slopeUpper8 = 70;
            slopeLower43 = 5;
            slopeUpper43 = 70;
        end
        lapseUpper = 0.01;
    case 'SlopeFixed'
        if (convertToDb)
            slopeLower8 = 1.08;
            slopeUpper8 = 1.08;
            slopeLower43 = 2.69;
            slopeUpper43 = 2.69;
        else
            slopeLower8 = 10.8;
            slopeUpper8 = 10.8;
            slopeLower43 = 26.9;
            slopeUpper43 = 26.9;
        end
        lapseUpper = 0.01;
end
guessUpper = 0.05;

%% Path to analysis output dir
analysisDir = getpref('AOMicroRepeat','analysisDir');

%% Read table of fit information
wState = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet'); %,'VariableNamingRule','preserve');
warning(wState);

%% Get all MOCS beta for full sessions at 8 arcmin
theDiameter = 8;
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = find(theDiameters == theDiameter)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'MOCS') & strcmp(dataTable.Split,'All'));
            betaMOCS(pp,dd,ss) = dataTable.Beta(index);
            lapseFromFitMOCS(pp,dd,ss) = dataTable.Lapse(index);
            lapseFromDataMOCS(pp,dd,ss) = dataTable.lapseEstFromData(index);
            lapseNTrialsFromDataMOCS(pp,dd,ss) = dataTable.nTrialsLapseEstFromData(index);
            guessFromFitMOCS(pp,dd,ss) = dataTable.Guess(index);
            guessFromDataMOCS(pp,dd,ss) = dataTable.guessEstFromData(index);
            guessNTrialsFromDataMOCS(pp,dd,ss) = dataTable.nTrialsGuessEstFromData(index);

            if (lapseFromFitMOCS(pp,dd,ss) >= lapseUpper)
                fprintf('*** MOCS lapse rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitMOCS(pp,dd,ss),lapseUpper,pp,dd,ss);
            end
            if (guessFromFitMOCS(pp,dd,ss) >= guessUpper)
                fprintf('*** MOCS guess rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitMOCS(pp,dd,ss),guessUpper,pp,dd,ss);
            end
            if (betaMOCS(pp,dd,ss) <= slopeLower8)
                fprintf('*** MOCS slope size 8 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeLower8,pp,dd,ss);
            end
            if (betaMOCS(pp,dd,ss) >= slopeUpper8)
                fprintf('*** MOCS slope size 8 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeUpper8,pp,dd,ss);
            end
        end
    end
end
fprintf('MOCS, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n\n', ...
    theDiameter,mean(betaMOCS(:)),std(betaMOCS(:)), ...
    mean(lapseFromFitMOCS(:)),mean(lapseFromDataMOCS(:)),mean(lapseNTrialsFromDataMOCS(:)), ...
    mean(guessFromFitMOCS(:)),mean(guessFromDataMOCS(:)),mean(guessNTrialsFromDataMOCS(:)));

%% Get all QUEST beta for full sessions at 8 arcmin
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = find(theDiameters == theDiameter)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'QUEST') & strcmp(dataTable.Split,'All'));
            betaQUEST(pp,dd,ss) = dataTable.Beta(index);
            lapseFromFitQUEST(pp,dd,ss) = dataTable.Lapse(index);
            lapseFromDataQUEST(pp,dd,ss) = dataTable.lapseEstFromData(index);
            lapseNTrialsFromDataQUEST(pp,dd,ss) = dataTable.nTrialsLapseEstFromData(index);
            guessFromFitQUEST(pp,dd,ss) = dataTable.Guess(index);
            guessFromDataQUEST(pp,dd,ss) = dataTable.guessEstFromData(index);
            guessNTrialsFromDataQUEST(pp,dd,ss) = dataTable.nTrialsGuessEstFromData(index);

            if (lapseFromFitQUEST(pp,dd,ss) >= lapseUpper)
                fprintf('*** QUEST lapse rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitQUEST(pp,dd,ss),lapseUpper,pp,dd,ss);
            end
            if (guessFromFitQUEST(pp,dd,ss) >= guessUpper)
                fprintf('*** QUEST guess rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitQUEST(pp,dd,ss),guessUpper,pp,dd,ss);
            end
            if (betaQUEST(pp,dd,ss) <= slopeLower8)
                fprintf('*** QUEST slope size 8 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeLower8,pp,dd,ss);
            end
            if (betaQUEST(pp,dd,ss) >= slopeUpper8)
                fprintf('*** QUEST slope size 8 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeUpper8,pp,dd,ss);
            end
        end
    end
end
index = find(lapseNTrialsFromDataQUEST(:) >= 10);
fprintf('QUEST, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n\n', ...
    theDiameter,mean(betaQUEST(:)),std(betaQUEST(:)), ...
    mean(lapseFromFitQUEST(:)),mean(lapseFromDataQUEST(index)),mean(lapseNTrialsFromDataQUEST(index)), ...
    mean(guessFromFitQUEST(:)),mean(guessFromDataQUEST(:)),mean(guessNTrialsFromDataQUEST(:)));

%% Get all MOCS beta for full sessions at 43 arcmin
theDiameter = 43;
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = find(theDiameters == theDiameter)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'MOCS') & strcmp(dataTable.Split,'All'));
            betaMOCS(pp,dd,ss) = dataTable.Beta(index);
            lapseFromFitMOCS(pp,dd,ss) = dataTable.Lapse(index);
            lapseFromDataMOCS(pp,dd,ss) = dataTable.lapseEstFromData(index);
            lapseNTrialsFromDataMOCS(pp,dd,ss) = dataTable.nTrialsLapseEstFromData(index);
            guessFromFitMOCS(pp,dd,ss) = dataTable.Guess(index);
            guessFromDataMOCS(pp,dd,ss) = dataTable.guessEstFromData(index);
            guessNTrialsFromDataMOCS(pp,dd,ss) = dataTable.nTrialsGuessEstFromData(index);

            if (lapseFromFitMOCS(pp,dd,ss) >= lapseUpper)
                fprintf('*** MOCS lapse rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitMOCS(pp,dd,ss),lapseUpper,pp,dd,ss);
            end
            if (guessFromFitMOCS(pp,dd,ss) >= guessUpper)
                fprintf('*** MOCS guess rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitMOCS(pp,dd,ss),guessUpper,pp,dd,ss);
            end
            if (betaMOCS(pp,dd,ss) <= slopeLower43)
                fprintf('*** MOCS slope size 43 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeLower43,pp,dd,ss);
            end
            if (betaMOCS(pp,dd,ss) >= slopeUpper43)
                fprintf('*** MOCS slope size 43 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeUpper43,pp,dd,ss);
            end
        end
    end
end
fprintf('MOCS, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n\n', ...
    theDiameter,mean(betaMOCS(:)),std(betaMOCS(:)), ...
    mean(lapseFromFitMOCS(:)),mean(lapseFromDataMOCS(:)),mean(lapseNTrialsFromDataMOCS(:)), ...
    mean(guessFromFitMOCS(:)),mean(guessFromDataMOCS(:)),mean(guessNTrialsFromDataMOCS(:)));

%% Get all QUEST beta for full sessions at 43 arcmin
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = find(theDiameters == theDiameter)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'QUEST') & strcmp(dataTable.Split,'All'));
            betaQUEST(pp,dd,ss) = dataTable.Beta(index);
            lapseFromFitQUEST(pp,dd,ss) = dataTable.Lapse(index);
            lapseFromDataQUEST(pp,dd,ss) = dataTable.lapseEstFromData(index);
            lapseNTrialsFromDataQUEST(pp,dd,ss) = dataTable.nTrialsLapseEstFromData(index);
            guessFromFitQUEST(pp,dd,ss) = dataTable.Guess(index);
            guessFromDataQUEST(pp,dd,ss) = dataTable.guessEstFromData(index);
            guessNTrialsFromDataQUEST(pp,dd,ss) = dataTable.nTrialsGuessEstFromData(index);

            if (lapseFromFitQUEST(pp,dd,ss) >= lapseUpper)
                fprintf('*** QUEST lapse rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitQUEST(pp,dd,ss),lapseUpper,pp,dd,ss);
            end
            if (guessFromFitQUEST(pp,dd,ss) >= guessUpper)
                fprintf('*** QUEST guess rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitQUEST(pp,dd,ss),guessUpper,pp,dd,ss);
            end
            if (betaQUEST(pp,dd,ss) <= slopeLower43)
                fprintf('*** QUEST slope size 43 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeLower43,pp,dd,ss);
            end
            if (betaQUEST(pp,dd,ss) >= slopeUpper43)
                fprintf('*** QUEST slope size 43 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeUpper43,pp,dd,ss);
            end
        end
    end
end
index = find(lapseNTrialsFromDataQUEST(:) >= 10);
fprintf('QUEST, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n\n', ...
    theDiameter,mean(betaQUEST(:)),std(betaQUEST(:)), ...
    mean(lapseFromFitQUEST(:)),mean(lapseFromDataQUEST(index)),mean(lapseNTrialsFromDataQUEST(index)), ...
    mean(guessFromFitQUEST(:)),mean(guessFromDataQUEST(:)),mean(guessNTrialsFromDataQUEST(:)));
