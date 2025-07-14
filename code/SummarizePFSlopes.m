
%% Clear
clear; close all;

%% Output variant
outputVariant = 'SlopeFree';

%% Parameters
axisFontSize = 14;
labelFontSize = 18;
titleFontSize = 20;
legendFontSize = 14;
markerSize = 12;

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
        end
    end
end
fprintf('MOCS, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n', ...
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
        end
    end
end
index = find(lapseNTrialsFromDataQUEST(:) >= 10);
fprintf('QUEST, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n', ...
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
        end
    end
end
fprintf('MOCS, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n', ...
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
        end
    end
end
index = find(lapseNTrialsFromDataQUEST(:) >= 10);
fprintf('QUEST, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n', ...
    theDiameter,mean(betaQUEST(:)),std(betaQUEST(:)), ...
    mean(lapseFromFitQUEST(:)),mean(lapseFromDataQUEST(index)),mean(lapseNTrialsFromDataQUEST(index)), ...
    mean(guessFromFitQUEST(:)),mean(guessFromDataQUEST(:)),mean(guessNTrialsFromDataQUEST(:)));
