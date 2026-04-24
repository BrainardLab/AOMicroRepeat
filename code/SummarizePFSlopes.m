
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
    nPFs(pp) = 0;
    nGuessesAtLimit(pp) = 0;
    nLapsesAtLimit(pp) = 0;
    nBetasAtLowerLimit(pp) = 0;
    nBetasAtUpperLimit(pp) = 0;
end

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

            nPFs(pp) = nPFs(pp) + 1;
            if (lapseFromFitMOCS(pp,dd,ss) >= lapseUpper)
                fprintf('*** MOCS lapse rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitMOCS(pp,dd,ss),lapseUpper,pp,dd,ss);
                nLapsesAtLimit(pp) = nLapsesAtLimit(pp) + 1;
            end
            if (guessFromFitMOCS(pp,dd,ss) >= guessUpper)
                fprintf('*** MOCS guess rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitMOCS(pp,dd,ss),guessUpper,pp,dd,ss);
                nGuessesAtLimit(pp) = nGuessesAtLimit(pp) + 1;
            end
            if (betaMOCS(pp,dd,ss) <= slopeLower8)
                fprintf('*** MOCS slope size 8 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeLower8,pp,dd,ss);
                nBetasAtLowerLimit(pp) = nBetasAtLowerLimit(pp) + 1;
            end
            if (betaMOCS(pp,dd,ss) >= slopeUpper8)
                fprintf('*** MOCS slope size 8 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeUpper8,pp,dd,ss);
                nBetasAtUpperLimit(pp) = nBetasAtUpperLimit(pp) + 1;
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

            nPFs(pp) = nPFs(pp) + 1;
            if (lapseFromFitQUEST(pp,dd,ss) >= lapseUpper)
                fprintf('*** QUEST lapse rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitQUEST(pp,dd,ss),lapseUpper,pp,dd,ss);
                nLapsesAtLimit(pp) = nLapsesAtLimit(pp) + 1;
            end
            if (guessFromFitQUEST(pp,dd,ss) >= guessUpper)
                fprintf('*** QUEST guess rate size 8 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitQUEST(pp,dd,ss),guessUpper,pp,dd,ss);
                nGuessesAtLimit(pp) = nGuessesAtLimit(pp) + 1;
            end
            if (betaQUEST(pp,dd,ss) <= slopeLower8)
                fprintf('*** QUEST slope size 8 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeLower8,pp,dd,ss);
                 nBetasAtLowerLimit(pp) = nBetasAtLowerLimit(pp) + 1;
            end
            if (betaQUEST(pp,dd,ss) >= slopeUpper8)
                fprintf('*** QUEST slope size 8 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeUpper8,pp,dd,ss);
                nBetasAtUpperLimit(pp) = nBetasAtUpperLimit(pp) + 1;
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

            nPFs(pp) = nPFs(pp) + 1;
            if (lapseFromFitMOCS(pp,dd,ss) >= lapseUpper)
                fprintf('*** MOCS lapse rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitMOCS(pp,dd,ss),lapseUpper,pp,dd,ss);
                nLapsesAtLimit(pp) = nLapsesAtLimit(pp) + 1;
            end
            if (guessFromFitMOCS(pp,dd,ss) >= guessUpper)
                fprintf('*** MOCS guess rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitMOCS(pp,dd,ss),guessUpper,pp,dd,ss);
                nGuessesAtLimit(pp) = nGuessesAtLimit(pp) + 1;
            end
            if (betaMOCS(pp,dd,ss) <= slopeLower43)
                fprintf('*** MOCS slope size 43 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeLower43,pp,dd,ss);
                 nBetasAtLowerLimit(pp) = nBetasAtLowerLimit(pp) + 1;
            end
            if (betaMOCS(pp,dd,ss) >= slopeUpper43)
                fprintf('*** MOCS slope size 43 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaMOCS(pp,dd,ss),slopeUpper43,pp,dd,ss);
                nBetasAtUpperLimit(pp) = nBetasAtUpperLimit(pp) + 1;
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

            nPFs(pp) = nPFs(pp) + 1;
            if (lapseFromFitQUEST(pp,dd,ss) >= lapseUpper)
                fprintf('*** QUEST lapse rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    lapseFromFitQUEST(pp,dd,ss),lapseUpper,pp,dd,ss);
                nLapsesAtLimit(pp) = nLapsesAtLimit(pp) + 1;
            end
            if (guessFromFitQUEST(pp,dd,ss) >= guessUpper)
                fprintf('*** QUEST guess rate size 43 %0.2f at or greater than bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    guessFromFitQUEST(pp,dd,ss),guessUpper,pp,dd,ss);
                nGuessesAtLimit(pp) = nGuessesAtLimit(pp) + 1;
            end
            if (betaQUEST(pp,dd,ss) <= slopeLower43)
                fprintf('*** QUEST slope size 43 %0.2f at or less than lower bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeLower43,pp,dd,ss);
                 nBetasAtLowerLimit(pp) = nBetasAtLowerLimit(pp) + 1;
            end
            if (betaQUEST(pp,dd,ss) >= slopeUpper43)
                fprintf('*** QUEST slope size 43 %0.2f at or greater than upper bound of %0.2f, pp = %d, dd = %d, ss = %d\n', ...
                    betaQUEST(pp,dd,ss),slopeUpper43,pp,dd,ss);
                nBetasAtUpperLimit(pp) = nBetasAtUpperLimit(pp) + 1;
            end
        end
    end
end
index = find(lapseNTrialsFromDataQUEST(:) >= 10);
fprintf('QUEST, size %d, mean beta %0.2f, stdev %0.3f, lapse from fit %0.3f, lapse from data %0.3f (%0.1f), guess from fit %0.3f, guess from data %0.3f (%0.1f)\n\n', ...
    theDiameter,mean(betaQUEST(:)),std(betaQUEST(:)), ...
    mean(lapseFromFitQUEST(:)),mean(lapseFromDataQUEST(index)),mean(lapseNTrialsFromDataQUEST(index)), ...
    mean(guessFromFitQUEST(:)),mean(guessFromDataQUEST(:)),mean(guessNTrialsFromDataQUEST(:)));

for pp = 1:length(theSubjects)
    fprintf('Subject %d\n',pp)
    fprintf('\tLapses at limit %0.2f%%\n',round(100*nLapsesAtLimit(pp)/nPFs(pp),2));
    fprintf('\tGuesses at limit %0.2f%%\n',round(100*nGuessesAtLimit(pp)/nPFs(pp),2));
    fprintf('\tBetas at lower %0.2f%%\n',round(100*nBetasAtLowerLimit(pp)/nPFs(pp),2));
    fprintf('\tBetas at upper %0.2f%%\n',round(100*nBetasAtUpperLimit(pp)/nPFs(pp),2));
end

% Write out guess table values in form they can go into Excel table
rateTablefile = fopen(fullfile(analysisDir,outputVariant,'tableS2.txt'),"w");
tempMOCS = round(guessFromDataMOCS.*guessNTrialsFromDataMOCS);
tempQUEST = round(guessFromDataQUEST.*guessNTrialsFromDataQUEST);
for ss = 1:size(tempMOCS,2)
    fprintf(rateTablefile,'\nFalse Positive Rates, Session %d\n',ss);
    for pp = 1:length(theSubjects)
        fprintf(rateTablefile,'%s\t%d/%d\t%d/%d\t%d/%d\t%d/%d\t%d/%d (%0.1f%%)\n', ...
            theSubjects{pp}, ...
            tempMOCS(pp,1,ss), guessNTrialsFromDataMOCS(pp,1,ss), ...
            tempQUEST(pp,1,ss), guessNTrialsFromDataQUEST(pp,1,ss), ...
            tempMOCS(pp,2,ss), guessNTrialsFromDataMOCS(pp,2,ss), ...
            tempQUEST(pp,2,ss), guessNTrialsFromDataQUEST(pp,2,ss), ...
            sum(tempMOCS(pp,:,ss)+tempQUEST(pp,:,ss)), ...
            sum(guessNTrialsFromDataMOCS(pp,:,ss)+guessNTrialsFromDataQUEST(pp,:,ss)), ...
            100*sum(tempMOCS(pp,:,ss)+tempQUEST(pp,:,ss)) / ...
            sum(guessNTrialsFromDataMOCS(pp,:,ss)+guessNTrialsFromDataQUEST(pp,:,ss)) ...
            );
    end
    tempTempMOCS = squeeze(sum(tempMOCS(:,:,:)));
    tempTempQUEST = squeeze(sum(tempQUEST(:,:,:)));
    tempNMOCS = squeeze(sum(guessNTrialsFromDataMOCS(:,:,:)));
    tempNQUEST = squeeze(sum(guessNTrialsFromDataQUEST(:,:,:)));
    fprintf(rateTablefile,'%s\t%d/%d\t%d/%d\t%d/%d\t%d/%d\t%d/%d (%0.1f%%)\n', ...
            'Total', ...
            tempTempMOCS(1,ss), tempNMOCS(1,ss), ...
            tempTempQUEST(1,ss), tempNQUEST(1,ss), ...
            tempTempMOCS(2,ss), tempNMOCS(2,ss), ...
            tempTempQUEST(2,ss), tempNQUEST(2,ss), ...
            sum(tempTempMOCS(:,ss)+tempTempQUEST(:,ss)), ...
            sum(tempNMOCS(:,ss)+tempNQUEST(:,ss)), ...
            round(100*sum(tempTempMOCS(:,ss)+tempTempQUEST(:,ss)) / ...
            sum(tempNMOCS(:,ss)+tempNQUEST(:,ss)),1) ...
            );
    fprintf(rateTablefile,'(n = %d)\t(%0.1f%%)\t(%0.1f%%)\t(%0.1f%%)\t(%0.1f%%)\n',length(theSubjects), ...
        round(100*tempTempMOCS(1,ss)/tempNMOCS(1,ss),1), ...
        round(100*tempTempQUEST(1,ss)/tempNQUEST(1,ss),1), ...
        round(100*tempTempMOCS(2,ss)/tempNMOCS(2,ss),1), ...
        round(100*tempTempQUEST(2,ss)/tempNQUEST(2,ss),1));    
end

% Write out lapse table values in form they can go into Excel table
% Not enough QUEST trials for QUEST estimates to be meaningful.
tempMOCS = round(lapseFromDataMOCS.*lapseNTrialsFromDataMOCS);
for ss = 1:size(tempMOCS,2)
    fprintf(rateTablefile,'\nLapse Rates, Session %d\n',ss);
    for pp = 1:length(theSubjects)
        fprintf(rateTablefile,'%s\t%d/%d\t%d/%d\t%d/%d (%0.1f%%)\n', ...
            theSubjects{pp}, ...
            tempMOCS(pp,1,ss), lapseNTrialsFromDataMOCS(pp,1,ss), ...
            tempMOCS(pp,2,ss), lapseNTrialsFromDataMOCS(pp,2,ss), ...
            nansum(tempMOCS(pp,:,ss)), ...
            nansum(lapseNTrialsFromDataMOCS(pp,:,ss)), ...
            100*sum(tempMOCS(pp,:,ss)) / ...
            sum(lapseNTrialsFromDataMOCS(pp,:,ss)) ...
            );
    end
    tempTempMOCS = squeeze(sum(tempMOCS(:,:,:)));
    tempNMOCS = squeeze(sum(lapseNTrialsFromDataMOCS(:,:,:)));
    fprintf(rateTablefile,'%s\t%d/%d\t%d/%d\t%d/%d (%0.1f%%)\n', ...
            'Total', ...
            tempTempMOCS(1,ss), tempNMOCS(1,ss), ...
            tempTempMOCS(2,ss), tempNMOCS(2,ss), ...
            nansum(tempTempMOCS(:,ss)), ...
            nansum(tempNMOCS(:,ss)), ...
            round(100*sum(tempTempMOCS(:,ss)) / ...
            sum(tempNMOCS(:,ss)),1));
    fprintf(rateTablefile,'(n = %d)\t(%0.1f%%)\t(%0.1f%%)\n',length(theSubjects), ...
        round(100*tempTempMOCS(1,ss)/tempNMOCS(1,ss),1), ...
        round(100*tempTempMOCS(2,ss)/tempNMOCS(2,ss),1));    
end
fclose(rateTablefile);
        
