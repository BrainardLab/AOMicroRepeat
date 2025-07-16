
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

%% Path to analysis output dir
analysisDir = getpref('AOMicroRepeat','analysisDir');

%% Read table of fit information
wState = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet'); %,'VariableNamingRule','preserve');
warning(wState);

%% Get all MOCS data for full sessions
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'MOCS') & strcmp(dataTable.Split,'All'));
            sensitivityMOCS(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
        end
    end
end

%% Get all QUEST data for full sessions
theSubjects = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
theSessions = unique(dataTable.Session);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        for ss = 1:length(theSessions)
            index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & (dataTable.Session == theSessions(ss) & ...
                strcmp(dataTable.Method,'QUEST') & strcmp(dataTable.Split,'All'));
            sensitivityQUEST(pp,dd,ss) = -dataTable.CorrectedThreshold_dB_(index);
        end
    end
end

%% Make Figure 3a
limMin = 18; limMax = 30;
f = figure; clf; hold on
set(gca,'FontName', 'Helvetica','FontSize',axisFontSize);
plot(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1),'bo','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2),'b^','MarkerFaceColor','b','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1),'ro','MarkerFaceColor','r','MarkerSize',markerSize);
plot(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2),'r^','MarkerFaceColor','r','MarkerSize',markerSize);
plot([limMin limMax],[limMin limMax],'k:','LineWidth',1);
xlabel('MOCS sensitivity (dB)','FontSize',labelFontSize);
ylabel('QUEST sensitivity (dB)','FontSize',labelFontSize);
legend( ...
    {sprintf('Session %d, %d arcmin',theSessions(1),theDiameters(1)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(2),theDiameters(1)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(1),theDiameters(2)) ; ...
    sprintf('Session %d, %d arcmin',theSessions(2),theDiameters(2)) ; ...
    ''}, ...
    'Location','SouthEast');
axis('square');
axis([limMin limMax limMin limMax]);
saveas(gcf,fullfile(analysisDir,outputVariant,'Figure3a.pdf'),'pdf');

%% t-tests
[~,p(1,1)] = ttest(sensitivityMOCS(:,1,1),sensitivityQUEST(:,1,1));
[~,p(1,1)] = ttest(sensitivityMOCS(:,1,2),sensitivityQUEST(:,1,2));
[~,p(2,1)] = ttest(sensitivityMOCS(:,2,1),sensitivityQUEST(:,2,1));
[~,p(2,2)] = ttest(sensitivityMOCS(:,2,2),sensitivityQUEST(:,2,2));
fprintf('MOCS vs QUEST t-test p values\n');
for dd = 1:length(theDiameters)
    for ss = 1:length(theSessions)
        fprintf('\t%d arcmin, session %d, p = %0.3f\n',theDiameters(dd),theSessions(ss),p(dd,ss));
    end
end
