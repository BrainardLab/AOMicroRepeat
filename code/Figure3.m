

%% Path to analysis output dir
analysisDir = getpref('AOMicroRepeat','analysisDir');

%% Read table of fit information
dataTable = readtable(fullfile(analysisDir,'fitTable.xlsx'),'FileType','spreadsheet','VariableNamingRule','preserve');

%% Get all MOCS data
theParticipants = unique(dataTable.Subject);
theDiameters = unique(dataTable.Diameter);
for pp = 1:length(theSubjects)
    for dd = 1:length(theDiameters)
        index = strcmp(dataTable.Subject,theSubjects{pp}) & (dataTable.Diameter == theDiameters(dd)) & strcmp(dataTable.Method,'MOCS') & strcmp(dataTable.Split,'All');
        dataTable(index,:).Criterion
    end
end
