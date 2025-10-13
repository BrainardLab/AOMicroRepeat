%% FigureSetup
%
% Set some parameters that we want to keep the same across figures and read in the output
% of the analysis.  Saves redoing this for each figure script.

%% Parameters
plotFont = 'Times';
axisFontSize = 14;
labelFontSize = 18;
titleFontSize = 20;
legendFontSize = 14;
markerSize = 12;
plotSize = [100, 100,800,800];

%% Path to analysis output dir
analysisDir = getpref('New_analysis_20250912','analysisDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\New_analysis_20250912\analysisDir');

%% Read table of fit information
wState = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet'); %,'VariableNamingRule','preserve');
warning(wState);