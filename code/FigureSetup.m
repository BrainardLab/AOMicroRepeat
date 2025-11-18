%% FigureSetup
%
% Set some parameters that we want to keep the same across figures and read in the output
% of the analysis.  Saves redoing this for each figure script.

%% Parameters
plotFont = 'Times';
axisFontSize = 12;
labelFontSize = 12;
titleFontSize = 20;
legendFontSize = 12;
markerSize = 5;
% plotSize = [100, 100,600,600];
limMin = 16; limMax = 32;

%% Path to analysis output dir

 analysisDir = getpref('New_analysis_20250912','analysisDir');

% analysisDir = getpref('AOMicroRepeat','analysisDir');
% setpref('AOMicroRepeat','setpref','dataDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\New_analysis_20250912\dataDir');
setpref('New_analysis_20250912','analysisDir','C:\Users\niveg\Aguirre-Brainard Lab Dropbox\Nivedhitha Govindasamy\AO-Microperimetry Repeatability Paper\Data_for_paper\David_code_analysis\New_analysis_20250912\analysisDir');
analysisDir = getpref('New_analysis_20250912','analysisDir');

%% Read table of fit information
wState = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
dataTable = readtable(fullfile(analysisDir,outputVariant,'fitTable.xlsx'),'FileType','spreadsheet'); %,'VariableNamingRule','preserve');
warning(wState);