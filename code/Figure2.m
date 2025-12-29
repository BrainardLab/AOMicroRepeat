% Copy .png figures for Figure2 over to the top level.

%% Path to analysis output dir
 analysisDir = getpref('AOMicroRepeat','analysisDir');
 outputVariant = 'SlopeFree1';

 %% QUEST
copyfile(fullfile(analysisDir,outputVariant,'11002','Session2','Size8','QUEST','All','staircasePlotNoTitle_run1.png'),fullfile(analysisDir,outputVariant,'figure2a.png'));
copyfile(fullfile(analysisDir,outputVariant,'11002','Session2','Size8','QUEST','All','psychometricFcn_NoTitle.png'),fullfile(analysisDir,outputVariant,'figure2b.png'));

%% MOCS 
copyfile(fullfile(analysisDir,outputVariant,'11002','Session2','Size8','MOCS','All','staircasePlotNoTitle_run1.png'),fullfile(analysisDir,outputVariant,'figure2c.png'));
copyfile(fullfile(analysisDir,outputVariant,'11002','Session2','Size8','MOCS','All','psychometricFcn_NoTitle.png'),fullfile(analysisDir,outputVariant,'figure2d.png'));