function AOCompObserverLocalHook
% AOCompObserverLocalHook
%
% Configure things for working on the AOCompObserver project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by defalut,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUseProject('AOCompObserver') to set up for
% this project.  You then edit your local copy to match your configuration.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

%% Say hello
theProject = 'AOCompObserver';
fprintf('Running %s local hook\n',theProject);

%% Remove old preferences
if (ispref(theProject))
    rmpref(theProject);
end

%% Put project toolbox onto path
%
% Specify project name and location
projectName = theProject;
projectBaseDir = tbLocateProject(theProject);

%% Figure out where baseDir for other kinds of data files is.
sysInfo = GetComputerInfo();
switch (sysInfo.localHostName)
    case 'eagleray'
        % DHB's desktop
        baseDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
 
    otherwise
        % Some unspecified machine, try user specific customization
        switch(sysInfo.userShortName)
            % Could put user specific things in, but at the moment generic
            % is good enough.
            otherwise
                baseDir = fullfile('/Users',sysInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
        end
end

%% Set preferences for project i/o

% This is where the psychophysical data is stored
dataDir = fullfile(baseDir,'AOPY_data',theProject);

% This is where the output goes
analysisDir = fullfile(baseDir,'AOPY_analysis',theProject);

% Set the preferences
setpref(theProject,'dataDir',dataDir);
setpref(theProject,'analysisDir',analysisDir);
