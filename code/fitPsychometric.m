% fitPsychometric
%
% Wrapper into Palamedes that fits the logistic with parameters reasonable
% for our data.
function [fit_log_intensity,fit_psychometric,corrected_fit_psychometric, ...
    log_threshold,corrected_log_threshold,psiParamsValues] = fitPsychometric(log_intensity,num_pos,out_of_num,threshold_criterion,convert_to_db, ...
    guessUpper,lapseUpper,slopeLower,slopeUpper)
%
% Usage:
%   [fit_log_intensity,fit_psychometric,corrected_fit_psychonmetric, ...
%        log_threshold,corrected_log_threshold,psiParamsValues] = fitPsychometric(log_intensity,num_pos,out_of_num,threshold_criterion)
%
% Inputs:
%   log_intensity : (vector) log10 stimulus levels tested
%   num_pos : (vector) for each stimulus level tested, the number of trials a positive response was given
%   out_of_num : (vector) for each stimulus level tested, the total number of trials
%   threshold_criterion : [scalar] criterion for threshold
%
% History:
%   07/05/21  amn  Wrote it.
%   06/28/25  dhb  Adopted for this project

% Calculate x-axis values to plot
%
% Stimulus values to plot for a smooth fit.
if (convert_to_db)
    fit_log_intensity = linspace(-35,0,1000);
else
    fit_log_intensity = linspace(-3.5,0,1000);
end

% Psychometric function form (alternative: PAL_Gumbel).
PF = @PAL_Logistic;

% 'paramsFree' is a boolean vector that determines what parameters get
% searched over (1: free parameter, 0: fixed parameter).
paramsFree = [1 1 1 1];

% Set up starting points:
if (convert_to_db)
    searchGrid.alpha = -35:.1:0;
    searchGrid.beta = slopeLower:0.1:slopeUpper;
else
    searchGrid.alpha = -3.5:.01:0;
    searchGrid.beta = slopeLower:1:slopeUpper;
end
searchGrid.gamma = linspace(0,guessUpper,3);
searchGrid.lambda = linspace(0,lapseUpper,6);

% Set up lapse limits.
guessLimits = [0 guessUpper];
lapseLimits = [0 lapseUpper];

% Set up standard options for Palamedes search.
options = PAL_minimize('options');

% Fit with Palemedes Toolbox.
[psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
    'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
    'checkLimits',false);

% Check whether beta is bigger or smaller than we would like,
% constrain it if so.
if (psiParamsValues(2) > max(searchGrid.beta))
    searchGrid.beta = max(searchGrid.beta);
    paramsFree = [1 0 1 1];
    [psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
        'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
        'checkLimits',false);
    % disp('Overrun on beta')
    % psiParamsValues
elseif (psiParamsValues(2) < min(searchGrid.beta))
    searchGrid.beta = min(searchGrid.beta);
    paramsFree = [1 0 1 1];
    [psiParamsValues] = PAL_PFML_Fit(log_intensity,num_pos,out_of_num,searchGrid,paramsFree,PF, ...
        'lapseLimits',lapseLimits,'guessLimits',guessLimits,'searchOptions',options,'gammaEQlambda',false, ...
        'checkLimits',false);
    % disp('Underunrun on beta')
    % psiParamsValues
end

% Get fitted curve values.
fit_psychometric = PF(psiParamsValues,fit_log_intensity);

% Calculate threshold: the difference between the stimulus levels at
% performances equal to 0.7602 and 0.5.
log_threshold = PF(psiParamsValues,threshold_criterion,'inverse');

% Now correct for guessing and lapse rate
psiParamsValuesCorrect = psiParamsValues;
psiParamsValuesCorrect(3) = 0;
psiParamsValuesCorrect(4) = 0;
corrected_fit_psychometric = PF(psiParamsValuesCorrect,fit_log_intensity);
corrected_log_threshold = PF(psiParamsValuesCorrect,threshold_criterion,'inverse');

end