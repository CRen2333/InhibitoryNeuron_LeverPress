function [dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign, fPath] = delayDecRecordings
%overview over recordings that are used in the delayed detection study.
%Also includes some other basic variables to ensure that different codes
%will use the same values.

dataOverview = {
    
'mSM30' 'Visual' '10-Oct-2017'; ...
'mSM30' 'Visual' '12-Oct-2017'; ...
'mSM34' 'Visual' '01-Dec-2017'; ...
'mSM34' 'Visual' '02-Dec-2017'; ...
'mSM36' 'Visual' '05-Dec-2017'; ...
'mSM36' 'Visual' '07-Dec-2017'; ...
'mSM46' 'Visual' '01-Dec-2017'; ...
'mSM46' 'Visual' '13-Jan-2018'; ...
'mSM49' 'Visual' '19-Dec-2017'; ...
'mSM49' 'Visual' '12-Mar-2018'; ...
'mSM57' 'Visual' '02-Feb-2018'; ...
'mSM57' 'Visual' '08-Feb-2018'; ...

'mSM43' 'Audio' '21-Nov-2017'; ...
'mSM43' 'Audio' '23-Nov-2017'; ...
'mSM44' 'Audio' '21-Nov-2017'; ...
'mSM44' 'Audio' '29-Nov-2017'; ...
'mSM53' 'Audio' '14-Mar-2018'; ...
'mSM53' 'Audio' '21-Mar-2018'; ...
'mSM55' 'Audio' '13-Feb-2018'; ...
'mSM55' 'Audio' '16-Feb-2018'; ...
'mSM56' 'Audio' '22-Feb-2018'; ...
'mSM56' 'Audio' '27-Feb-2018'; ...
'mSM65' 'Audio' '08-Sep-2018'; ...
'mSM65' 'Audio' '10-Sep-2018'; ...
'mSM66' 'Audio' '08-Sep-2018'; ...
'mSM66' 'Audio' '10-Sep-2018'; ...

};

%regressors for motor-based reconstruction
motorLabels = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'fastPupil' 'slowPupil'... 
    'piezo' 'whisk' 'nose' 'face' 'body' 'BaselineMove' 'HandleMove' 'StimulusMove' 'WaitMove' 'Move' 'bhvVideo'}; 

%regressors for sensory-based reconstruction
sensorLabels = {'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'}; 

%regressors for cognitive reconstruction
cogLabels = {'time' 'Choice' 'reward' 'prevReward' 'prevChoice' 'prevMod' 'water'}; 

%segment indices to compute changes per task episodes.
segIdx = {1:54 55:81 89:99 122:132 140:162 170:189}; %segment index for handel aligned data
segLabels = {'Baseline' 'Handle' 'Stim1' 'Stim2' 'Wait' 'Response'};
segIdxRealign = {1:54 55:79 80:98 113:131 132:162 163:179}; %segment index after realignment on both handle and stimulus

% file path on grid server
fPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
