%% Run Function
% Main function that is called by external wrapper run files for each task.
% Inputs: pt (runs on one patient based on subject ID) and task (Resting,
% CATDI, AmpMod, AmpModControl, Ablation). Make sure if running Amp Mod,
% stim is set to 1.

function runfunction(pt,task,stim)

switch nargin
    case 2
        stim = 0;
    case 3
        stim = stim;
end

switch stim
    case 0
        Step01_Preprocessing;
        %
        % clearvars -except pt task stim
        redopoweranalysis = 1;
        Step02_ChannelSelection;
        Step02_ROI_Activity_Reconstructions;
        Step02_PowerAnalysis;
        clearvars -except pt task
        Step03_Analysis
        Step04_SharedModel_and_Temperature;
        Step04_AperiodicComp;
        Step05_GatherParams;
        clearvars -except pt task

    case 1  
        Step01_PreprocessingStim;
        %
        clearvars -except pt task
        redopoweranalysis = 1;
        Step02_ROI_Activity_ReconstructionsStim;
        Step02_PowerAnalysisStim;
        clearvars -except pt task
        Step03_AnalysisStim
        Step04_SharedModel_and_TemperatureStim;
        Step04_AperiodicCompStim;
        Step05_GatherParamsStim;
        clearvars -except pt task
end

end
