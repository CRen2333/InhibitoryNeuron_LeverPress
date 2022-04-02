%% CR_4042832-O/L/R
clear
IN = 'VIP_SOM';
Animal = 'CR_4042832-O';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200821'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200821'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200822'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200822'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4113794-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4113794-L';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'201003'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'201003'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'201004'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'201004'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153298-R/LR
clear
IN = 'VIP_SOM';
Animal = 'CR_4153298-LR';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'201025'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'201025'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'201026'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'201026'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3974574-R/3974573-O
clear
IN = 'VIP_SOM';
Animal = 'CR_3974573-O';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200404'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200404'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200405'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200405'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-L/R
clear
IN = 'VIP_SOM';
Animal = 'CR_4017421-R';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200606'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200606'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200607'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200607'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4041831-R/LR
clear
IN = 'VIP_SOM';
Animal = 'CR_4042831-LR';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200706'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200706'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200707'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200707'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4041830-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4042830-L';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200718'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200718'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200719'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200719'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017386-O
clear
IN = 'VIP_SOM';
Animal = 'CR_4017386-O';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200413'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200413'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200414'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200414'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017386-L
clear
IN = 'VIP_SOM';
Animal = 'CR_4017386-L';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200523'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200523'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200524'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200524'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-O
clear
IN = 'VIP_SOM';
Animal = 'CR_4017421-O';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200527'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200527'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200528'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200528'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-LR
clear
IN = 'VIP_SOM';
Animal = 'CR_4017421-LR';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200615'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200615'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200616'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200616'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042831-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4042831-L';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200620'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200620'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200621'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200621'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042830-R
clear
IN = 'VIP_SOM';
Animal = 'CR_4042830-R';
Stage = 'Nai';
Imaging_Fields{1}.Date = {'200720'};
Imaging_Fields{1}.Piezo = ones(1,1);
Imaging_Fields{2}.Date = {'200720'};
Imaging_Fields{2}.Piezo = ones(1,1).*2;
Imaging_Fields{3}.Date = {'200721'};
Imaging_Fields{3}.Piezo = ones(1,1);
Imaging_Fields{4}.Date = {'200721'};
Imaging_Fields{4}.Piezo = ones(1,1).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4113794-R
clear
IN = 'VIP_SOM';
Animal = 'CR_4113794-R';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201023','201025','201027'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201023','201025','201027'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201024','201026','201028'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201024','201026','201028'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153298-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4153298-L';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201106','201108','201110'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201106','201108','201110'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201107','201109','201111'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201107','201109','201111'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153299-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4153299-L';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201120','201122','201124'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201120','201122','201124'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201121','201123','201125'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201121','201123','201125'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153369-O
clear
IN = 'VIP_SOM';
Animal = 'CR_4153369-O';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201128','201130','201202'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201128','201130','201202'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201129','201201','201203'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201129','201201','201203'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153382-O
clear
IN = 'VIP_SOM';
Animal = 'CR_4153382-O';
Stage = 'Exp';
% Imaging_Fields{1}.Date = {'201213','201215','201217'};
% Imaging_Fields{1}.Piezo = ones(1,3);
% Imaging_Fields{2}.Date = {'201213','201215','201217'};
% Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{1}.Date = {'201214','201216','201218'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201214','201216','201218'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153382-L
clear
IN = 'VIP_SOM';
Animal = 'CR_4153382-L';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201213','201215','201217'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201213','201215','201217'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201214','201216','201218'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201214','201216','201218'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4113794-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4113794-L';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201019','201021','201023'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201019','201021','201023'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201020','201022','201024'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201020','201022','201024'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153298-R/LR
clear
IN = 'VIP_SOM';
Animal = 'CR_4153298-LR';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201110','201112','201114'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201110','201112','201114'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201111','201113','201115'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201111','201113','201115'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153299-R
clear
IN = 'VIP_SOM';
Animal = 'CR_4153299-R';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201120','201122','201124'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201120','201122','201124'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201121','201123','201125'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201121','201123','201125'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153369-L
clear
IN = 'VIP_SOM';
Animal = 'CR_4153369-L';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201128','201130','201202'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201128','201130','201202'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201129','201201','201203'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201129','201201','201203'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4153382-R
clear
IN = 'VIP_SOM';
Animal = 'CR_4153382-R';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'201213','201215','201217'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'201213','201215','201217'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'201214','201216','201218'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'201214','201216','201218'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4259761-O/L
clear
IN = 'VIP_SOM';
Animal = 'CR_4259761-O';
% Stage = 'Nai';
% Imaging_Fields{1}.Date = {'210409'};
% Imaging_Fields{1}.Piezo = ones(1,1);
% Imaging_Fields{2}.Date = {'210409'};
% Imaging_Fields{2}.Piezo = ones(1,1).*2;
% Imaging_Fields{3}.Date = {'210410'};
% Imaging_Fields{3}.Piezo = ones(1,1);
% Imaging_Fields{4}.Date = {'210410'};
% Imaging_Fields{4}.Piezo = ones(1,1).*2;

Stage = 'Exp';
Imaging_Fields{1}.Date = {'210425','210427','210429'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'210425','210427','210429'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'210426','210428','210430'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'210426','210428','210430'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4259761-R/LR
clear
IN = 'VIP_SOM';
Animal = 'CR_4259761-LR';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'210429','210501','210503'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'210429','210501','210503'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'210430','210502','210504'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'210430','210502','210504'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4259762-O/L/LL
clear
IN = 'VIP_SOM';
Animal = 'CR_4259762-LL';
Stage = 'Exp';
Imaging_Fields{1}.Date = {'210517','210519','210521'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{2}.Date = {'210517','210519','210521'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{3}.Date = {'210518','210520','210522'};
Imaging_Fields{3}.Piezo = ones(1,3);
Imaging_Fields{4}.Date = {'210518','210520','210522'};
Imaging_Fields{4}.Piezo = ones(1,3).*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep Stage filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');