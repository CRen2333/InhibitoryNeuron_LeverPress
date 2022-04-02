%% MOM
%% 1) Copy bpod to server people folder
Animals = {'CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'};
for curr_animal = 1:length(Animals)
    if ismember(Animals{curr_animal},{'CR_3526643-R','WL_3526641-R','WL_3526642-L','WL_3526642-R','CR_3702608-O','CR_3702608-LR','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3887044-O','CR_3887044-L','CR_3786160-LR','CR_3887043-O','CR_3887043-L','CR_3887041-O','CR_4259757-R','CR_4302983-LL','CR_4302983-R','CR_4302984-O','CR_4302984-R','CR_4302984-LR','CR_4303048-L','CR_4303048-R','CR_4303048-LR','CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{'WL_3526549-R','WL_3526580-L','WL_3526580-O','WL_3526578-O','WL_3526578-R','WL_3547272-O','WL_3547272-L','WL_3547273-R','WL_3547272-LR','CR_3672035-L','CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O','CR_3786142-O','CR_3786142-L','CR_3786142-R','CR_3886982-O','CR_3886982-L','CR_3886982-LR','CR_3886982-R','CR_3887040-O','CR_3887040-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3561044-O','CR_3561044-L'})
        IN = 'ChAT';
    elseif ismember(Animals{curr_animal},{'CR_3575263-L','CR_3575263-R'})
        IN = 'DBH';
    elseif ismember(Animals{curr_animal},{'CR_3672031-L','CR_3786160-LL','CR_3786159-O'})
        IN = 'VIP_ChAT_Cas';
    elseif ismember(Animals{curr_animal},{'CR_3974574-O','CR_3974574-L','CR_3974574-R','CR_3974573-O','CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-L','CR_4017421-R','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042830-R','CR_4042832-O','CR_4042832-L','CR_4042832-R','CR_4113794-O','CR_4113794-L','CR_4113794-R','CR_4153298-O','CR_4153298-L','CR_4153298-R','CR_4153298-LR',...
            'CR_4153299-O','CR_4153299-L','CR_4153369-O','CR_4153382-O','CR_4153382-L','CR_4153299-R','CR_4153369-L','CR_4153382-R','CR_4153382-LR','CR_4259761-O','CR_4259761-L','CR_4259761-R','CR_4259761-LR','CR_4259762-O','CR_4259762-L','CR_4259762-LL'})
        IN = 'VIP_SOM';
    elseif ismember(Animals{curr_animal},{'CR_4113793-O','CR_4113793-L','CR_4153383-O','CR_4153383-L','CR_4153383-R','CR_4153383-LR'})
        IN = 'VIP_hM4Di';
    end
    DataPath = fullfile('Z:\Data\ImagingRig1\Behavior',Animals{curr_animal},'Lever\Session Data');
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'Bpod');
    if ~exist(TargetPath)
        mkdir(TargetPath);
    end
    copyfile(DataPath,TargetPath);
end
clear all

%% 2) Read Ephus trace
Animals = {'CR_4365807-O'};
Dates = {'210530','210531'};
for curr_animal = 1:length(Animals)
    if ismember(Animals{curr_animal},{'CR_3526643-R','WL_3526641-R','WL_3526642-L','WL_3526642-R','CR_3702608-O','CR_3702608-LR','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3887044-O','CR_3887044-L','CR_3786160-LR','CR_3887043-O','CR_3887043-L','CR_3887041-O','CR_4259757-R','CR_4302983-LL','CR_4302983-R','CR_4302984-O','CR_4302984-R','CR_4302984-LR','CR_4303048-L','CR_4303048-R','CR_4303048-LR','CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{'WL_3526549-R','WL_3526580-L','WL_3526580-O','WL_3526578-O','WL_3526578-R','WL_3547272-O','WL_3547272-L','WL_3547273-R','WL_3547272-LR','CR_3672035-L','CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O','CR_3786142-O','CR_3786142-L','CR_3786142-R','CR_3886982-O','CR_3886982-L','CR_3886982-LR','CR_3886982-R','CR_3887040-O','CR_3887040-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3575263-L','CR_3575263-R'})
        IN = 'DBH';
    elseif ismember(Animals{curr_animal},{'CR_3672031-L','CR_3786160-LL','CR_3786159-O'})
        IN = 'VIP_ChAT_Cas';
    elseif ismember(Animals{curr_animal},{'CR_3974574-O','CR_3974574-L','CR_3974574-R','CR_3974573-O','CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-L','CR_4017421-R','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042830-R','CR_4042832-O','CR_4042832-L','CR_4042832-R','CR_4113794-O','CR_4113794-L','CR_4113794-R','CR_4153298-O','CR_4153298-L','CR_4153298-R','CR_4153298-LR',...
            'CR_4153299-O','CR_4153299-L','CR_4153369-O','CR_4153382-O','CR_4153382-L','CR_4153299-R','CR_4153369-L','CR_4153382-R','CR_4153382-LR','CR_4259761-O','CR_4259761-L','CR_4259761-R','CR_4259761-LR','CR_4259762-O','CR_4259762-L','CR_4259762-LL'})
        IN = 'VIP_SOM';
    elseif ismember(Animals{curr_animal},{'CR_4113793-O','CR_4113793-L','CR_4153383-O','CR_4153383-L','CR_4153383-R','CR_4153383-LR'})
        IN = 'VIP_hM4Di';
    end
    for curr_date = 1:length(Dates)
        Sessions = dir(fullfile('Z:\Data\ImagingRig1',Dates{curr_date},Animals{curr_animal}));
        Sessions = {Sessions(cellfun(@(x) ~isempty(strfind(x, 'Rec')), {Sessions.name})).name};
        for curr_session = 1:length(Sessions)
            DataPath = fullfile('Z:\Data\ImagingRig1',Dates{curr_date},Animals{curr_animal},Sessions{curr_session},[Animals{curr_animal}(1:2) '0001']);
            TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace',Dates{curr_date});
            if ~exist(TargetPath)
                mkdir(TargetPath);
            end
            data = MomEphus_load_xsg_continuous(DataPath);
            if isempty(data.channels)
                warning('No xsg files saved');
                No_xsg = true;
                continue
            else
                No_xsg = false;
            end
            xsg_sample_rate = 10000;
            for curr_channel = 1:length(data.channel_names)
                channel = data.channel_names{curr_channel};
                channel_index = cellfun(@(x) strcmp(x,channel),data.channel_names);
                temp_data = data.channels(:,channel_index);
                if sum(isnan(temp_data))>0
                    disp('Damaged data');
                    temp_data = temp_data(~isnan(temp_data));
                end
                switch channel
                    case 'Frame' % 2p imaging
                        % Get frame times (in seconds) from frame trigger trace                        
                        Frame_times{curr_session} = (find(temp_data(2:end) > 2.5 & temp_data(1:end-1) < 2.5) + 1)/xsg_sample_rate;
                        % check
                        diff_Frame_times = diff(Frame_times{curr_session});
                        jump_index = find(diff_Frame_times>0.3);
                        if ~isempty(jump_index)
                            if length(jump_index) == 1
                                warning('Correcting Frame_times');
                                Frame_times{curr_session} = Frame_times{curr_session}(jump_index+1:end);
                            else
                                warning('Check frame times!!!');
%                                 break
                            end
                        end
                        clear temp_data diff_Frame_times jump_index
                    case 'Lever'
                        [Lever_active{curr_session}, Lever_force_resample{curr_session}, lever_force_smooth, lever_velocity_envelope_smooth,]...
                            = CR_parseLeverMovement_Updated(temp_data, 500, 150); % 1K Hz
                        clear temp_data lever_force_smooth lever_velocity_envelope_smooth
                    case 'Lick'
                        [Licking_active{curr_session}, Licking_bout{curr_session}, Licking_bout_switch{curr_session}] = CR_GetLickingState(temp_data,3);
                        clear temp_data
                    case 'PupilFrame'
                        PupilFrame_times{curr_session} = (find(temp_data(2:end) > 2.5 & temp_data(1:end-1) < 2.5) + 1)/xsg_sample_rate;
                        clear temp_data
                    case 'Trial_number'
                        [TrialInfo_xsg{curr_session}] = CR_ReadBitCodeShift(temp_data, 0);
                        clear temp_data
                end
            end
        clear data
        end
        if No_xsg
            continue
        end
        saveFile = [TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'];
        save(saveFile,'Frame_times','Lever_active','Lever_force_resample','Licking_active','Licking_bout','Licking_bout_switch','PupilFrame_times','TrialInfo_xsg','-v7.3');
        clear Frame_times Lever_active Lever_force_resample Licking_active Licking_bout Licking_bout_switch PupilFrame_times TrialInfo_xsg
    end
end
clear all

%% check frame times
Animals = {'CR_4302984-O','CR_4302984-R','CR_4302984-LR'};
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    if ismember(Animals{curr_animal},{'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','CR_3526643-R','WL_3526641-R','WL_3526642-L','WL_3526642-R','CR_3702608-O','CR_3702608-LR','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3887044-O','CR_3887044-L','CR_3786160-LR','CR_3887043-O','CR_3887043-L','CR_3887041-O','CR_4259757-R','CR_4302983-LL','CR_4302983-R','CR_4302984-O','CR_4302984-R','CR_4302984-LR'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{'KP_3459921_1','KP_3461990_1','WL_3526549-R','WL_3526580-L','WL_3526580-O','WL_3526578-O','WL_3526578-R','WL_3547272-O','WL_3547272-L','WL_3547273-R','WL_3547272-LR','CR_3672035-L','CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O','CR_3786142-O','CR_3786142-L','CR_3786142-R','CR_3886982-O','CR_3886982-L','CR_3886982-LR','CR_3886982-R','CR_3887040-O','CR_3887040-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3575263-L','CR_3575263-R'})
        IN = 'DBH';
    elseif ismember(Animals{curr_animal},{'CR_3672031-L','CR_3786160-LL','CR_3786159-O'})
        IN = 'VIP_ChAT_Cas';
    elseif ismember(Animals{curr_animal},{'CR_3974574-O','CR_3974574-L','CR_3974574-R','CR_3974573-O','CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-L','CR_4017421-R','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042830-R','CR_4042832-O','CR_4042832-L','CR_4042832-R','CR_4113794-O','CR_4113794-L','CR_4113794-R','CR_4153298-O','CR_4153298-L','CR_4153298-R','CR_4153298-LR',...
            'CR_4153299-O','CR_4153299-L','CR_4153369-O','CR_4153382-O','CR_4153382-L','CR_4153299-R','CR_4153369-L','CR_4153382-R','CR_4153382-LR'})
        IN = 'VIP_SOM';
    elseif ismember(Animals{curr_animal},{'CR_4113793-O','CR_4113793-L','CR_4153383-O','CR_4153383-L','CR_4153383-R','CR_4153383-LR'})
        IN = 'VIP_hM4Di';
    end
    Dates = dir(fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace'));
    Dates = {Dates.name};
    Dates = Dates(3:end);
    for curr_date = 1:length(Dates)
        TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace',Dates{curr_date});
        if ~exist([TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'])
            continue
        end
        load([TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'],'Frame_times','-mat');
        for curr_session = 1:length(Frame_times)
            diff_Frame_times = diff(Frame_times{curr_session});
            jump_index = find(diff_Frame_times>0.3);
            if ~isempty(jump_index)                           
                disp([Animal ', ' Dates{curr_date} ', '  num2str(curr_session)]);
            end
            clear temp_data diff_Frame_times jump_index
        end
    end
end                    

%% check frame times, Pharm
IN = 'SOM';
% SessionType = 'Nai_Ant';
SessionType = 'Exp_Ago';
switch SessionType
    case 'Nai_Ant'
        Animals = {'CR_3702608-O','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3887044-L','CR_3786160-LR','CR_3887043-O','CR_3887043-L','CR_3887041-O'};
        Drug = 'Ant';
    case 'Exp_Ago'
        Animals = {'WL_3526642-L','WL_3526642-R','CR_3619106-R','CR_3619106-LR','CR_3633192-L','CR_3633192-R','CR_3702608-O','CR_3702608-LR','CR_3658844-R','CR_3672031-O','CR_3633193-L'};           
        Drug = 'Ago';
end

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    switch SessionType
        case 'Exp_Ago'
            if ismember(Animal,{'CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_3633193-L'})
                SessionType_2 = 'Exp_Ago_Ant';
            else
                SessionType_2 = SessionType;
            end
        otherwise
            SessionType_2 = SessionType;
    end
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'Pharm',SessionType_2,'LeverTrace');
    Dates = dir(DataPath);
    Dates = {Dates.name};
    Dates = Dates(3:end);
    for curr_date = 1:length(Dates)
        TargetPath = fullfile(DataPath,Dates{curr_date});
        load([TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'],'Frame_times','-mat');
        for curr_session = 1:length(Frame_times)
            diff_Frame_times = diff(Frame_times{curr_session});
            jump_index = find(diff_Frame_times>0.3);
            if ~isempty(jump_index)                           
                disp([Animal ', ' Dates{curr_date} ', '  num2str(curr_session)]);
            end
            clear temp_data diff_Frame_times jump_index
        end
    end
end 

%% Bscope1 
%% 1) Copy bpod to people folder, BScope1
Animals = {'CR_3672020-R'};
Dates = {'190511','190512','190513'};

for curr_animal = 1:length(Animals)
    if ismember(Animals{curr_animal},{'CR_3619106-R','CR_3619106-LR','CR_3633192-R','CR_3633192-L','CR_3702608-O','CR_3702608-LR','CR_3633193-L','CR_3633193-R','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3619073-L','CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R'})
        IN = 'ChAT';
    elseif ismember(Animals(curr_animal),{'CR_3575265-LR'})
        IN = 'DBH';
    elseif ismember(Animals{curr_animal},{'CR_3672031-O','CR_3786160-LL','CR_3786159-O'})
        IN = 'VIP_ChAT_Cas';
    end
    for ii = 1:length(Dates)
        DataPath = fullfile('Z:\Data\ImagingRig3\',Dates{ii},Animals{curr_animal},'Rec_1',['CR_' Dates{ii} '_' Animals{curr_animal}(4:end) '_Rec1.mat']);
        TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'Bpod');
        if ~exist(TargetPath)
            mkdir(TargetPath);
        end
        copyfile(DataPath,TargetPath);
    end
end
clear all

%% 2) Read Ephus trace, BScope1
Animals = {'CR_3672020-R'};
Dates = {'190511','190512','190513'};

active_threshold = 1; % CR_4383143-R, CR_4429262-O
for curr_animal = 1:length(Animals)
    if ismember(Animals{curr_animal},{'CR_3619106-R','CR_3619106-LR','CR_3633192-R','CR_3633192-L','CR_3633193-R','CR_3633193-L','CR_3702608-O','CR_3702608-LR','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3619073-L','CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R','CR_4383143-L','CR_4383143-R','CR_4429262-O','CR_4383144-R','CR_4412582-R','CR_4412583-LR'})
        IN = 'ChAT';
    elseif ismember(Animals{curr_animal},{'CR_3575265-O','CR_3575265-LR'})
        IN = 'DBH';
    end
    for curr_date = 1:length(Dates)
        Sessions = dir(fullfile('Z:\Data\ImagingRig3',Dates{curr_date},Animals{curr_animal}));

        Sessions = {Sessions(cellfun(@(x) ~isempty(strfind(x, 'Rec')), {Sessions.name})).name};
        for curr_session = 1:length(Sessions)
            OptoOn_times{curr_session} = [];
            OptoOff_times{curr_session} = [];
            DataPath = fullfile('Z:\Data\ImagingRig3',Dates{curr_date},Animals{curr_animal},Sessions{curr_session},[Animals{curr_animal}(1:2) '0001']);

            TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace',Dates{curr_date});
            if ~exist(TargetPath)
                mkdir(TargetPath);
            end
            data = MomEphus_load_xsg_continuous(DataPath);
            xsg_sample_rate = 10000;
            for curr_channel = 1:length(data.channel_names)
                channel = data.channel_names{curr_channel};
                channel_index = cellfun(@(x) strcmp(x,channel),data.channel_names);
                temp_data = data.channels(:,channel_index);
                switch channel
                    case 'Frame' % 2p imaging
                        % Get frame times (in seconds) from frame trigger trace                        
                        Frame_times{curr_session} = (find(temp_data(2:end) > 2.5 & temp_data(1:end-1) < 2.5) + 1)/xsg_sample_rate;
                        clear temp_data
                    case 'Lever'
                        [Lever_active{curr_session}, Lever_force_resample{curr_session}, lever_force_smooth, lever_velocity_envelope_smooth,]...
                            = CR_parseLeverMovement_Updated(temp_data, 500, 150, active_threshold); % 1K Hz
                        clear temp_data lever_force_smooth lever_velocity_envelope_smooth
                    case 'Lick'
                        [Licking_active{curr_session}, Licking_bout{curr_session}, Licking_bout_switch{curr_session}] = CR_GetLickingState(temp_data,3);
                        clear temp_data
                    case 'Opto'
                        OptoOn_times{curr_session} = (find(temp_data(2:end) > 2.5 & temp_data(1:end-1) < 2.5) + 1)/xsg_sample_rate;
                        OptoOff_times{curr_session} = (find(temp_data(2:end) < 2.5 & temp_data(1:end-1) > 2.5) + 1)/xsg_sample_rate;
                        clear temp_data
                    case 'Trial_number'
                        [TrialInfo_xsg{curr_session}] = CR_ReadBitCodeShift(temp_data, 0);
                        clear temp_data
                end
            end
        clear data
        if sum(diff(TrialInfo_xsg{curr_session}(:,2))~=1)>0
%             break
            TrialInfo_xsg{curr_session}(:,2) = 1:length(TrialInfo_xsg{curr_session});
            disp('Trial index Warning')
        end
        end
        saveFile = [TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'];
        save(saveFile,'Frame_times','Lever_active','Lever_force_resample','Licking_active','Licking_bout','Licking_bout_switch','OptoOn_times','OptoOff_times','TrialInfo_xsg','-v7.3');
        clear Frame_times Lever_active Lever_force_resample Licking_active Licking_bout Licking_bout_switch OptoOn_times OptoOff_times TrialInfo_xsg
    end
end
clear all

%% Double check trial index for bscope 1
Animals = {'CR_4412582-R','CR_4412583-LR'};
Dates = {'210911','210912','210913','210914','210915','210916','210917','210918','210919','210920','210921','210922','210923','210924','210925','210926','210927','210928','210929','210930','211001','211002'};
for curr_animal = 2
    if ismember(Animals{curr_animal},{'CR_3619073-L','CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R','CR_4383143-L','CR_4383143-R','CR_4429262-O','CR_4383144-R','CR_4412582-R','CR_4412583-LR'})
        IN = 'ChAT';
    end
    for curr_date = 1:length(Dates)
        Diff_Bitcode_Dispatcher = [];
        Diff_Bitcode_Xsg = [];
        Bitcode_Dispatcher = [];
        TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace',Dates{curr_date});
    	load([TargetPath filesep Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat'],'TrialInfo_xsg');
        BpodPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'Bpod');
        load([BpodPath filesep 'CR_' Dates{curr_date} Animals{curr_animal}(3:end) '.mat'],'saved_history');
        TrialInfo_Dispatcher = saved_history.ProtocolsSection_parsed_events;
        TrialNum_Total = length(TrialInfo_Dispatcher);
        for kk = 1:TrialNum_Total
            Bitcode_Dispatcher(kk,1) = TrialInfo_Dispatcher{kk,1}.states.bitcode(1);
        end
        Diff_Bitcode_Dispatcher(:,1) = 2:length(Bitcode_Dispatcher);
        Diff_Bitcode_Dispatcher(:,2) = diff(Bitcode_Dispatcher);
        Diff_Bitcode_Xsg(:,1) = TrialInfo_xsg{1}(2:end,2);
        Diff_Bitcode_Xsg(:,2)= diff(TrialInfo_xsg{1}(:,1));
        figure(curr_date); hold on;
        plot(Diff_Bitcode_Dispatcher(:,1),Diff_Bitcode_Dispatcher(:,2),'k*');
        plot(Diff_Bitcode_Xsg(:,1),Diff_Bitcode_Xsg(:,2),'ro');
    end
end

clear all

%% Check bad trial index
clear all
close all
clc

Animals = {'CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'};
for curr_animal = 1:length(Animals)
    if ismember(Animals{curr_animal},{'KP_3463808_1','KP_3463808_2','KP_3475729_LR','KP_3480351_1','CR_3526643-R','WL_3526641-R','WL_3526642-L','WL_3526642-R','CR_3526578-L','CR_3702608-O','CR_3702608-LR','CR_3658844-L','CR_3658844-R','CR_3658845-R','CR_3672031-O','CR_4303048-L','CR_4303048-R','CR_4303048-LR','CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{'KP_3459921_1','KP_3461990_1','WL_3526549-R','WL_3526580-L','WL_3526580-O','WL_3526578-O','WL_3526578-R','WL_3547272-O','WL_3547272-L','WL_3547273-R','WL_3547273-LR','CR_3672035-L','CR_3672035-R','CR_3702224-O','CR_3702224-L','CR_3702224-R','CR_3702226-O','CR_3786142-O','CR_3786142-L','CR_3786142-R','CR_3887041-L','CR_3887041-R','CR_3936483-O'})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3619073-L','CR_3619073-R','CR_3619073-LR','CR_3619074-O'})
        IN = 'ChAT';
    elseif ismember(Animals{curr_animal},{'CR_3575263-L','CR_3575263-R','CR_3575265-O','CR_3575265-LR'})
        IN = 'DBH';
    elseif ismember(Animals{curr_animal},{'CR_3672031-L','CR_3786160-LL','CR_3786159-O'})
        IN = 'VIP_ChAT_Cas';
    elseif ismember(Animals{curr_animal},{'CR_3974574-O','CR_3974574-L','CR_3974574-R','CR_3974573-O','CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-L','CR_4017421-R','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042830-R','CR_4042832-O','CR_4042832-L','CR_4042832-R','CR_4113794-O','CR_4113794-L','CR_4113794-R','CR_4153298-O','CR_4153298-L','CR_4153298-R','CR_4153298-LR',...
            'CR_4153299-O','CR_4153299-L','CR_4153369-O','CR_4153382-O','CR_4153382-L','CR_4153299-R','CR_4153369-L','CR_4153382-R','CR_4153382-LR','CR_4259761-O','CR_4259761-L','CR_4259761-R','CR_4259761-LR','CR_4259762-O','CR_4259762-L','CR_4259762-LL'})
        IN = 'VIP_SOM';
    elseif ismember(Animals{curr_animal},{'CR_4113793-O','CR_4113793-L','CR_4153383-O','CR_4153383-L','CR_4153383-R','CR_4153383-LR'})
        IN = 'VIP_hM4Di';
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN',IN,Animals{curr_animal},'LeverTrace');
    Dates = dir(TargetPath);
    Dates = {Dates.name};
    Dates = Dates(3:end);
    for curr_date = 1:length(Dates)
        if ~exist(fullfile(TargetPath,Dates{curr_date},[Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat']))
            continue
        end
        load(fullfile(TargetPath,Dates{curr_date},[Animals{curr_animal} '_' Dates{curr_date} '_EphusTraces.mat']),'TrialInfo_xsg');
        for curr_session = 1:length(TrialInfo_xsg)
            tempIndex = diff(TrialInfo_xsg{curr_session});
            if any(tempIndex(:,2)~=1)
                disp([Animals{curr_animal} '_' Dates{curr_date} '_' num2str(curr_session)]);
            end
        end
    end
    disp('Done');
    clear tempIndex
end