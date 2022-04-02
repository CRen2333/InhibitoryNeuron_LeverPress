%% ***** Formal Coding V1.0 *****
clear all
close all
clc

INs = {'VIP'};
Animals = {'CR_3438544-O','CR_3438544-L','CR_3438544-R'};


for curr_IN = 1:length(INs)
    IN = INs{curr_IN};
    disp(IN);
    cd(['Z:\People\Chi\WFLP_IN' filesep IN]);
     
    for curr_animal = 1:length(Animals)
        
        %% Clean
        clearvars -except Animals INs IN Initial curr_animal;
        clc;
        Animal = Animals{curr_animal};

        %% Load data path
        if ispc
            Data_path = ['Z:\People\Chi\WFLP_IN' filesep IN filesep Animal];
        elseif ismac
            Data_path = ['/Volumes/lab/People/Chi/WFLP_IN' filesep IN filesep Animal];
        end
        
        %% Find out data path and dat
        % For xsg
        Xsg_path = [Data_path filesep 'LeverTrace'];
        Xsg_folder_list = dir(Xsg_path);
        Xsg_folder_list = Xsg_folder_list(3:end);
        Date_list = cellfun(@(x) x(1:6), {Xsg_folder_list.name}, 'UniformOutput' ,false);
        
        %% Load Dispatcher and xsg files for current day
        disp([IN ' ' Animal]);
        
        ii = 0;
        Not_consRecord = {};
        Not_cons = {};
        % Basically go throuh trial by trial
        for curr_day = 1:length(Date_list)
            
            Date = Date_list{curr_day};
            disp(Date);

            %% Load xsg file list
            Xsg_curr_path = [Data_path filesep 'LeverTrace' filesep Date];
            Xsg_file_list = dir(Xsg_curr_path);
            % Should have multiple filenames, use cell, otherwise will only save the
            % first one.
            Xsg_filename_list = {Xsg_file_list(cellfun(@(x) ~isempty(strfind(x, 'CR0001AAAA')), {Xsg_file_list.name})).name};
            if isempty(Xsg_filename_list)
                warning(['No xsg file saved for ' Animal ' on ' Date '.']);
                continue
            end
            Xsg_filename_list = sort(Xsg_filename_list); % from 1st to end

            % Initialization for each day
            TrialInfo_xsg = cell(size(Xsg_filename_list));

            %% Load xsg file sequentially

            % User interaction
            disp(['Running ' Animal ' on ' Date]);

            for curr_xsg = 1:length(Xsg_filename_list)
                if strcmp(Animal,'3233232-R') && strcmp(Date,'170723') && curr_xsg == 1
                    continue
                end
                % Loading
                load([Xsg_curr_path filesep Xsg_filename_list{curr_xsg}],'-MAT'); % Because cd is not in the LeverTrace folder.
                % Get xsg sample rate
                Sample_rate_xsg = header.acquirer.acquirer.sampleRate; % Should be 10000
                % Get recording of bitcode, lever trace and licking
                Bitcode_xsg = data.acquirer.trace_1;
                % Get trial index(2nd column) and bitcode time(1st column) from xsg
                if strcmp(IN,'ChAT_eOPN3') && str2num(Date)<210623
                    if length(strfind(Xsg_filename_list{curr_xsg},'0001')) == 2
                        TrialInfo_xsg{curr_xsg} = [];
                    elseif length(strfind(Xsg_filename_list{curr_xsg},'0002')) == 1 % first xsg CR0001AAAA0001.xsg
                        [TrialInfo_xsg{curr_xsg}, shiftTimePoint] = CR_ReadBitCodeShift_firstxsg(Bitcode_xsg);
                    else
                        [TrialInfo_xsg{curr_xsg}] = CR_ReadBitCodeShift(Bitcode_xsg, shiftTimePoint);
                    end
                else
                    if length(strfind(Xsg_filename_list{curr_xsg},'0001')) == 2 % first xsg CR0001AAAA0001.xsg
                        [TrialInfo_xsg{curr_xsg}, shiftTimePoint] = CR_ReadBitCodeShift_firstxsg(Bitcode_xsg);
                    else
                        [TrialInfo_xsg{curr_xsg}] = CR_ReadBitCodeShift(Bitcode_xsg, shiftTimePoint);
                    end
                end
                % Ignore if no trials exist
                if isempty(TrialInfo_xsg{curr_xsg})
                    disp(['No Bitecode in ' Xsg_filename_list{curr_xsg} '.'])
                    continue
                end
                % For bad linux machine: trial index error reading
                TrialInfo_xsg{curr_xsg} = TrialInfo_xsg{curr_xsg}(TrialInfo_xsg{curr_xsg}(:,2)<300,:); % We never go beyond 300 trials
                % If trials are not consecutive, delete the offending trial
                consec_trials = diff([0; TrialInfo_xsg{curr_xsg}(:,2)]);
                if sum((consec_trials(2:end) ~= 1)) > 0
                    Not_cons{curr_day,curr_xsg} = sum((consec_trials(2:end) ~= 1));
                    temp_record = [IN '_' Animal '_' Date '_' Xsg_filename_list{curr_xsg}];
                    ii = ii+1;
                    Not_consRecord{ii,1} = temp_record;
                else
                    Not_cons{curr_day,curr_xsg} = 0;
                end
            end
            TrialInfo_xsg_all{curr_day} = TrialInfo_xsg;
        end
        TargetPath = ['Z:\People\Chi\WFLP_IN\CheckNonconsTrials'];
        cd(TargetPath);
                
        disp('Saving...');
        TargetPath = ['Z:\People\Chi\WFLP_IN\CheckNonconsTrials'];
        if ~isdir(TargetPath)
            mkdir(TargetPath);
        end
        save([TargetPath filesep IN '_' Animal '_CheckNonconsTrials.mat'], 'Date_list', 'Not_cons', 'Not_consRecord', 'TrialInfo_xsg_all','-mat');
    end
end

%% Manual inspection and correction
% load the file
clear;
IN = 'VIP';
Animal = 'CR_3438544-R';
load([IN '_' Animal '_CheckNonconsTrials.mat']);
if exist([IN '_' Animal '_CheckNonconsTrials_Corrected.mat'])
    aaa = load([IN '_' Animal '_CheckNonconsTrials_Corrected.mat']);
end
Corrected_date = length(aaa.Date_list);
TrialInfo_xsg_all(1:Corrected_date) = aaa.TrialInfo_xsg_all;
disp(['Corrected:' num2str(Corrected_date)]);
clear aaa;

%% Eyeballing 
clear all
close all
clc

cd('Z:\People\Chi\WFLP_IN\CheckNonconsTrials');
Allfiles = dir('*.mat');
Allfiles = {Allfiles.name}';
Allfiles = Allfiles(cellfun(@(x) ~isempty(strfind(x, 'Corrected')), Allfiles));
Allfiles = sort(Allfiles);
Checkfiles = Allfiles(68:69);
% Manual control which animal to check
ii = 2; %1

load(Checkfiles{ii});
disp(Checkfiles{ii});
IN = Checkfiles{ii}(1:10);
% IN = Checkfiles{ii}(1:9);
switch IN
    case 'PV'
        Animal = Allfiles{ii}(4:15);
    case 'CHAT'
        Animal = Allfiles{ii}(6:17);
    case 'ChAT_Cas'
        Animal = Allfiles{ii}(10:21);
    case 'VIP_DREADD'
        Animal = Allfiles{ii}(12:23);
    case 'AChR_Ago_Nai'
        Animal = Checkfiles{ii}(14:25);
    case 'ChAT_eOPN3'
        Animal = Checkfiles{ii}(12:23);
    case 'VIP_ArchT'
        Animal = Checkfiles{ii}(11:22);
    otherwise
        Animal = Allfiles{ii}(5:16);
end

% Select date
for jj = [1:30]
    clear Bitcode_Dispatcher Diff_Bitcode_Dispatcher TrialInfo_Dispatcher
%     Data_path = ['Z:\People\Chi\WFLP_IN' filesep IN filesep Animal filesep 'ForBhv' filesep 'Dispatcher'];
    Data_path = ['Z:\People\Chi\WFLP_IN' filesep IN filesep Animal filesep 'Dispatcher'];
    warning off
    if exist([Data_path filesep 'CR_' Date_list{jj} '_' Animal(4:end) '.mat'],'file') == 2
        load([Data_path filesep 'CR_' Date_list{jj} '_' Animal(4:end) '.mat'])
    else
        disp([Date_list{jj} ': No dispatcher'])
        continue
    end
    warning on
    TrialInfo_Dispatcher = saved_history.ProtocolsSection_parsed_events;
    TrialNum_Total = length(TrialInfo_Dispatcher);
    for kk = 1:TrialNum_Total
        Bitcode_Dispatcher(kk,1) = TrialInfo_Dispatcher{kk,1}.states.bitcode(1);
    end
    Diff_Bitcode_Dispatcher(:,1) = 2:TrialNum_Total;
    Diff_Bitcode_Dispatcher(:,2) = diff(Bitcode_Dispatcher);
    figure(jj)
    hold on
    plot(Diff_Bitcode_Dispatcher(:,1),Diff_Bitcode_Dispatcher(:,2),'k*');

    for curr_xsg = 1:length(TrialInfo_xsg_all{jj})
        if isempty(TrialInfo_xsg_all{jj}{curr_xsg})
            continue
        end
        TrialInfo_xsg_all{jj}{curr_xsg}(2:end,3) = diff(TrialInfo_xsg_all{jj}{curr_xsg}(:,1));
        TrialInfo_xsg_all{jj}{curr_xsg}(1,3) = nan;
        plot(TrialInfo_xsg_all{jj}{curr_xsg}(2:end,2),diff(TrialInfo_xsg_all{jj}{curr_xsg}(:,1)),'ro');
    end
end

close all

save(Checkfiles{ii},'TrialInfo_xsg_all','-append');
clear all
        
                
    

