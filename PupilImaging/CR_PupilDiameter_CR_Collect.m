%% Collect Pupil Dia
clear all
close all
clc

IN = 'VIP';
Initial = 'CR';
Animals = {'3438544-O','3438544-L','3438544-R'};

for curr_animal = 1:length(Animals)
    clearvars -except IN Initial Animals curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Pupil']);
    Date_list = dir(cd);
    Date_list = {Date_list.name};
    Date_list = Date_list(3:end);
    Date_list = sort(Date_list)';
    
    if ismember(Animal,{'3271320-L','3333408-O','3333408-L','3333408-R'})
        task_start_date = '171105';        
    elseif ismember(Animal, {'3258531-O','3258531-L'})
        task_start_date = '171129';
    elseif ismember(Animal, {'3358884-L','3358884-R'})
        task_start_date = '171213';
    elseif ismember(Animal, {'3358883-O','3358883-R'})
        task_start_date = '171214';
    elseif ismember(Animal, {'3373693-O','3373693-L'})
        task_start_date = '180111';
    elseif ismember(Animal, {'3373693-R','3373693-LR'})
        task_start_date = '180117';
    elseif ismember(Animal, {'3420509-L','3420509-R'})
        task_start_date = '180206';
    elseif ismember(Animal, {'3358845-O','3358845-LR'})
        task_start_date = '180313';
    elseif ismember(Animal, {'3373785-O','3373785-L','3373785-R'})
        task_start_date = '180313';
    elseif ismember(Animal, {'3495100-R'})
        task_start_date = '180920';
    else
        task_start_date = Date_list{1};
    end
    index = find(ismember(Date_list,task_start_date));
    Date_list = Date_list(index:end);
    
    for curr_day = 1:length(Date_list)
        Date = Date_list{curr_day};
        tic
        disp(Date);
        TargetPath = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Pupil' filesep Date];
        file_list = dir([TargetPath]);
        file_list = file_list(3:end);
        file_list = {file_list.name}';
        file_list = file_list(cellfun(@(x) contains(x, '_Fitting')&&contains(x, '_Seg_'), file_list));
        file_list = sort(file_list);
        if isempty(file_list)
            PupilDia_refine_all{curr_day} = [];
            PupilDia_refine_all_medfilter{curr_day} = [];
            PupilDia_refine_all_norm{curr_day} = [];
            PupilDia_refine_all_medfilter_norm{curr_day} = [];
            PupilDia_refine_all_norm_butterworth{curr_day} = [];
            Derivative_PupilDia_butterworth{curr_day} = [];
            dil_start_end{curr_day} = [];
            cons_start_end{curr_day} = [];
            sigphase{curr_day} = [];
            continue
        end
        
        % Baseline
        load([Date filesep Initial '_' Date '_' Animal '_Baseline_PupilDia_CR.mat'],'Baseline_Pupil_Dia');
        for curr_file = 1:length(file_list)
            curr_file_name = file_list{curr_file};
            load([Date filesep curr_file_name],'PupilDia_refine','bestFits_refine');
            if ismember(Animal, {'3373693-LR','3358845-LR','3491479-LR','3438521-LR','3547207-LR','3702169-LR'})
                xsg_index = str2num(curr_file_name(21:24));
            else
                xsg_index = str2num(curr_file_name(20:23));
            end
            PupilDia_refine_all{curr_day}{xsg_index}{curr_file,1} = PupilDia_refine;
            temp = cell2mat(bestFits_refine);
            PupilPos_refine_all{curr_day}{xsg_index}{curr_file,1} = temp(:,1:2);
            clear PupilDia_refine
        end
        for ii = 1:xsg_index
            PupilDia_refine_all{curr_day}{ii} = cell2mat(PupilDia_refine_all{curr_day}{ii});
            PupilDia_refine_all_medfilter{curr_day}{ii} = medfilt1(PupilDia_refine_all{curr_day}{ii},5);
            PupilDia_refine_all_norm{curr_day}{ii} = PupilDia_refine_all{curr_day}{ii}./Baseline_Pupil_Dia;
            PupilDia_refine_all_medfilter_norm{curr_day}{ii} = PupilDia_refine_all_medfilter{curr_day}{ii}./Baseline_Pupil_Dia;
            [PupilDia_refine_all_norm_butterworth{curr_day}{ii},Derivative_PupilDia_butterworth{curr_day}{ii},dil_start_end{curr_day}{ii},cons_start_end{curr_day}{ii},sigphase{curr_day}{ii}] = CR_GetPupil(PupilDia_refine_all_norm{curr_day}{ii},1,15,10);
            PupilPos_refine_all{curr_day}{ii} = cell2mat(PupilPos_refine_all{curr_day}{ii});
        end
        toc
    end
    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'MovAnalysis']);
    if exist([Initial '_' Animal '_Traces.mat'],'file')
        save([Initial '_' Animal '_Traces.mat'],'PupilDia_refine_all','PupilDia_refine_all_norm','PupilDia_refine_all_medfilter','PupilDia_refine_all_medfilter_norm','PupilPos_refine_all','-append');
        save([Initial '_' Animal '_Traces.mat'],'PupilDia_refine_all_norm_butterworth','Derivative_PupilDia_butterworth','dil_start_end','cons_start_end','sigphase','-append');
    else
        save([Initial '_' Animal '_Traces.mat'],'PupilDia_refine_all','PupilDia_refine_all_norm','PupilDia_refine_all_medfilter','PupilDia_refine_all_medfilter_norm','PupilPos_refine_all','-v7.3');
        save([Initial '_' Animal '_Traces.mat'],'PupilDia_refine_all_norm_butterworth','Derivative_PupilDia_butterworth','dil_start_end','cons_start_end','sigphase','-append');
    end
end

%% Collect pupil diameter aligned to cued-rewarded movement onset
clear all
close all
clc

IN = 'VIP';
Initial = 'CR';
Animals = {'3438544-O','3438544-L','3438544-R'};

sample_window = [-7:90]; % movonset is 8
for curr_animal = 1:length(Animals)
    clearvars -except IN Initial Animals curr_animal sample_window
    Animal = Animals{curr_animal};
    disp(Animal);
    
    % Load index first
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Movanalysis']);
    load([Initial '_' Animal '_CuedMov.mat'],'CuedMov_SingleAnimal','-mat');
    load([Initial '_' Animal '_Traces.mat'],'PupilDia_refine_all_medfilter','PupilDia_refine_all_medfilter_norm','Derivative_PupilDia_butterworth','PupilPos_refine_all','-mat');
    
    for ii = 1:length(CuedMov_SingleAnimal)
        if isempty(CuedMov_SingleAnimal{ii})
            continue
        end
        if isempty(PupilDia_refine_all_medfilter{ii})
            continue
        end
        for jj = 1:min(length(PupilDia_refine_all_medfilter{ii}),length(CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info))            
            if isempty(CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj})
                continue
            end
            temp = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info{jj};
            Reward_index = logical(temp(:,5));
            movonset_frame = round(temp(:,2).*15);
            PupilDia_MovOnset_Aligned_all{ii}{1,jj} = [];
            PupilDia_MovOnset_Aligned_all_norm{ii}{1,jj} = [];
            PupilSpeed_MovOnset_Aligned_all_norm{ii}{1,jj} = [];
            for kk = 1:length(movonset_frame)
                curr_sample_window = movonset_frame(kk)+sample_window;
                if curr_sample_window(1) <= 0
                    temp_trace = nan(size(sample_window))';
                    temp_trace_norm = nan(size(sample_window))';
                    temp_trace_der_norm = nan(size(sample_window))';
                elseif curr_sample_window(end) > length(PupilDia_refine_all_medfilter{ii}{jj})
                    temp_trace = nan(size(sample_window))';
                    temp_trace_norm = nan(size(sample_window))';
                    temp_trace_der_norm = nan(size(sample_window))';
                else
                    temp_trace = PupilDia_refine_all_medfilter{ii}{jj}(curr_sample_window);
                    temp_trace_norm = PupilDia_refine_all_medfilter_norm{ii}{jj}(curr_sample_window);
                    temp_trace_der_norm = Derivative_PupilDia_butterworth{ii}{jj}(curr_sample_window);
                end
                PupilDia_MovOnset_Aligned_all{ii}{1,jj} = [PupilDia_MovOnset_Aligned_all{ii}{1,jj},temp_trace];
                PupilDia_MovOnset_Aligned_all_norm{ii}{1,jj} = [PupilDia_MovOnset_Aligned_all_norm{ii}{1,jj},temp_trace_norm];
                PupilSpeed_MovOnset_Aligned_all_norm{ii}{1,jj} = [PupilSpeed_MovOnset_Aligned_all_norm{ii}{1,jj},temp_trace_der_norm];
                clear temp_trace temp_trace_norm temp_trace_der_norm
            end
            PupilDia_MovOnset_Aligned_all_Rewarded{ii}{1,jj} = PupilDia_MovOnset_Aligned_all{ii}{1,jj}(:,Reward_index);
            PupilDia_MovOnset_Aligned_all_norm_Rewarded{ii}{1,jj} = PupilDia_MovOnset_Aligned_all_norm{ii}{1,jj}(:,Reward_index);
            PupilDia_MovOnset_Aligned_all_unRewarded{ii}{1,jj} = PupilDia_MovOnset_Aligned_all{ii}{1,jj}(:,~Reward_index);
            PupilDia_MovOnset_Aligned_all_norm_unRewarded{ii}{1,jj} = PupilDia_MovOnset_Aligned_all_norm{ii}{1,jj}(:,~Reward_index);
            PupilSpeed_MovOnset_Aligned_all_norm_Rewarded{ii}{1,jj} = PupilSpeed_MovOnset_Aligned_all_norm{ii}{1,jj}(:,Reward_index);
            PupilSpeed_MovOnset_Aligned_all_norm_unRewarded{ii}{1,jj} = PupilSpeed_MovOnset_Aligned_all_norm{ii}{1,jj}(:,~Reward_index);
        end
    end
    
    save([Initial '_' Animal '_CuedMov.mat'],'PupilDia_MovOnset_Aligned_all_Rewarded','PupilDia_MovOnset_Aligned_all_norm_Rewarded',...
        'PupilDia_MovOnset_Aligned_all_unRewarded','PupilDia_MovOnset_Aligned_all_norm_unRewarded','PupilDia_MovOnset_Aligned_all',...
        'PupilDia_MovOnset_Aligned_all_norm','PupilSpeed_MovOnset_Aligned_all_norm','PupilSpeed_MovOnset_Aligned_all_norm_Rewarded',...
        'PupilSpeed_MovOnset_Aligned_all_norm_unRewarded','-append');
    
    clear PupilDia_MovOnset_Aligned_all_Rewarded PupilDia_MovOnset_Aligned_all_norm_Rewarded PupilDia_MovOnset_Aligned_all_unRewarded ...
        PupilDia_MovOnset_Aligned_all_norm_unRewarded PupilDia_MovOnset_Aligned_all PupilDia_MovOnset_Aligned_all_norm...
        PupilSpeed_MovOnset_Aligned_all_norm PupilSpeed_MovOnset_Aligned_all_norm_Rewarded PupilSpeed_MovOnset_Aligned_all_norm_unRewarded
    clear CuedMov_SingleAnimal
    
end


