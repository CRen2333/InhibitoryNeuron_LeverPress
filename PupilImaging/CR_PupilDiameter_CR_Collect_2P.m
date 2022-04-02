%% Collect Pupil Dia
clear all
close all
clc

Animals = {'CR_3633170-L'};

for curr_animal = 1:length(Animals)
    clearvars -except Animals curr_animal
    Animal = Animals{curr_animal};
    if ismember(Animals{curr_animal},{'CR_3526643-R','WL_3526641-R','WL_3526642-L','WL_3526642-R'})
        IN = 'VIP';
    elseif ismember(Animals{curr_animal},{'WL_3526549-R','WL_3526580-L','WL_3526580-O','WL_3526578-O'})
        IN = 'SOM';
    elseif ismember(Animals{curr_animal},{'CR_3575265-O'})
        IN = 'DBH';
    end
    SourcePath = fullfile('Z:\People\Chi\TwoP_IN\PupilFitting',Animal);
    Date_list = dir(SourcePath);
    Date_list = {Date_list.name};
    Date_list = Date_list(3:end);
    Date_list = sort(Date_list)';
    
    task_start_date = Date_list{1};
    index = find(ismember(Date_list,task_start_date));
    Date_list = Date_list(index:end);
    
    for curr_day = 1:length(Date_list)
        Date = Date_list{curr_day};
        tic
        Sessions = dir([SourcePath filesep Date]);
        Sessions = {Sessions.name};
        Sessions = Sessions(cellfun(@(x) contains(x, 'Rec'), Sessions));
        Sessions = sort(Sessions);
        for curr_session = 1:length(Sessions)
            Session = Sessions{curr_session};
            disp([Animal ' ' Date ' ' Session]);
            file_list = dir([SourcePath filesep Date filesep Session]);
            file_list = {file_list.name};
            Baseline_file = file_list(cellfun(@(x) contains(x, '_Baseline_PupilDia'), file_list));
            file_list = file_list(cellfun(@(x) contains(x, '_Seg')&&contains(x, '_Fitting'), file_list));
            file_list = sort(file_list);
            if isempty(file_list)
                PupilDia_refine_all{curr_day,1}{curr_session,1} = [];
                PupilDia_refine_all_medfilter{curr_day,1}{curr_session,1} = [];
                PupilDia_refine_all_norm{curr_day,1}{curr_session,1} = [];
                PupilDia_refine_all_medfilter_norm{curr_day,1}{curr_session,1} = [];
                PupilDia_refine_all_norm_butterworth{curr_day,1}{curr_session,1} = [];
                Derivative_PupilDia_butterworth{curr_day,1}{curr_session,1} = [];
                dil_start_end{curr_day,1}{curr_session,1} = [];
                cons_start_end{curr_day,1}{curr_session,1} = [];
                sigphase{curr_day,1}{curr_session,1} = [];
                PupilPos_refine_all{curr_day,1}{curr_session,1} = [];
                PupilPos_refine_all_medfilter{curr_day,1}{curr_session,1} = [];
                Laser{curr_day,1}{curr_session,1} = [];
                continue
            end
        
            % Baseline
            load([SourcePath filesep Date filesep Session filesep Baseline_file{1}],'Baseline_Pupil_Dia');
            for curr_file = 1:length(file_list)
                curr_file_name = file_list{curr_file};
                load([SourcePath filesep Date filesep Session filesep curr_file_name],'PupilDia_refine','bestFits_refine','Laser_syn','-mat');                 
                PupilDia_refine_all{curr_day,1}{curr_session,1}{curr_file,1} = PupilDia_refine;
                % [x0 y0 a b alpha score]
                for ii = 1:length(bestFits_refine)
                    if isempty(bestFits_refine{ii})
                        PupilPos_refine_all{curr_day,1}{curr_session,1}{curr_file,1}(ii,1:2) = nan;
                    else
                        PupilPos_refine_all{curr_day,1}{curr_session,1}{curr_file,1}(ii,:) = bestFits_refine{ii}(1,1:2);
                    end
                end
                if size(Laser_syn,2) > 1 % coding by mistake
                    Laser_syn = Laser_syn(1,:)';
                end
                Laser{curr_day,1}{curr_session,1}{curr_file,1} = Laser_syn;
                clear PupilDia_refine bestFits_refine
            end
        	
            PupilDia_refine_all{curr_day,1}{curr_session,1} = cell2mat(PupilDia_refine_all{curr_day,1}{curr_session,1});
            PupilDia_refine_all_medfilter{curr_day,1}{curr_session,1} = medfilt1(PupilDia_refine_all{curr_day,1}{curr_session,1},5);
            PupilDia_refine_all_norm{curr_day,1}{curr_session,1} = PupilDia_refine_all{curr_day,1}{curr_session,1}./Baseline_Pupil_Dia;
            PupilDia_refine_all_medfilter_norm{curr_day,1}{curr_session,1} = PupilDia_refine_all_medfilter{curr_day,1}{curr_session,1}./Baseline_Pupil_Dia;
            [PupilDia_refine_all_norm_butterworth{curr_day,1}{curr_session,1},Derivative_PupilDia_butterworth{curr_day,1}{curr_session,1},dil_start_end{curr_day,1}{curr_session,1},cons_start_end{curr_day,1}{curr_session,1},sigphase{curr_day,1}{curr_session,1}] = CR_GetPupil(PupilDia_refine_all_norm{curr_day,1}{curr_session,1},1,15,10);
            PupilPos_refine_all{curr_day,1}{curr_session,1} = cell2mat(PupilPos_refine_all{curr_day,1}{curr_session,1});
            PupilPos_refine_all_medfilter{curr_day,1}{curr_session,1}(:,1) = medfilt1(PupilPos_refine_all{curr_day,1}{curr_session,1}(:,1),5);
            PupilPos_refine_all_medfilter{curr_day,1}{curr_session,1}(:,2) = medfilt1(PupilPos_refine_all{curr_day,1}{curr_session,1}(:,2),5);
            Laser{curr_day,1}{curr_session,1} = cell2mat(Laser{curr_day,1}{curr_session,1});
        toc
        end
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN', IN, Animal, 'MovAnalysis');
    if exist([TargetPath filesep Animal '_PupilTraces.mat'],'file')
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilDia_refine_all','PupilDia_refine_all_norm','PupilDia_refine_all_medfilter','PupilDia_refine_all_medfilter_norm','-append');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilDia_refine_all_norm_butterworth','Derivative_PupilDia_butterworth','dil_start_end','cons_start_end','sigphase','-append');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilPos_refine_all','PupilPos_refine_all_medfilter','-append');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'Laser','-append');
    else
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilDia_refine_all','PupilDia_refine_all_norm','PupilDia_refine_all_medfilter','PupilDia_refine_all_medfilter_norm','-v7.3');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilDia_refine_all_norm_butterworth','Derivative_PupilDia_butterworth','dil_start_end','cons_start_end','sigphase','-append');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'PupilPos_refine_all','PupilPos_refine_all_medfilter','-append');
        save([TargetPath filesep Animal '_PupilTraces.mat'],'Laser','-append');
    end
end

