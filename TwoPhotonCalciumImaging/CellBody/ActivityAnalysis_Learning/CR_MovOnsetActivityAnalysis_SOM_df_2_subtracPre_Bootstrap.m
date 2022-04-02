%% Neurons activated or suppressed during movement
clear all;
close all;
clc;

IN = 'SOM';
Animals = {'KP_3459921_1','KP_3461990_1','WL_3526578-O','WL_3526580-O','WL_3547273-LR','WL_3547273-R','CR_3786142-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'};

% Classify neuron base on bootstrap
nboot = 1000;
alpha_threshold = [0.05];
frame_threshold = 0.25;

for curr_animal = 1:2%length(Animals)
    clear df_f_MovOnset_Alinged_sub
    Animal = Animals{curr_animal};    
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    DataPath_boot = [DataPath filesep 'Boot'];
    if ~exist(DataPath_boot)
        mkdir(DataPath_boot);
    end
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','-mat');
    load([DataPath filesep Animal '_ImagingInfo.mat'],'Imaging_Fields','-mat');
    if strcmp(Animal,'KP_3459921_1')
        for curr_roi = 1:size(df_f_MovOnset_Alinged{3,8},1)
            df_f_MovOnset_Alinged{3,8}{curr_roi,1} = df_f_MovOnset_Alinged{3,8}{curr_roi,1}([1:2:72],:);
        end
        for curr_roi = 1:size(df_f_MovOnset_Alinged{4,8},1)
            df_f_MovOnset_Alinged{4,8}{curr_roi,1} = df_f_MovOnset_Alinged{4,8}{curr_roi,1}([2:2:72],:);
        end
    end
    
    for curr_field = 1:size(df_f_MovOnset_Alinged,1)
        for curr_session = 1:min(11,size(df_f_MovOnset_Alinged(curr_field,:),2))
            if isempty(df_f_MovOnset_Alinged{curr_field,curr_session})
                 df_f_MovOnset_Alinged_sub{curr_field,curr_session} = [];
                 continue
            end
            if size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},1) == 36
                MovOnset_Frame = 8;
                Baseline_Frame = [3:5];
            elseif size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},1) == 72
                MovOnset_Frame = 15;
                Baseline_Frame = [7:10];
            elseif size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},1) == 76
                MovOnset_Frame = 16;
                Baseline_Frame = [8:11];
            end
            
            for curr_roi = 1:size(df_f_MovOnset_Alinged{curr_field,curr_session},1)
                temp = df_f_MovOnset_Alinged{curr_field,curr_session}{curr_roi};
                df_f_MovOnset_Alinged_sub{curr_field,curr_session}{curr_roi,1} = temp - repmat(nanmean(temp(Baseline_Frame,:),1),size(temp,1),1);                                
            end
        end
    end
    clear df_f_MovOnset_Alinged
    temp = df_f_MovOnset_Alinged_sub;
    clear df_f_MovOnset_Alinged_sub
    for curr_field = 1:size(temp,1)
        for curr_roi = 1:size(temp{curr_field,1},1)
            df_f_MovOnset_Alinged_sub{curr_field,1}{curr_roi,1} = temp{curr_field,1}{curr_roi,1};
            df_f_MovOnset_Alinged_sub{curr_field,2}{curr_roi,1} = [];
            for ii = 2:4
                if isempty(temp{curr_field,ii})
                    continue
                else
                    df_f_MovOnset_Alinged_sub{curr_field,2}{curr_roi,1} = [df_f_MovOnset_Alinged_sub{curr_field,2}{curr_roi,1},temp{curr_field,ii}{curr_roi,1}];
                end
            end
            df_f_MovOnset_Alinged_sub{curr_field,3}{curr_roi,1} = [];
            for ii = 5:8
                if isempty(temp{curr_field,ii})
                    continue
                else
                    df_f_MovOnset_Alinged_sub{curr_field,3}{curr_roi,1} = [df_f_MovOnset_Alinged_sub{curr_field,3}{curr_roi,1},temp{curr_field,ii}{curr_roi,1}];
                end
            end
            df_f_MovOnset_Alinged_sub{curr_field,4}{curr_roi,1} = [];
            for ii = 9:size(temp,2)
                if isempty(temp{curr_field,ii})
                    continue
                else
                    df_f_MovOnset_Alinged_sub{curr_field,4}{curr_roi,1} = [df_f_MovOnset_Alinged_sub{curr_field,4}{curr_roi,1},temp{curr_field,ii}{curr_roi,1}];
                end
            end
        end
    end
    clear temp
    for curr_session = 1:size(df_f_MovOnset_Alinged_sub,2)
        for curr_field = 1:size(df_f_MovOnset_Alinged_sub,1)
            save_bootstat_compare_matrix = {};
            if isempty(df_f_MovOnset_Alinged_sub{curr_field,curr_session}) || size(df_f_MovOnset_Alinged_sub{curr_field,curr_session}{1,1},2)<10
                Label.MvmModuLabel_Act{curr_field,1}(:,curr_session) = nan;
                Label.MvmModuLabel_Sup{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_ActOnset{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_ActTHalf{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_ActTPeak{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_SupOnset{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_SupTHalf{curr_field,1}(:,curr_session) = nan;
                Time.MvmModuLabel_SupTPeak{curr_field,1}(:,curr_session) = nan;
                save([DataPath_boot filesep Animal '_Field' num2str(curr_field) '_Session' num2str(curr_session) '_BootMatrix.mat'],'save_bootstat_compare_matrix','-v7.3');
                continue
            end
            frame_rate = Imaging_Fields{curr_field}.FrameRate(curr_session);
            for curr_roi = 1:length(df_f_MovOnset_Alinged_sub{curr_field,curr_session})
                disp([Animal ' Session' num2str(curr_session) ' Field' num2str(curr_field) ' ROI' num2str(curr_roi) '/' num2str(length(df_f_MovOnset_Alinged_sub{curr_field,curr_session}))]);
                temp_matrix = df_f_MovOnset_Alinged_sub{curr_field,curr_session}{curr_roi,1};
                if sum(isnan(sum(temp_matrix))) == size(temp_matrix,2)
                    save_bootstat_compare_matrix{curr_roi,1} = [];
                    Label.MvmModuLabel_Act{curr_field,1}(curr_roi,curr_session) = nan;
                    Label.MvmModuLabel_Sup{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_ActOnset{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_ActTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_ActTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_SupOnset{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_SupTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_SupTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                continue
                end
                if size(temp_matrix,1) == 36
                    MovOnset_Frame = 8;
                    Baseline_Frame = [3:5];
                elseif size(temp_matrix,1) == 72
                    MovOnset_Frame = 15;
                    Baseline_Frame = [7:10];
                elseif size(temp_matrix,1) == 76
                    MovOnset_Frame = 16;
                    Baseline_Frame = [8:11];
                end
                baseline_matrix = temp_matrix(Baseline_Frame,:);
                baseline_matrix = baseline_matrix(:);
                baseline_matrix = baseline_matrix(~isnan(baseline_matrix));
                [bootstat_baseline,~] = bootstrp(nboot,@nanmean,baseline_matrix);
                % No upsampling
                for curr_frame = 1:size(temp_matrix,1)
                    [bootstat_frame(:,curr_frame),~] = bootstrp(nboot,@nanmean,temp_matrix(curr_frame,:));
                end
                compare_matrix = bootstat_frame-repmat(bootstat_baseline,1,size(temp_matrix,1));
                lower_bound = prctile(compare_matrix,alpha_threshold/2*100,1);
                upper_bound = prctile(compare_matrix,100-alpha_threshold/2*100,1);
                binary_label = abs(lower_bound>0)+abs(upper_bound<0);
                binary_label(1:Baseline_Frame(end)) = 0;
                lower_bound = lower_bound(Baseline_Frame(end)+1:end);
                upper_bound = upper_bound(Baseline_Frame(end)+1:end);
                if length(binary_label) == 36
                    binary_label_2 =  binary_label+[binary_label(2:end),0]+[binary_label(3:end),0,0];
                    temp_onset = find(binary_label_2==3,1,'first');
                else
                    binary_label_2 =  binary_label+[binary_label(2:end),0]+[binary_label(3:end),0,0]+[binary_label(4:end),0,0,0]+[binary_label(5:end),0,0,0,0];
                    temp_onset = find(binary_label_2==5,1,'first');
                end
                % Activated
                if sum(lower_bound>0)>=length(lower_bound)*frame_threshold
                    Label.MvmModuLabel_Act{curr_field,1}(curr_roi,curr_session) = 1;
                    if isempty(temp_onset)
                        Time.MvmModuLabel_ActOnset{curr_field,1}(curr_roi,curr_session) = nan;
                        Time.MvmModuLabel_ActTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                        Time.MvmModuLabel_ActTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                    else
                        Time.MvmModuLabel_ActOnset{curr_field,1}(curr_roi,curr_session) = (temp_onset-MovOnset_Frame)/frame_rate;
                        temp_mean = nanmean(temp_matrix,2);
                        temp_mean(1:temp_onset-1) = nan;
                        [temp_max,temp_max_I] = nanmax(temp_mean);
                        temp_trace_diff = abs(temp_mean(1:temp_max_I)-(temp_max/2+temp_mean(temp_onset)/2));
                        [~,temp_I_half] = sort(temp_trace_diff,'ascend');
                        temp_I_half(temp_I_half<temp_onset) = [];
                        if length(temp_I_half)<2
                            Time.MvmModuLabel_ActTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                        else
                            Time.MvmModuLabel_ActTHalf{curr_field,1}(curr_roi,curr_session) = (min(temp_I_half(1:2))-MovOnset_Frame)/frame_rate;
                        end
                        Time.MvmModuLabel_ActTPeak{curr_field,1}(curr_roi,curr_session) = (temp_max_I-MovOnset_Frame)/frame_rate;
                    end
                else
                    Label.MvmModuLabel_Act{curr_field,1}(curr_roi,curr_session) = 0;
                    Time.MvmModuLabel_ActOnset{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_ActTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_ActTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                end
                if sum(upper_bound<0)>=length(upper_bound)*frame_threshold
                    Label.MvmModuLabel_Sup{curr_field,1}(curr_roi,curr_session) = -1;
                    if isempty(temp_onset)
                        Time.MvmModuLabel_SupOnset{curr_field,1}(curr_roi,curr_session) = nan;
                        Time.MvmModuLabel_SupTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                        Time.MvmModuLabel_SupTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                    else
                        Time.MvmModuLabel_SupOnset{curr_field,1}(curr_roi,curr_session) = (temp_onset-MovOnset_Frame)/frame_rate;
                        temp_mean = nanmean(temp_matrix,2);
                        temp_mean(1:temp_onset-1) = nan;
                        [temp_min,temp_min_I] = nanmax(temp_mean);
                        temp_trace_diff = abs(temp_mean(1:temp_min_I)-(temp_min/2+temp_mean(temp_onset)/2));
                        [~,temp_I_half] = sort(temp_trace_diff,'ascend');
                        temp_I_half(temp_I_half<temp_onset) = [];
                        if length(temp_I_half)<2
                            Time.MvmModuLabel_SupTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                        else
                            Time.MvmModuLabel_SupTHalf{curr_field,1}(curr_roi,curr_session) = (min(temp_I_half(1:2))-MovOnset_Frame)/frame_rate;
                        end
                        Time.MvmModuLabel_SupTPeak{curr_field,1}(curr_roi,curr_session) = (temp_min_I-MovOnset_Frame)/frame_rate;
                    end
                else
                    Label.MvmModuLabel_Sup{curr_field,1}(curr_roi,curr_session) = 0;
                    Time.MvmModuLabel_SupOnset{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_SupTHalf{curr_field,1}(curr_roi,curr_session) = nan;
                    Time.MvmModuLabel_SupTPeak{curr_field,1}(curr_roi,curr_session) = nan;
                end
                save_bootstat_compare_matrix{curr_roi,1} = compare_matrix;
                clear binary_label binary_label_2 temp_matrix basline_matrix bootstat_frame compare_matrix temp_mean temp_I_half temp_max_I temp_min_I
            end
            save([DataPath_boot filesep Animal '_Field' num2str(curr_field) '_Stage' num2str(curr_session) '_BootMatrix.mat'],'save_bootstat_compare_matrix','-v7.3');
        end
    end
    save([DataPath filesep Animal '_Activity_bootstats_stage.mat'],'Label','Time','-v7.3');
    clear save_bootstat_compare_matrix Label Time Imaging_Fields
end

  