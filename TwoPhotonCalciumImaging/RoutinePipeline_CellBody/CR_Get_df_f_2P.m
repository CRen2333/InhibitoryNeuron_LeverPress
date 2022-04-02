%% Get df_f after mannullay checking rois
clear all;
close all;
clc;
 
IN = 'VIP';
Animals = {'WL_3526642-R'};
Rig = 'MOM';

% Get shrinked roi and background roi
for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal Rig
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['C:\Lab\Temp\ROI\' Animal];
%     General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Summed_2P_image'];
    cd([General_path filesep 'ROI']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(3:end);
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        cd([General_path filesep 'ROI' filesep Field filesep 'Fixed']);
        fixed_ROI_list = dir(cd);
        fixed_ROI_list = {fixed_ROI_list.name};
        fixed_ROI_list = sort(fixed_ROI_list);
        fixed_ROI_list = fixed_ROI_list(3:end);
        % check whether shrinked roi has been made or not
        Target_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Summed_2P_image' filesep 'ROI' filesep Field filesep 'Shrinked'];
        if ~exist(Target_path)
            mkdir(Target_path);
        end
        shrinked_ROI_list = dir(Target_path);
        shrinked_ROI_list = {shrinked_ROI_list.name};
        shrinked_ROI_list = sort(shrinked_ROI_list);
        if length(shrinked_ROI_list) > 2
            shrinked_ROI_list = shrinked_ROI_list(3:end);
        else
            shrinked_ROI_list = {};
        end
        for curr_roi_file = 1:length(fixed_ROI_list)
            if curr_roi_file <= length(shrinked_ROI_list)
                if strcmp(shrinked_ROI_list{curr_roi_file}(1:end-19),fixed_ROI_list{curr_roi_file}(1:end-10))
                    disp(['ROI file #' num2str(curr_roi_file) ' shrinked']);
                    continue
                end
            end
            roi_file = [General_path filesep 'ROI' filesep Field filesep 'Fixed' filesep fixed_ROI_list{curr_roi_file}];
            sum_file = [General_path filesep Field filesep fixed_ROI_list{curr_roi_file}(1:end-19) '.tif'];
            bg_roi_inner = 2;
            bg_roi_outer = 8;
            switch Rig
                case 'MOM'
                    roi_shrink_r = 0.8;
                case 'BScope1'
                    roi_shrink_r = 0.9; % MOM 0.8, BS1 0.9
            end
            [polygon, im, N, M] = make_bg_ring(roi_file, sum_file, bg_roi_inner, bg_roi_outer, roi_shrink_r);
            Target_filename = [fixed_ROI_list{curr_roi_file}(1:end-10) '_shrinked.roi'];
            save([Target_path filesep Target_filename],'polygon','N','M');
            clear polygon im N M
        end
    end
end

%% Get roi traces & df/f
clear all;
close all;
clc;

IN = 'VIP';
Animals = {'WL_3526642-R'};
Rig = 'MOM';

% Get shrinked roi and background roi
for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal Rig
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'Summed_2P_image' filesep 'ROI']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(3:end);
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        roi_file_folder = [General_path filesep 'Summed_2P_image' filesep 'ROI' filesep Field filesep 'Shrinked'];
        shrinked_ROI_list = dir(roi_file_folder);
        shrinked_ROI_list = {shrinked_ROI_list.name};
        shrinked_ROI_list = sort(shrinked_ROI_list);
        shrinked_ROI_list = shrinked_ROI_list(3:end);
        % Get imaging dates
        Dates = cellfun(@(x) x(4:9), shrinked_ROI_list, 'UniformOutput', false);
        Dates = unique(Dates);
        for curr_day = 1:length(Dates)
            Date = Dates{curr_day};
            disp(Date);            
            Target_path = [General_path filesep 'df_f' filesep Field];
            if ~exist(Target_path)
                mkdir(Target_path);
            end
            roi_files_curr_day = shrinked_ROI_list(cellfun(@(x) contains(x, Date), shrinked_ROI_list));
            image_folder_path = ['F:\Data\MotionCorrection\' Animal filesep Date];
            tiff_folders = dir(image_folder_path);
            tiff_folders = {tiff_folders.name};
            tiff_folders = sort(tiff_folders);
            tiff_folders = tiff_folders(cellfun(@(x) contains(x, 'Rec'), tiff_folders));
            for curr_session = 1:length(roi_files_curr_day)
                roi_filename = [roi_file_folder filesep roi_files_curr_day{curr_session}];
                switch Rig
                    case 'BScope1'
                        zPlane = 'Z1';
                    case 'MOM'
                        if ismember(Animal,{'CR_3526643-R'})
                            zPlane = 'red_Z1';
                        else
                            if ismember(Field,{'Field_1','Field_3'})
                                zPlane = 'Z1';
                            elseif ismember(Field,{'Field_2','Field_4'})
                                zPlane = 'Z2';
                            end
                        end
                end
                
                tiff_path = [image_folder_path filesep tiff_folders{curr_session} filesep zPlane filesep 'motioncorrected_tiff'];
                interp_sub = 0;
                local_comp = 0;
                switch Rig
                    case 'MOM'
                        framerate = 28; 
                        [roi_trace{curr_session},roi_trace_bg{curr_session}] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp);  
                    case 'BScope1'
                        framerate = 30;
                        [roi_trace{curr_session},roi_trace_bg{curr_session}] = CR_AP_Rig3_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp);                        
                end
                roi_trace_long{curr_session} = roi_trace{curr_session} - roi_trace_bg{curr_session};
                
                [roi_trace_df{curr_session}, roi_trace_baseline{curr_session}] = AP_baselineEstimation(roi_trace_long{curr_session},framerate);
                % select 10 rois for qualirt check
                if size(roi_trace_df{curr_session},1) < 5
                    selected_roi_index = [1:size(roi_trace_df{curr_session},1)];
                elseif size(roi_trace_df{curr_session},1) < 10                    
                    selected_roi_index = [1:5];
                elseif size(roi_trace_df{curr_session},1) < 20                    
                    selected_roi_index = [1:2:10];
                else
                    selected_roi_index = [1:2:20];
                end
                roi_for_qualitycheck{curr_day}{curr_session} = roi_trace{curr_session}(selected_roi_index,:);
            end
            Target_filename = [Animal '_' Date '_ROI_Traces.mat'];
            save([Target_path filesep Target_filename],'roi_trace','roi_trace_bg','roi_trace_df','roi_trace_baseline','-v7.3');
            clear roi_trace roi_trace_bg roi_trace_df roi_trace_baseline
        end
        save([General_path filesep 'df_f' filesep Field filesep Animal '_' Field '_ROI_Traces_Check.mat'],'roi_for_qualitycheck','-v7.3');
        clear roi_for_qualitycheck
    end
end

%% For new animals with piezo, the frame tag is not renewed --> NAN in alternative frames
clear all;
close all;
clc;
% 
IN = 'VIP';
Animals = {'WL_3526642-R'};
Rig = 'MOM';

for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'df_f']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(cellfun(@(x) ~isempty(strfind(x, 'Field')), Fields));
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        df_f_folder = [General_path filesep 'df_f' filesep Field];
        Target_path = [General_path filesep 'df_f' filesep Field];
        df_f_list = dir(df_f_folder);
        df_f_list = {df_f_list.name};
        df_f_list = sort(df_f_list);
        df_f_list = df_f_list(3:end-1);       
        for curr_day = 1:length(df_f_list)
            curr_df_f = df_f_list{curr_day};
            disp(curr_df_f);            
            load([Target_path filesep curr_df_f],'roi_trace','roi_trace_bg','-mat');
            Date = curr_df_f(14:19);
            for curr_session = 1:length(roi_trace)
                frame_num = size(roi_trace{curr_session},2);
                if ismember(Field, {'Field_1','Field_3'})
                    frame_index = [1:2:frame_num];
                elseif ismember(Field, {'Field_2','Field_4'})
                    frame_index = [2:2:frame_num];
                else
                    disp('error!');
                end
                roi_trace{curr_session} = roi_trace{curr_session}(:,frame_index);
                roi_trace_bg{curr_session} = roi_trace_bg{curr_session}(:,frame_index);
                roi_trace_long{curr_session} = roi_trace{curr_session} - roi_trace_bg{curr_session};
                % *** not used, re-estimate later ***
                framerate = 14;
                [roi_trace_df{curr_session}, roi_trace_baseline{curr_session}] = AP_baselineEstimation(roi_trace_long{curr_session},framerate);    
                % ***********************************
                if size(roi_trace_df{curr_session},1) < 5
                    selected_roi_index = [1:size(roi_trace_df{curr_session},1)];
                elseif size(roi_trace{curr_session},1) < 10                    
                    selected_roi_index = [1:5];
                elseif size(roi_trace{curr_session},1) < 20                    
                    selected_roi_index = [1:2:10];
                else
                    selected_roi_index = [1:2:20];
                end
                roi_for_qualitycheck{curr_day}{curr_session} = roi_trace{curr_session}(selected_roi_index,:);
            end
            Target_filename = [Animal '_' Date '_ROI_Traces.mat'];
            save([Target_path filesep Target_filename],'roi_trace','roi_trace_bg','roi_trace_df','roi_trace_baseline','-v7.3');
            clear roi_trace roi_trace_bg roi_trace_df roi_trace_baseline
        end
        save([General_path filesep 'df_f' filesep Field filesep Animal '_' Field '_ROI_Traces_Check.mat'],'roi_for_qualitycheck','-v7.3');
        clear roi_for_qualitycheck
    end
end

%% ************************************************************************************* %%
%% Check truncate point, MOM has off edge artifact
clear all;
close all;
clc;

IN = 'VIP';
Animals = {'WL_3526642-R'};
Rig = 'MOM';

for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'df_f']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(cellfun(@(x) ~isempty(strfind(x, 'Field')), Fields));
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI')) && isempty(strfind(x, 'Check')), Dates));
        Dates = cellfun(@(x) x(14:19), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        load([General_path filesep 'df_f' filesep Field filesep Animal '_' Field '_ROI_Traces_Check.mat']);
        for curr_day = 1:length(Dates)
            for curr_session = 1:length(roi_for_qualitycheck{curr_day})
                window_start = 1:5000;
                window_end = size(roi_for_qualitycheck{curr_day}{curr_session},2)-5000:size(roi_for_qualitycheck{curr_day}{curr_session},2);
                figure; set(gcf,'position',[400,400,700,400]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck{curr_day}{curr_session},1)
                    plot(roi_for_qualitycheck{curr_day}{curr_session}(curr_roi,window_start)-1000*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session) ' Start']);
                [temp_start,~] = ginput(1);
                temp_start = round(temp_start);
                close(gcf);
                figure; set(gcf,'position',[400,400,700,400]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck{curr_day}{curr_session},1)
                    plot(roi_for_qualitycheck{curr_day}{curr_session}(curr_roi,window_end)-1000*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session) ' End']);
                [temp_end,~] = ginput(1);
                temp_end = temp_end + size(roi_for_qualitycheck{curr_day}{curr_session},2) - 5000;
                temp_end = round(temp_end);
                close(gcf);
                figure; set(gcf,'position',[400,400,700,400]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck{curr_day}{curr_session},1)
                    plot(roi_for_qualitycheck{curr_day}{curr_session}(curr_roi,:)-1000*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session)]);
                [temp_all,~] = ginput(2);
                temp_all = round(temp_all);
                close(gcf)
                % start
                truncatePoint{curr_session}(1) = max([temp_start,temp_all(1),1]);
                % end
                truncatePoint{curr_session}(2) = min([temp_end,temp_all(2),size(roi_for_qualitycheck{curr_day}{curr_session},2)]);
                if truncatePoint{curr_session}(2) == size(roi_for_qualitycheck{curr_day}{curr_session},2)
                    istruncate{curr_session} = true;
                    truncatePoint{curr_session}(2) = truncatePoint{curr_session}(2)-20; % last 20 data points are often bad;
                else
                    istruncate{curr_session} = true;
                end
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'truncatePoint','istruncate','-append');
            clear truncatePoint istruncate temp_start temp_end temp_all
        end
    end
end

%% Baseline f estimation, still AP method, change baseline smooth window to 4 min
clear all;

IN = 'VIP';
Animals = {'WL_3526642-R'};
Rig = 'MOM';

df_threshold = 50;
for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal df_threshold
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'df_f']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(cellfun(@(x) ~isempty(strfind(x, 'Field')), Fields));
    load([Animal '_ImagingInfo.mat']);
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI_Traces')) && isempty(strfind(x, 'Check')), Dates));
        if ismember(Animal,{'KP_3475729_LR','CR_3619106-LR','WL_3547273-LR','CR_4042831-LR','CR_4017421-LR'})
            Dates = cellfun(@(x) x(15:20), Dates, 'UniformOutput', false);
        else
            Dates = cellfun(@(x) x(14:19), Dates, 'UniformOutput', false);
        end
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace roi_trace_bg roi_trace_baseline roi_trace_df roi_trace_df_2 roi_trace_baseline_2 roi_trace_long    
            clear roi_trace_origin roi_trace_bg_origin roi_trace roi_trace_bg roi_trace_df roi_trace_baseline truncatePoint                
            clear istruncate
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'-mat');
            if ~exist('roi_trace_origin')
                disp('Saving original traces...');
                roi_trace_origin = roi_trace;
                roi_trace_bg_origin = roi_trace_bg;
            end
            for curr_session = 1:length(roi_trace) 
                disp([Animal ' Field' num2str(curr_field) ' ' Dates{curr_day} ' Session' num2str(curr_session)]);
                if truncatePoint{curr_session}(2)-truncatePoint{curr_session}(1)+1 ~= size(roi_trace{curr_session},2)
                    roi_trace{curr_session} = roi_trace_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2));
                    roi_trace_bg{curr_session} = roi_trace_bg_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2)); 
                end
                roi_trace_long{curr_session} = roi_trace{curr_session} - roi_trace_bg{curr_session};
                framerate = Imaging_Fields{curr_field}.FrameRate(curr_day);                
                roi_trace_df_2{curr_session} = nan(size(roi_trace_long{curr_session}));
                roi_trace_baseline_2{curr_session} = nan(size(roi_trace_long{curr_session}));
                for roi = 1:size(roi_trace_long{curr_session},1)
                    temp_trace = roi_trace_long{curr_session}(roi,:);
                    nan_index = ~isnan(temp_trace);
                    if sum(~nan_index)
                        disp([num2str(roi) ' nan detected']);
                        if sum(~nan_index) == length(temp_trace)
                            disp([num2str(roi) ' all nan, ignore']);
                            continue
                        end
                    end
                    [temp_df,temp_baseline] = CR_AP_baselineEstimation(temp_trace(nan_index),framerate,0,4);
                    roi_trace_df_2{curr_session}(roi,nan_index) = temp_df;
                    roi_trace_baseline_2{curr_session}(roi,nan_index) =temp_baseline;
                    clear temp_trace temp_df temp_baseline
                end
                
                % get rid of noise
                roi_trace_df_2{curr_session}(abs(roi_trace_df_2{curr_session}) >= df_threshold) = nan;
                % Plot for post check
                FigTargetPath = [General_path filesep 'df_f' filesep 'TraceToCheck' filesep Field];
                if ~exist(FigTargetPath)
                    mkdir(FigTargetPath)
                end
                close all;
                figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 900]);
                hold on;
                for curr_roi = 1:size(roi_trace_long{curr_session},1)
                    tt = [1:length(roi_trace_long{curr_session})]/framerate;
                    plot(tt,roi_trace_df_2{curr_session}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
                end
                set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');
                savefig(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) '.fig']);
                saveas(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) '.png']);
                pause(0.1);
                close(gcf);
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_origin','roi_trace_bg_origin','roi_trace','roi_trace_bg','roi_trace_df_2','roi_trace_baseline_2','-append'); 
        end
    end
end

%% Please eyeball traces and nan bad values (usually show huge noise due to dim signal)
% WL_3526642-R
clear all
IN = 'VIP';
Animal = 'WL_3526642-R';
RejectedROI.Field_1 = [14];
RejectedROI.Field_2 = [];
RejectedROI.Field_3 = [9,23,43,53,70,90,95,103];
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f' filesep Animal '_ImagingInfo.mat'];
save(TargetPath,'RejectedROI','-append');

%% Reject bad cells
clear all
close all
clc

IN = 'VIP';
Animals = {'WL_3526642-R'};

for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal df_threshold
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'df_f']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(cellfun(@(x) ~isempty(strfind(x, 'Field')), Fields));
    load([Animal '_ImagingInfo.mat'],'-mat');
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI_Traces')) && isempty(strfind(x, 'Check')), Dates));
        if ismember(Animal,{'KP_3475729_LR','WL_3547273-LR','CR_4042831-LR','CR_4017421-LR'})
            Dates = cellfun(@(x) x(15:20), Dates, 'UniformOutput', false);
        else
            Dates = cellfun(@(x) x(14:19), Dates, 'UniformOutput', false);
        end
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace_baseline_2 roi_trace_df_2 roi_trace_df_2_correction truncatePoint istruncate
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_df_2','roi_trace_baseline_2','truncatePoint','istruncate','-mat');
            for curr_session = 1:length(roi_trace_df_2) 
                disp([Animal ' Field' num2str(curr_field) ' ' Dates{curr_day} ' Session' num2str(curr_session)]);
                disp(['Rejected ROI: ' num2str(RejectedROI.(Field))]);
                roi_trace_df_2{curr_session}(RejectedROI.(Field),:) = [];
                roi_trace_baseline_2{curr_session}(RejectedROI.(Field),:) = [];
                roi_trace_df_2_correction{curr_session} = true;
                disp(['Original truncatePoint: ' num2str(truncatePoint{curr_session})]);
                truncatePoint{curr_session}(2) = size(roi_trace_df_2{curr_session},2) + truncatePoint{curr_session}(1) - 1;
                disp(['Updated truncatePoint: ' num2str(truncatePoint{curr_session})]);
                % figure
                FigTargetPath = [General_path filesep 'df_f' filesep 'TraceToCheck' filesep Field];
                if ~exist(FigTargetPath)
                    mkdir(FigTargetPath)
                end
                framerate = Imaging_Fields{curr_field}.FrameRate(curr_day);
                close all;
                figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 900]);
                hold on;
                for curr_roi = 1:size(roi_trace_df_2{curr_session},1)
                    tt = [1:length(roi_trace_df_2{curr_session})]/framerate;
                    plot(tt,roi_trace_df_2{curr_session}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
                end
                set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');
                savefig(gcf,[FigTargetPath filesep 'New_' Field '_' Dates{curr_day} '_' num2str(curr_session) '.fig']);
                saveas(gcf,[FigTargetPath filesep 'New_' Field '_' Dates{curr_day} '_' num2str(curr_session) '.png']);
                pause(0.1);
                close(gcf);
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_df_2','roi_trace_baseline_2','truncatePoint','istruncate','roi_trace_df_2_correction','-append'); 
        end
    end
end

%% Get CaEvents and z-scored traces
clear all
close all
clc

IN = 'VIP';
Animals = {'WL_3526642-R'};

for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal df_threshold
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal];
    cd([General_path filesep 'df_f']);
    Fields = dir(cd);
    Fields = {Fields.name};
    Fields = sort(Fields);
    Fields = Fields(cellfun(@(x) ~isempty(strfind(x, 'Field')), Fields));
    load([Animal '_ImagingInfo.mat'],'-mat');
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI_Traces')) && isempty(strfind(x, 'Check')), Dates));
        if ismember(Animal,{'KP_3475729_LR','WL_3547273-LR','CR_4042831-LR','CR_4017421-LR'})
            Dates = cellfun(@(x) x(15:20), Dates, 'UniformOutput', false);
        else
            Dates = cellfun(@(x) x(14:19), Dates, 'UniformOutput', false);
        end
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace_df_2 CaEvents_2 ZScore_2
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_df_2','-mat');
            for curr_session = 1:length(roi_trace_df_2) 
                disp([Animal ' Field' num2str(curr_field) ' ' Dates{curr_day} ' Session' num2str(curr_session)]);
                thresh = 2;
                method = 2;
                framerate = Imaging_Fields{curr_field}.FrameRate(curr_day);
                noise_est_window_min = 6;
                figureOn = false;
                CaEvents_2{curr_session} = CR_AP_caEvents_thresh(roi_trace_df_2{curr_session},thresh,method,framerate,noise_est_window_min,figureOn);
                if any(isnan(roi_trace_df_2{curr_session}(:)))
                    temp_matrix = bsxfun(@minus, roi_trace_df_2{curr_session}', nanmean(roi_trace_df_2{curr_session}'));
                    ZScore_2{curr_session} = bsxfun(@rdivide, temp_matrix, nanstd(roi_trace_df_2{curr_session}'));
                    ZScore_2{curr_session} = ZScore_2{curr_session}';
                    clear temp_matrix
                else
                    ZScore_2{curr_session} = zscore(roi_trace_df_2{curr_session},[],2);
                end 
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'CaEvents_2','ZScore_2','-append'); 
        end
    end
end


