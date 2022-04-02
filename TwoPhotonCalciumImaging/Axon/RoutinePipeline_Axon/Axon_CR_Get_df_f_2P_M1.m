%% Get roi traces & df/f, for M1
clear all;
close all;
clc;

IN = 'ChAT';
Animals = {'CR_3672020-R'};
Rig = 'BScope1';
Dates = {'190511'};

% Get axon roi and background roi
for curr_animal = 1:length(Animals)
    clearvars -except IN Animals curr_animal Rig Dates
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal);
    Local_path = fullfile('F:\Data\MotionCorrection\',[Animal '_warp']);
       
    for curr_day = 1:length(Dates)
        Date = Dates{curr_day};
        disp(Date);            
        Target_path = [General_path filesep 'df_f' filesep 'Field_1'];
        if ~exist(Target_path)
            mkdir(Target_path);
        end
        
        image_folder_path = fullfile(Local_path,Date);
        tiff_folders = dir(image_folder_path);
        tiff_folders = {tiff_folders.name};
        tiff_folders = sort(tiff_folders);
        tiff_folders = tiff_folders(3:end);
        for curr_session = 1:length(tiff_folders)
            roi_filename = fullfile(image_folder_path,tiff_folders{curr_session},'Z1','motioncorrected_tiff','summed',[Animal(1:3) Date Animal(3:end) '_roi.mat']);
            bg_roi_file = fullfile(image_folder_path,tiff_folders{curr_session},'Z1','motioncorrected_tiff','summed',[Animal(1:3) Date Animal(3:end) '_Rec' num2str(curr_session) '_00001_00001_summed_50_warped_bg.roi']);    
            tiff_path = fullfile(image_folder_path,tiff_folders{curr_session},'Z1','motioncorrected_tiff');
            interp_sub = 0;
            local_comp = 0;
            switch Rig
                case 'MOM'
                    framerate = 28; 
                    [roi_trace{curr_session},roi_trace_bg{curr_session}] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp);                 
                case 'BScope1'
                    framerate = 30;
                    [roi_trace{curr_session},roi_trace_bg{curr_session},polygon{curr_session}] = CR_Suite2P_Axon_Rig3_getConcatTrace_continuous_batch(roi_filename,bg_roi_file,tiff_path,interp_sub,local_comp);                        
            end
                
        end
        Target_filename = [Animal '_' Date '_ROI_Traces.mat'];
        save([Target_path filesep Target_filename],'roi_trace','roi_trace_bg','polygon','-v7.3');
        clear roi_trace roi_trace_bg polygon
    end

end

%% Check truncate point 
clear all;
close all;
clc;

IN = 'ChAT';
Animals = {'CR_3672020-R'};

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
        cd(Field);
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI')) && isempty(strfind(x, 'Check')), Dates));
        Dates = cellfun(@(x) x(length(Animal)+2:length(Animal)+7), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            Date = Dates{curr_day};
            load([Animal '_' Date '_ROI_Traces.mat'],'-mat');            

            for curr_session = 1:length(roi_trace)
                if size(roi_trace{curr_session},1) < 10                    
                    selected_roi_index = [1:5];
                elseif size(roi_trace{curr_session},1) < 20                    
                    selected_roi_index = [1:2:10];
                else
                    selected_roi_index = [1:2:20];
                end
                roi_for_qualitycheck = roi_trace{curr_session}(selected_roi_index,:);
                window_start = 1:5000;
                window_end = size(roi_for_qualitycheck,2)-5000:size(roi_for_qualitycheck,2);
                figure; set(gcf,'position',[0,0,1200,600]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck,1)
                    plot(roi_for_qualitycheck(curr_roi,window_start)-1000*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session) ' Start']);
                [temp_start,~] = ginput(1);
                temp_start = round(temp_start);
                close(gcf);
                figure; set(gcf,'position',[50,50,1400,800]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck,1)
                    plot(roi_for_qualitycheck(curr_roi,window_end)-500*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session) ' End']);
                [temp_end,~] = ginput(1);
                temp_end = temp_end + size(roi_for_qualitycheck,2) - 5000;
                temp_end = round(temp_end);
                close(gcf);
                figure; set(gcf,'position',[50,50,1400,800]); hold on;
                for curr_roi = 1:size(roi_for_qualitycheck,1)
                    plot(roi_for_qualitycheck(curr_roi,:)-1000*(curr_roi-1));
                end
                title([Animal ' ' Field ' ' num2str(curr_day) '/' num2str(length(Dates)) ' ' num2str(curr_session)]);
                [temp_all,~] = ginput(2);
                temp_all = round(temp_all);
                close(gcf)
                % start
                truncatePoint{curr_session}(1) = max([temp_start,temp_all(1),1]);
                % end
                truncatePoint{curr_session}(2) = min([temp_end,temp_all(2),size(roi_for_qualitycheck,2)]);
                if truncatePoint{curr_session}(2) == size(roi_for_qualitycheck,2)
                    istruncate{curr_session} = true;
                    truncatePoint{curr_session}(2) = truncatePoint{curr_session}(2)-20; % last 20 data points are often bad;
                else
                    istruncate{curr_session} = true;
                end
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'truncatePoint','istruncate','-append');
            clear roi_trace roi_trace_bg polygon roi_for_qualitycheck truncatePoint istruncate temp_start temp_end temp_all
        end
    end
end

%% Two baseline estimation methods: 1) 30% in 4 min moving window; 2) Still AP method, change baseline smooth window to 4 min
clear all;

IN = 'ChAT';
Animals = {'CR_3672020-R'};
Rig = 'BScope1';

df_threshold = 100;
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
        Dates = cellfun(@(x) x(length(Animal)+2:length(Animal)+7), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace_origin roi_trace_bg_origin roi_trace roi_trace_bg roi_trace_baseline roi_trace_df roi_trace_df_2 roi_trace_baseline_2 roi_trace_long    
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'-mat');
            if ~exist('roi_trace_origin')
                disp('Saving original traces...');
                roi_trace_origin = roi_trace;
                roi_trace_bg_origin = roi_trace_bg;
            end
            for curr_session = 1:length(roi_trace) 
                disp([Animal ' Field' num2str(curr_field) ' ' Dates{curr_day} ' Session' num2str(curr_session)]);
                roi_trace{curr_session} = roi_trace_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2));
                roi_trace_bg{curr_session} = roi_trace_bg_origin{curr_session}(:,truncatePoint{curr_session}(1):truncatePoint{curr_session}(2));
                roi_trace_bg{curr_session} = nanmean(roi_trace_bg{curr_session},1);
                roi_trace_long{curr_session} = roi_trace{curr_session} - repmat(roi_trace_bg{curr_session},size(roi_trace{curr_session},1),1);
                framerate = Imaging_Fields{curr_field}.FrameRate(curr_day);
                roi_trace_df_2{curr_session} = nan(size(roi_trace_long{curr_session}));
                roi_trace_baseline_2{curr_session} = nan(size(roi_trace_long{curr_session}));
                roi_trace_df_3{curr_session} = nan(size(roi_trace_long{curr_session}));
                roi_trace_baseline_3{curr_session} = nan(size(roi_trace_long{curr_session}));
                all_nan_index = sum(isnan(roi_trace_long{curr_session}),2);
                baseline_window = round(framerate*2*60);
                prctile_thre = 30;
                tic
                for frame = 1:size(roi_trace_long{curr_session},2)
                    if frame <= baseline_window
                        roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,1:frame+baseline_window),prctile_thre,2);
                    elseif frame >= size(roi_trace_long{curr_session},2)-baseline_window
                        roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,frame-baseline_window:size(roi_trace_long{curr_session},2)),prctile_thre,2);
                    else
                        roi_trace_baseline_2{curr_session}(:,frame) = prctile(roi_trace_long{curr_session}(:,frame-baseline_window:frame+baseline_window),prctile_thre,2);
                    end
                end
                if any(all_nan_index)
                    disp('Calculating for nan containing traces');
                    nan_rois = find(all_nan_index);
                    for roi = nan_rois
                        temp_trace = roi_trace_long{curr_session}(roi,:);
                        nan_index = ~isnan(temp_trace);
                        if sum(~nan_index) == length(temp_trace)
                            disp([num2str(roi) ' all nan, ignore']);
                            continue
                        end
                        temp_trace = temp_trace(nan_index);
                        clear temp_baseline
                        for frame = 1:length(temp_trace)
                            if frame <= baseline_window
                                temp_baseline(frame) = prctile(temp_trace(1:frame+baseline_window),prctile_thre);
                            elseif frame >= size(roi_trace_long{curr_session},2)-baseline_window
                                temp_baseline(frame) = prctile(temp_trace(frame-baseline_window:size(roi_trace_long{curr_session},2)),prctile_thre);
                            else
                                temp_baseline(frame) = prctile(temp_trace(frame-baseline_window:frame+baseline_window),prctile_thre);
                            end
                        end
                        roi_trace_baseline_2{curr_session}(roi,nan_index) = temp_baseline;
                    end
                end
                toc
                roi_trace_df_2{curr_session} = (roi_trace_long{curr_session}-roi_trace_baseline_2{curr_session})./roi_trace_baseline_2{curr_session};
                % get rid of noise
                roi_trace_df_2{curr_session}(abs(roi_trace_df_2{curr_session}) >= df_threshold) = nan;
                % AP baseline
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
                    [temp_df,temp_baseline] = CR_AP_baselineEstimation_Axon(temp_trace(nan_index),framerate,0,4);
                    roi_trace_df_3{curr_session}(roi,nan_index) = temp_df;
                    roi_trace_baseline_3{curr_session}(roi,nan_index) =temp_baseline;
                    clear temp_trace temp_df temp_baseline
                end
                % get rid of noise
                roi_trace_df_3{curr_session}(abs(roi_trace_df_3{curr_session}) >= df_threshold) = nan;
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
%                 saveas(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) '.png']);
                pause(0.1);
                close(gcf);
                
                figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 900]);
                hold on;
                for curr_roi = 1:size(roi_trace_long{curr_session},1)
                    tt = [1:length(roi_trace_long{curr_session})]/framerate;
                    plot(tt,roi_trace_df_3{curr_session}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
                end
                set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');
                savefig(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) '_AP.fig']);
%                 saveas(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) '.png']);
                pause(0.1);
                close(gcf);
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_origin','roi_trace_bg_origin','roi_trace_df_2','roi_trace_baseline_2','roi_trace_df_3','roi_trace_baseline_3','-append'); 
        end
    end
end

%% Reject bad ROIs (if too dim or out of focal plain in the middle)
load('CR_3672020-R_190511_ROI_Traces.mat');
Rejected_ROIs = [];
roi_trace_df_2{1}(Rejected_ROIs,:) = [];
figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 800]);
hold on;
for curr_roi = 1:size(roi_trace_df_2{1},1)
    tt = [1:length(roi_trace_df_2{1})]/30.3;
    plot(tt,roi_trace_df_2{1}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
end
set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');

close all
save('CR_3672020-R_190511_ROI_Traces.mat','roi_trace_df_2','Rejected_ROIs','-append');
clear


%% Re plot
clear all;

IN = 'ChAT';
Animals = {'CR_3672020-R'};

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
        Dates = cellfun(@(x) x(length(Animal)+2:length(Animal)+7), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace_df_2 Rejected_ROIs roi_trace_origin
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_df_2','roi_trace_origin','Rejected_ROIs','-mat');
            framerate = Imaging_Fields{curr_field}.FrameRate(curr_day); 
            for curr_session = 1:length(roi_trace_df_2) 
                
                FigTargetPath = [General_path filesep 'df_f' filesep 'TraceToCheck' filesep Field];
                if ~exist(FigTargetPath)
                    mkdir(FigTargetPath)
                end
                close all;
                figure; set(gcf,'color',[0.9,0.9,0.9],'position',[50 50 1200 900]);
                hold on;
                for curr_roi = 1:size(roi_trace_df_2{curr_session},1)
                    tt = [1:length(roi_trace_df_2{curr_session})]/framerate;
                    plot(tt,roi_trace_df_2{curr_session}(curr_roi,:)-5*(curr_roi-1)); xlim([1 tt(end)]);    
                end
                set(gca,'color',[0.9,0.9,0.9]); box off; ylabel('df/f'); xlabel('Time (sec)');
                savefig(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) 'New.fig']);
                saveas(gcf,[FigTargetPath filesep Field '_' Dates{curr_day} '_' num2str(curr_session) 'New.png']);
                pause(0.1);
                close(gcf);                              
            end
        end
    end
end


%% Check after rejection
clear all;

IN = 'ChAT';
Animals = {'CR_3633170-R'};

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
    for curr_field = 1:length(Fields)
        Field = Fields{curr_field};
        Dates = dir([General_path filesep 'df_f' filesep Field]);
        Dates = {Dates.name};
        Dates = Dates(cellfun(@(x) ~isempty(strfind(x, 'ROI_Traces')) && isempty(strfind(x, 'Check')), Dates));
        Dates = cellfun(@(x) x(length(Animal)+2:length(Animal)+7), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear Rejected_ROIs
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'-mat');
            for curr_session = 1:length(roi_trace_df_2) 
                if exist('Rejected_ROIs')
                    if size(roi_trace_df_2{curr_session},1)+length(Rejected_ROIs)~=size(roi_trace{curr_session},1)
                        disp([Animal '_' Dates{curr_day} '_ROI_Traces.mat']);
                    end
                else
                   if size(roi_trace_df_2{curr_session},1)~=size(roi_trace{curr_session},1)
                        disp([Animal '_' Dates{curr_day} '_ROI_Traces.mat']);
                   end 
                end
            end
        end
    end
end
                
%% Get z-scored traces
clear all
close all
clc

IN = 'ChAT';
Animals = {'CR_3672020-R'};

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
        Dates = cellfun(@(x) x(length(Animal)+2:length(Animal)+7), Dates, 'UniformOutput', false);
        Dates = sort(Dates);
        for curr_day = 1:length(Dates)
            clear roi_trace_df_2 CaEvents_2 ZScore_2
            load([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'roi_trace_df_2','-mat');
            for curr_session = 1:length(roi_trace_df_2) 
                disp([Animal ' Field' num2str(curr_field) ' ' Dates{curr_day} ' Session' num2str(curr_session)]);
                if any(isnan(roi_trace_df_2{curr_session}(:)))
                    temp_matrix = bsxfun(@minus, roi_trace_df_2{curr_session}', nanmean(roi_trace_df_2{curr_session}'));
                    ZScore_2{curr_session} = bsxfun(@rdivide, temp_matrix, nanstd(roi_trace_df_2{curr_session}'));
                    ZScore_2{curr_session} = ZScore_2{curr_session}';
                    clear temp_matrix
                else
                    ZScore_2{curr_session} = zscore(roi_trace_df_2{curr_session},[],2);
                end 
            end
            save([General_path filesep 'df_f' filesep Field filesep Animal '_' Dates{curr_day} '_ROI_Traces.mat'],'ZScore_2','-append'); 
        end
    end
end

