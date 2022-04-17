%% 2P Data ChAT
%% Averaged activity after movement onset, combine trials across sessions first then average
clear all;
close all;
clc;

IN = 'ChAT';
% M1
Animals = {'CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R'};
Fields = {'Field_1'};

for curr_animal = 1:length(Animals)
    clearvars -except Animals Drug IN curr_animal Fields
    Animal = Animals{curr_animal};
    LeverTrace_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'LeverTrace');
    Dates = dir(LeverTrace_Path);
    Dates = {Dates.name};
    Dates = sort(Dates(cellfun(@(x) contains(x, '17')||contains(x, '18')||contains(x, '19')||contains(x, '21'), Dates)))';
    % load imaging info
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'df_f',[Animal '_ImagingInfo.mat']),'Imaging_Fields','-mat');
    % load mov-onset frame index
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'MovAnalysis',[Animal '_CuedRewardedMov_FrameIndex.mat']),'CuedRewardedFrameIndex','CuedRewardedFrameIndex_z1','CuedRewardedFrameIndex_z2','-mat');
    for curr_field = 1:length(Imaging_Fields)
        Field = Fields{curr_field};
        Field_Dates = Imaging_Fields{curr_field}.Date;
        for curr_date = 1:length(Field_Dates)
            Date = Field_Dates{curr_date};
            df_f_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'df_f',Field);
            date_index = ismember(Dates,Date);
            FrameOnset_index = CuedRewardedFrameIndex{date_index};

            % load df/f
            if ~exist([df_f_Path filesep Animal '_' Date '_ROI_Traces.mat'])
                df_f_MovOnset_Alinged{curr_field,curr_date} = [];
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date} = [];
                continue
            end
            load([df_f_Path filesep Animal '_' Date '_ROI_Traces.mat'],'roi_trace_df_2','truncatePoint','-mat');

            for curr_session = 1:length(FrameOnset_index)
                disp([Animal ' ' Field '_' Date ' Session' num2str(curr_session)]);
                end_control = FrameOnset_index{curr_session}(end,:) <= round(truncatePoint{curr_session}(2));
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session}(:,end_control);
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session} - round(truncatePoint{curr_session}(1))+1;
                start_control = FrameOnset_index{curr_session}(1,:) > 0;
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session}(:,start_control);
                for roi = 1:size(roi_trace_df_2{curr_session},1)
                    temp_trace = roi_trace_df_2{curr_session}(roi,:);
                    for trial = 1:size(FrameOnset_index{curr_session},2)
                        df_f_MovOnset_Alinged{curr_field,curr_date}{curr_session,1}{roi,1}(:,trial) = temp_trace(FrameOnset_index{curr_session}(:,trial));
                    end
                    clear temp_trace temp_CaEvents temp_ZScore
                end           
            end
            if size(df_f_MovOnset_Alinged{curr_field,curr_date},1)>1
                temp = [df_f_MovOnset_Alinged{curr_field,curr_date}{1,1};df_f_MovOnset_Alinged{curr_field,curr_date}{2,1}];
                df_f_MovOnset_Alinged{curr_field,curr_date} = temp;
   
            else
                df_f_MovOnset_Alinged{curr_field,curr_date} = df_f_MovOnset_Alinged{curr_field,curr_date}{1};
            end
            for roi = 1:size(df_f_MovOnset_Alinged{curr_field,curr_date},1)
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);
            end
        end        
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    save([TargetPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','df_f_MovOnset_Alinged_EachROI_MeanAXTrial','-v7.3');
end

%%
clear all;
close all;
clc;

IN = 'ChAT';
Animals = {'CR_3619073-R','CR_3619073-LR','CR_3619074-O','CR_3633170-L','CR_3672020-O','CR_3672020-L','CR_3672020-R'};
Fields = {'Field_1'};

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','-mat');
    for curr_field = 1:size(df_f_MovOnset_Alinged,1)
        for curr_session = 1:size(df_f_MovOnset_Alinged(curr_field,:),2)
            if strcmp(Animal,'CR_3672020-L') && curr_session == 1
                Post_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session} = [];
                Post_MovOnset_Aligend_df_field_sub{curr_animal}{curr_field,curr_session} = [];
                continue
            end
            if isempty(df_f_MovOnset_Alinged{curr_field,curr_session})
                Post_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session} = [];
                Post_MovOnset_Aligend_df_field_sub{curr_animal}{curr_field,curr_session} = [];
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
                Post_MovOnset_Aligend_df_field_sub{curr_animal}{curr_field,curr_session}(curr_roi,:) = nanmean(temp(MovOnset_Frame:end,:)-...
                    repmat(nanmean(temp(Baseline_Frame,:),1),size(temp(MovOnset_Frame:end,:),1),1),1);
                Post_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session}(curr_roi,:) = nanmean(df_f_MovOnset_Alinged{curr_field,curr_session}{curr_roi}(MovOnset_Frame:end,:),1);
            end
            AxonNum(curr_animal,curr_session) = size(df_f_MovOnset_Alinged{curr_field,curr_session},1);
        end
        % 4 stages
        temp = cellfun(@(x) nanmean(x,2),Post_MovOnset_Aligend_df_field{curr_animal},'UniformOutput',false);        
        temp_size = cellfun(@(x) size(x,2),Post_MovOnset_Aligend_df_field{curr_animal},'UniformOutput',false);
        temp(cell2mat(temp_size)<3) = {[]};
        temp = temp';
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1} = temp(1:2);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,2} = temp(3:4);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,3} = temp(5:6);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,4} = temp(7:end);
        
        temp = cellfun(@(x) nanmean(x,2),Post_MovOnset_Aligend_df_field_sub{curr_animal},'UniformOutput',false);
        temp_size = cellfun(@(x) size(x,2),Post_MovOnset_Aligend_df_field_sub{curr_animal},'UniformOutput',false);
        temp(cell2mat(temp_size)<3) = {[]};
        temp = temp';
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,1} = temp(1:2);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,2} = temp(3:4);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,3} = temp(5:6);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,4} = temp(7:end);
                
    end
    % Average
    for ii = 1:4 % average within each session first
        Post_MovOnset_Aligend_df_stage_Mean(curr_animal,ii) = nanmean(cellfun(@mean, Post_MovOnset_Aligend_df_field_stage{curr_animal}{1,ii}));
        Post_MovOnset_Aligend_df_stage_sub_Mean(curr_animal,ii) = nanmean(cellfun(@mean, Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{1,ii}));
    end
    for ii = 1:4 % pool across sessions then average
        Post_MovOnset_Aligend_df_stage_Mean_2(curr_animal,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field_stage{curr_animal}{1,ii}));
        Post_MovOnset_Aligend_df_stage_sub_Mean_2(curr_animal,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{1,ii}));
    end
    clear df_f_MovOnset_Alinged ZScore_MovOnset_Alinged
end

[CX] = cbrewer('div','RdBu',10);
color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];   
color_line = nanmean(color_value(2:3,:));
stages = {'Naive','Early','Middle','Late'};

figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = Post_MovOnset_Aligend_df_stage_sub_Mean_2;
for curr_animal = 1:size(temp_var,1)
    plot(temp_var(curr_animal,:),'color',[0.5 0.5 0.5],'LineWidth',0.5);
end
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,1]); xticklabels(stages); ylabel('Mean df/f');
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,4);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.8,'p < 0.01','color','k','fontsize',8,'HorizontalAlignment','center');
xlabel('Stages'); ylabel('Mean df/f');
axis square

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep 'Figures_SubPre'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
saveas(gcf, [FigurePath filesep IN '_Mean_dfof.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof.pdf']); pause(1);

close all

% Heatmap
% 4 stages
clear df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged_EachROI_MeanAXTrial');
    if strcmp(Animal,'CR_3672020-L')
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,1} = [];
    end
    if strcmp(Animal,'CR_3633170-L')
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,1} = [];
    end
        
    for curr_field = 1:size(df_f_MovOnset_Alinged_EachROI_MeanAXTrial,1)
        for curr_session = 1:size(df_f_MovOnset_Alinged_EachROI_MeanAXTrial,2)
            if isempty(df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_session})
                temp_matrix{curr_session,1} = [];
                continue
            end
            temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_session};
            temp = temp-repmat(nanmean(temp(:,8:11),2),1,76);
            temp_matrix{curr_session,1} = temp; 
            clear temp
        end
    end
    df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,1} = cell2mat(temp_matrix(1:2,:));
    df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,2} = cell2mat(temp_matrix(3:4,:));
    df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,3} = cell2mat(temp_matrix(5:6,:));
    df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,4} = cell2mat(temp_matrix(7:12,:));
    clear temp_matrix
end

% Sort based on activity
for curr_animal = 1:length(Animals)
    for ii = 1:4
        temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,ii};
        mean_value = nanmean(temp,2);
        [~,I] = sort(mean_value,'descend');
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_sort{curr_animal,ii} = temp(I,:);
        AxonNum_stage(curr_animal,ii) = size(temp,1);
        clear temp
    end
end

figure; hold on; set(gcf,'color','w','position',[50 50 800 500]);
Stages = {'Naive','Early','Middle','Late'};
for ii = 1:4
    subplot(1,4,ii);
    imagesc(cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_sort(:,ii)),[-0.5 2.2]);
    line([16,16],ylim,'color','w','linestyle',':');
    cellnum = cumsum(AxonNum_stage(:,ii));
    for jj = 1:length(cellnum)-1
        line(xlim,[cellnum(jj)-0.5, cellnum(jj)-0.5],'color','k','linestyle',':');
    end
    ylabel('# of cells');
    xticks([16,46,76]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Stages{ii});
end
saveas(gcf, [FigurePath filesep IN '_MeanXTrial_dfof_conc.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_MeanXTrial_dfof_conc.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_MeanXTrial_dfof_conc.pdf']); pause(1);

% Average trace
clear temp_mean
for curr_animal = 1:length(Animals)
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,1};
    temp_mean{1}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,2};
    temp_mean{2}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,3};
    temp_mean{3}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,4};
    temp_mean{4}(curr_animal,:) = nanmean(temp,1);
end
% across animals
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = temp_mean{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
xlim([1 76]); ylim([-0.1 0.7])
xticks([16,46,76]);xticklabels({'0','1','2'});
xlabel('Time (sec)');
ylabel('Mean df/f');
line([16 16], ylim, 'color','k','linestyle',':','linewidth',1);
axis square;

saveas(gcf, [FigurePath filesep IN '_Mean_dfof_temporal.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_temporal.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_temporal.pdf']); pause(1);

save([FigurePath filesep IN '_ActivityAnalysis_SubPre_M1.mat'],'-append');

