%% 2P Data VIP
clear all;
close all;
clc;

IN = 'VIP_SOM';
Animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R',...
        'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
Fields = {'Field_1','Field_2','Field_3','Field_4'};
Stage = 'Nai';

for curr_animal = 1:length(Animals)
    clearvars -except Animals IN curr_animal Fields Stage
    Animal = Animals{curr_animal};
    LeverTrace_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,Stage,'LeverTrace');
    Dates = dir(LeverTrace_Path);
    Dates = {Dates.name};
    Dates = sort(Dates(cellfun(@(x) contains(x, '17')||contains(x, '18')||contains(x, '19')||contains(x, '20'), Dates)))';
    % load imaging info
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,Stage,'df_f',[Animal '_ImagingInfo.mat']),'Imaging_Fields','-mat');
    % load mov-onset frame index
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,Stage,'MovAnalysis',[Animal '_CuedRewardedMov_FrameIndex.mat']),'CuedRewardedFrameIndex','CuedRewardedFrameIndex_z1','CuedRewardedFrameIndex_z2','-mat');
    for curr_field = 1:length(Imaging_Fields)
        Field = Fields{curr_field};
        Field_Dates = Imaging_Fields{curr_field}.Date;
        for curr_date = 1:length(Field_Dates)
            Date = Field_Dates{curr_date};
            df_f_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,Stage,'df_f',Field);
            date_index = ismember(Dates,Date);
            if Imaging_Fields{curr_field}.Piezo(curr_date) == 0
                if ismember(Animal,{'KP_3463808_1','KP_3461990_1','KP_3459921_1',}) && str2double(Date) > 180310
                    if curr_field == 1 || curr_field == 3
                        disp([Animal Date ': two Sessions two Fields, Field_' num2str(curr_field)]);
                        FrameOnset_index = CuedRewardedFrameIndex{date_index}(1);
                    elseif curr_field == 2 || curr_field == 4
                        disp([Animal Date ': two Sessions two Fields, Field_' num2str(curr_field)]);
                        FrameOnset_index = CuedRewardedFrameIndex{date_index}(2);
                    end
                else
                    FrameOnset_index = CuedRewardedFrameIndex{date_index};
                end
            elseif Imaging_Fields{curr_field}.Piezo(curr_date) == 1
                FrameOnset_index = CuedRewardedFrameIndex_z1{date_index};
            elseif Imaging_Fields{curr_field}.Piezo(curr_date) == 2
                FrameOnset_index = CuedRewardedFrameIndex_z2{date_index};
            end
            % load df/f
            if ~exist([df_f_Path filesep Animal '_' Date '_ROI_Traces.mat'])
                df_f_MovOnset_Alinged{curr_field,curr_date} = [];
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date} = [];
                continue
            end
            load([df_f_Path filesep Animal '_' Date '_ROI_Traces.mat'],'roi_trace_df_2','CaEvents_2','ZScore_2','truncatePoint','-mat');
            
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
                        df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}{1,curr_session}(:,trial) = temp_trace(FrameOnset_index{curr_session}(:,trial));
                    end
                    clear temp_trace temp_CaEvents temp_ZScore
                end
            end
            for roi = 1:size(roi_trace_df_2{curr_session},1)
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1} = cell2mat(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1});
                % Get rid of noise, threshold = 20;
                temp_index = abs(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1})>=20;
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}(temp_index) = nan;
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}(:,sum(temp_index)>=(size(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},1)/2)) = nan;
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);
            end
        end        
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,Stage,'df_f');
    save([TargetPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','df_f_MovOnset_Alinged_EachROI_MeanAXTrial','-v7.3');
end

%% Averaged activity after movement onset, combine trials across sessions first then average
clear all;
close all;
clc;

IN = 'VIP_SOM';
Animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R',...
        'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
Fields = {'Field_1','Field_2','Field_3','Field_4'};
Stage = 'Nai';

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,Stage,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','-mat');
    load([DataPath filesep Animal '_ImagingInfo.mat'],'-mat');
    for curr_field = 1:size(df_f_MovOnset_Alinged,1)
        nan_cells{curr_animal,1}{curr_field,1} = [];
        for curr_session = 1:size(df_f_MovOnset_Alinged(curr_field,:),2)
            if isempty(df_f_MovOnset_Alinged{curr_field,curr_session})
                Post_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session} = [];
                continue
            end
            if size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},1) == 36
                MovOnset_Frame = 8;
                Baseline_Frame = [3:5];
            elseif size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},1) == 72
                MovOnset_Frame = 15;
                Baseline_Frame = [7:10];
            end
            for curr_roi = 1:size(df_f_MovOnset_Alinged{curr_field,curr_session},1)
                temp = df_f_MovOnset_Alinged{curr_field,curr_session}{curr_roi};
                Post_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session}(curr_roi,:) = nanmean(temp(MovOnset_Frame:end,:)-...
                    repmat(nanmean(temp(Baseline_Frame,:),1),size(temp(MovOnset_Frame:end,:),1),1),1);
                Pre_MovOnset_Aligend_df_field{curr_animal}{curr_field,curr_session}(curr_roi,:) = nanmean(temp(1:MovOnset_Frame-1,:)-...
                    repmat(nanmean(temp(Baseline_Frame,:),1),size(temp(1:MovOnset_Frame-1,:),1),1),1);               
            end
        end

        % Get stage, in this case, only naive
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,1) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,1)),2);
        if ismember(curr_animal,[12:20,22])
            Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(RejectedROI.(Fields{curr_field}),1) = nan;
        end
        
        [nan_index,~] = find(isnan(Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}));
        nan_cells{curr_animal,1}{curr_field,1} = [nan_cells{curr_animal,1}{curr_field,1}, nan_index];
        nan_cells{curr_animal,1}{curr_field,1} = unique(nan_cells{curr_animal,1}{curr_field,1});
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(nan_cells{curr_animal,1}{curr_field,1},:) = [];
    end
    % Pool fields
    Post_MovOnset_Aligend_df_stage{curr_animal} = cell2mat(Post_MovOnset_Aligend_df_field_stage{curr_animal});
    % Average
    Post_MovOnset_Aligend_df_stage_Mean(curr_animal,:) = nanmean(Post_MovOnset_Aligend_df_stage{curr_animal});
    clear df_f_MovOnset_Alinged
end

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep Stage filesep 'Figures_SubPre'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
save([FigurePath filesep IN '_ActivityAnalysis_SubPre.mat'],'-v7.3');

[CX] = cbrewer('seq','RdPu',9);
color_value = nanmean(CX(6:7,:));    
stages = {'Naive','Early','Middle','Late'};
hM4Di_animals = {'CR_3974574-R','CR_3974573-O','CR_4017421-L','CR_4017421-R','CR_4042831-R','CR_4042831-LR','CR_4042830-O','CR_4042830-L','CR_4042832-O','CR_4042832-L','CR_4042832-R'};
mCherry_animals = {'CR_4017386-O','CR_4017386-L','CR_4017421-O','CR_4017421-LR','CR_4042831-O','CR_4042831-L','CR_4042830-R','CR_4113794-O','CR_4113794-L','CR_4153298-R','CR_4153298-LR'};
hM4Di_index = ismember(Animals,hM4Di_animals);
mCherry_index = ismember(Animals,mCherry_animals);

All_cell_stage_conc_h = cell2mat(Post_MovOnset_Aligend_df_stage(1,hM4Di_index)');
All_cell_stage_conc_m = cell2mat(Post_MovOnset_Aligend_df_stage(1,~hM4Di_index)');

figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = All_cell_stage_conc_m;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
bar([1],temp_1,'EdgeColor',[0.5 0.5 0.5],'FaceColor','none','linewidth',1,'barwidth',0.7);
line([1,1],[temp_1(1)-temp_2(1),temp_1(1)+temp_2(1)],'color',[0.5 0.5 0.5],'LineWidth',1);
temp_var = All_cell_stage_conc_h;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
bar([2],temp_1,'EdgeColor',color_value,'FaceColor','none','linewidth',1,'barwidth',0.7);
line([2,2],[temp_1(1)-temp_2(1),temp_1(1)+temp_2(1)],'color',color_value,'LineWidth',1);
xlim([0.4,2.6]); ylim([0,0.3]); xticks([1:2]); xticklabels({'mCherry','hM4Di'}); ylabel('Mean df/f');
axis square;
% Stats
for ii = 1:length(Animals)
    NeuroNum{ii,1} = length(cell2mat(Post_MovOnset_Aligend_df_field_stage{ii}));
    Animals_tag{ii,1} = repmat(ii,NeuroNum{ii,1},1);
    hM4Di_test{ii,1} = repmat(hM4Di_index(ii),NeuroNum{ii,1},1);
end
Drug = cell2mat(hM4Di_test);
Animals_test = cell2mat(Animals_tag);
Animals_test = Animals_test(:);
y = cell2mat(Post_MovOnset_Aligend_df_stage');
tbl = table(Drug,Animals_test,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Drug = nominal(tbl.Drug);
lme = fitlme(tbl,'y ~ 1 + Drug + (1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;    
line([0.9,2.1],[0.25 0.25],'color','k')
text(1.5,0.27,[num2str(pValue(2))],'horizontalalignment','center','fontsize',8);

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep Stage filesep 'Figures_SubPre'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.pdf']); pause(1);

% Distribution
min(cell2mat(Post_MovOnset_Aligend_df_stage'))
max(cell2mat(Post_MovOnset_Aligend_df_stage'))
edges = [-1.4:0.1:2];
figure; set(gcf,'position',[200,200,200,200]); hold on;
histogram(All_cell_stage_conc_m,edges,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.3,'Normalization','probability');
histogram(All_cell_stage_conc_h,edges,'edgecolor','none','facecolor',color_value,'facealpha',0.3,'Normalization','probability');
xlim([-1.4 2]); xlabel('df/f'); ylabel('Prob.')
ylim([0 0.36]); axis square;

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_distribution.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_distribution.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_distribution.pdf']); pause(1);

animal_combo_h = nchoosek([1:length(hM4Di_animals)],round(length(hM4Di_animals)/2));
animal_combo_m = nchoosek([1:length(hM4Di_animals)],round(length(hM4Di_animals)/2));

for ii_combo = 1:size(animal_combo_h,1)
    curr_combo_h = animal_combo_h(ii_combo,:);
    curr_combo_m = animal_combo_m(ii_combo,:);
    curr_h_index = ismember(Animals,hM4Di_animals(curr_combo_h));
    curr_m_index = ismember(Animals,mCherry_animals(curr_combo_m));
    temp_matrix_h = cell2mat(Post_MovOnset_Aligend_df_stage(1,curr_h_index)');
    temp_matrix_m = cell2mat(Post_MovOnset_Aligend_df_stage(1,curr_m_index)');
    combo_matrix(ii_combo,1) = nanmean(temp_matrix_h);
    combo_matrix(ii_combo,2) = nanmean(temp_matrix_m);   
end

edges = [0.09:0.01:0.32];
figure; set(gcf,'position',[200,200,200,200]); hold on;
histogram(combo_matrix(:,2),edges,'edgecolor','none','facecolor',[0.5 0.5 0.5],'facealpha',0.3,'Normalization','probability');
histogram(combo_matrix(:,1),edges,'edgecolor','none','facecolor',color_value,'facealpha',0.3,'Normalization','probability');
xlim([0.08 0.33]); xlabel('df/f'); ylabel('Prob.')
ylim([0 0.2]); axis square;
[h,p] = kstest2(combo_matrix(:,1),combo_matrix(:,2));
text(0.2,0.17,[num2str(p)],'horizontalalignment','center','fontsize',8);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_Mean_dfof_Post.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_Mean_dfof_Post.png']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_Mean_dfof_Post.pdf']); pause(1);

save([FigurePath filesep IN '_ActivityAnalysis_SubPre.mat'],'-v7.3');
close all

% Heatmap
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,Stage,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged_EachROI_MeanAXTrial','df_f_MovOnset_Alinged','-mat');
    clear temp_matrix
    for curr_field = 1:size(df_f_MovOnset_Alinged_EachROI_MeanAXTrial,1)
        curr_session = 1;
        if isempty(df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_session})
            temp_matrix{curr_field,curr_session} = [];
        else
            temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_session};
            temp(nan_cells{curr_animal,1}{curr_field,1},:) = [];
            trial_num(curr_session) = size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},2);
            if size(temp,2) == 72
                temp = temp(:,[1:2:72]);
            end
            temp = temp-repmat(nanmean(temp(:,3:5),2),1,36);
            temp_matrix{curr_field,curr_session} = temp; 
            clear temp
        end
    end
    clear df_f_MovOnset_Alinged_EachROI_MeanAXTrial df_f_MovOnset_Alinged
    for ii = 1
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{curr_animal,ii} = cell2mat(temp_matrix(:,ii));
    end
end

df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_h = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages(hM4Di_index,:);
df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_m = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages(mCherry_index,:);
figure; hold on; set(gcf,'color','w','position',[50 50 400 300]);
% warning off
% color_map = cbrewer('div','RdBu',64);
% warning on
% color_map = flipud(color_map);
subplot(1,2,1);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_m);
[~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
temp = temp(I,:);
% imagesc(temp,[-2 2]);colormap(color_map);
imagesc(temp,[-0.8 2]);colormap(mycmap);
line([8,8],ylim,'color','k','linestyle',':');
ylabel('# of cells');
xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
title('mCherry');
subplot(1,2,2);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_h);
[~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
temp = temp(I,:);
% imagesc(temp,[-2 2]);colormap(color_map);
imagesc(temp,[-0.8 2]);colormap(mycmap);
line([8,8],ylim,'color','k','linestyle',':');
ylabel('# of cells');
xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
title('hM4Di');

% Manually adjusted color map
% mycmap = colormap(gca);

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_conc.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_conc.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_conc.pdf']); pause(1);

figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_m);
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5 0.5 0.5],'linewidth',1);
temp = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages_h);
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value,'FaceAlpha',0.3);
plot(temp_1,'color',color_value,'linewidth',2);
ylim([-0.05 0.45]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_Temporal.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_Temporal.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post_Temporal.pdf']); pause(1);

save([FigurePath filesep IN '_ActivityAnalysis_SubPre.mat'],'-append');

