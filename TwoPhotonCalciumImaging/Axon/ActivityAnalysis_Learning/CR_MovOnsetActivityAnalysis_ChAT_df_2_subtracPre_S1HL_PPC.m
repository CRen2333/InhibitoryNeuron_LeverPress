%% 2P Data ChAT
%% Averaged activity after movement onset, combine trials across sessions first then average
clear all;
close all;
clc;

IN = 'ChAT';
Animals = {'CR_4383143-L','CR_4383143-R','CR_4429262-O','CR_4383144-R','CR_4412582-R','CR_4412583-LR'};
Fields = {'Field_1','Field_2'};

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
Animals = {'CR_4383143-L','CR_4383143-R','CR_4429262-O','CR_4383144-R','CR_4412582-R','CR_4412583-LR'};
Fields = {'Field_1','Field_2'};

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
            AxonNum{curr_animal,curr_field}(curr_session) = size(df_f_MovOnset_Alinged{curr_field,curr_session},1);
            TrialNum{curr_animal,curr_field}(curr_session) = size(df_f_MovOnset_Alinged{curr_field,curr_session}{1},2);
        end
        % Get 4 stages
        temp = cellfun(@(x) nanmean(x,2),Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,:),'UniformOutput',false);        
        temp_size = cellfun(@(x) size(x,2),Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,:),'UniformOutput',false);
        temp = temp';
        temp(cell2mat(temp_size)<3) = {[]};
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1} = temp(1);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,2} = temp(2);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,3} = temp(3);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,4} = temp(4:end);
        
        temp = cellfun(@(x) nanmean(x,2),Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,:),'UniformOutput',false);
        temp_size = cellfun(@(x) size(x,2),Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,:),'UniformOutput',false);
        temp = temp';
        temp(cell2mat(temp_size)<3) = {[]};
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,1} = temp(1);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,2} = temp(2);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,3} = temp(3);
        Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,4} = temp(4:end);
        
        % Average, pool across sessions then average
        for ii = 1:4
            Post_MovOnset_Aligend_df_field_stage_Mean_2{curr_field}(curr_animal,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,ii}));
            Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{curr_field}(curr_animal,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field_stage_sub{curr_animal}{curr_field,ii}));
        end
    end
    
    clear df_f_MovOnset_Alinged ZScore_MovOnset_Alinged
end

[CX] = cbrewer('div','RdBu',10);
color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];   
color_line = nanmean(color_value(2:3,:));
stages = {'Naive','Early','Middle','Late'};

curr_field = 1; % 1: S1HL, 2: PPC
figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{curr_field};
for curr_animal = 1:size(temp_var,1)
    plot(temp_var(curr_animal,:),'color',[0.5 0.5 0.5],'LineWidth',0.5);
end
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,2]); xticklabels(stages); ylabel('Mean df/f');
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
text(3,1.5,['p = ' num2str(pValue(2), '%.4f')],'color','k','fontsize',8,'HorizontalAlignment','center')
axis square
title('S1HL','fontsize',10);
% title('PPC','fontsize',10);

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep 'Figures_SubPre_OtherRegions'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '.pdf']); pause(1);

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
    
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,1} = cell2mat(temp_matrix(1,:));
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,2} = cell2mat(temp_matrix(2,:));
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,3} = cell2mat(temp_matrix(3,:));
        df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,4} = cell2mat(temp_matrix(4:end,:));
        clear temp_matrix
    end
end

% Average trace
clear temp_mean
curr_field = 2;
for curr_animal = 1:length(Animals)
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,1};
    temp_mean{1}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,2};
    temp_mean{2}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,3};
    temp_mean{3}(curr_animal,:) = nanmean(temp,1);
    temp = df_f_MovOnset_Alinged_EachROI_MeanAXTrial_field_stages{curr_animal}{curr_field,4};
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
xlim([1 76]); ylim([-0.1 0.65])
xticks([16,46,76]);xticklabels({'0','1','2'});
xlabel('Time (sec)');
ylabel('Mean df/f');
line([16 16], ylim, 'color','k','linestyle',':','linewidth',1);
axis square;

saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '_temporal.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '_temporal.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_' num2str(curr_field) '_temporal.pdf']); pause(1);

save([FigurePath filesep IN '_ActivityAnalysis_SubPre_OtherRegion.mat'],'-v7.3');

%% Compare between regions
Activity.S1HL = Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{1};
Activity.PPC = Post_MovOnset_Aligend_df_field_stage_sub_Mean_2{2};
temp = load(['Z:\People\Chi\TwoP_IN\ChAT\Figures_SubPre\ChAT_ActivityAnalysis_SubPre.mat'],'Post_MovOnset_Aligend_df_stage_sub_Mean_2');
Activity.M1 = temp.Post_MovOnset_Aligend_df_stage_sub_Mean_2;
clear temp

Regions = {'M1','S1HL','PPC'};
colors = [color_line;0.31,0.4,0.58;0.2,0.45,0.3];
figure; set(gcf,'pos',[200 200 450 200],'color','w'); hold on;
subplot(1,2,1); hold on;
for curr_region = 1:3
    region = Regions{curr_region};
    temp_var = Activity.(region);
    temp_1 = nanmean(temp_var);
    temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
    plot(temp_1,'color',colors(curr_region,:),'LineWidth',2);
    for ii = 1:4
        line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors(curr_region,:),'LineWidth',2);
    end
end
xlim([0.5,4.5]); ylim([0,1]); xticks([1:4]); xticklabels(stages); ylabel('Mean df/f');
axis square;

subplot(1,2,2); hold on;
for curr_region = 1:3
    region = Regions{curr_region};
    temp_var = Activity.(region);
    temp_baseline = nanmean(temp_var(:,1));
    temp_var = temp_var/temp_baseline;
    temp_1 = nanmean(temp_var);
    temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
    plot(temp_1,'color',colors(curr_region,:),'LineWidth',2);
    for ii = 1:4
        line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors(curr_region,:),'LineWidth',2);
    end
end
xlim([0.5,4.5]); ylim([0,1.5]); xticks([1:4]); xticklabels(stages); ylabel('norm. Mean df/f');
axis square;

saveas(gcf, [FigurePath filesep IN '_Mean_dfof_MultiRegions.fig']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_MultiRegions.png']); pause(1);
saveas(gcf, [FigurePath filesep IN '_Mean_dfof_MultiRegions.pdf']); pause(1);

temp_var = [Activity.M1;Activity.S1HL;Activity.PPC];
Sessions = repmat([1:4],size(temp_var,1),1);
Animals_test = repmat([1:size(temp_var,1)]',1,4);
Region_test = [repmat('M',size(Activity.M1,1),4);repmat('S',size(Activity.S1HL,1),4);repmat('P',size(Activity.PPC,1),4)];
Sessions = Sessions(:);
Animals_test = Animals_test(:);
Region_test = Region_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y,Region_test);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Region_test = nominal(tbl.Region_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + Region_test*Sessions + (1|Animals_test) + (Sessions-1|Animals_test)','DummyVarCoding','reference')

region_combo = [1,2;1,3;2,3];
Regions = {'M1','S1HL','PPC'};
for ii = 1:size(region_combo,1)
    field_1 = Regions{region_combo(ii,1)};
    field_2 = Regions{region_combo(ii,2)};
    temp_var = [Activity.(field_1);Activity.(field_2)];
    Sessions = repmat([1:4],size(temp_var,1),1);
    Animals_test = repmat([1:size(temp_var,1)]',1,4);
    Region_test = [repmat(1,size(Activity.(field_1),1),4);repmat(2,size(Activity.(field_2),1),4)];
    Sessions = Sessions(:);
    Animals_test = Animals_test(:);
    Region_test = Region_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y,Region_test);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Region_test = nominal(tbl.Region_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + Region_test*Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue_regions(ii) = stats.pValue(4);
    beta_regions(ii) = beta(4);
end

