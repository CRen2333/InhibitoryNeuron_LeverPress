%% 2P Data VIP, mAChR and nAChR antagonist seperately
clear all;
close all;
clc;

IN = 'VIP';
SessionType = 'Nai_mAnt';
switch SessionType
    case 'Nai_nAnt'
        Animals = {'CR_4259757-R','CR_4302983-LL','CR_4302983-R','CR_4302984-O','CR_4302984-R','CR_4302984-LR'};
        Drug = 'nAnt';
    case 'Nai_mAnt'
        Animals = {'CR_4303048-L','CR_4303048-R','CR_4303048-LR','CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'};
        Drug = 'mAnt';
end
Fields = {'Field_1','Field_2'};

for curr_animal = 1:length(Animals)
    clearvars -except Animals Drug SessionType IN curr_animal Fields
    Animal = Animals{curr_animal};
    LeverTrace_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'Pharm',SessionType,'LeverTrace');
    Dates = dir(LeverTrace_Path);
    Dates = {Dates.name};
    Dates = sort(Dates(cellfun(@(x) contains(x, '21')||contains(x, '18')||contains(x, '19'), Dates)))';
    % load imaging info
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'Pharm',SessionType,'df_f',[Animal '_ImagingInfo.mat']),'Imaging_Fields','-mat');
    % load mov-onset frame index
    load(fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'Pharm',SessionType,'MovAnalysis',[Animal '_CuedRewardedMov_FrameIndex.mat']),'CuedRewardedFrameIndex','CuedRewardedFrameIndex_z1','CuedRewardedFrameIndex_z2','-mat');
    for curr_field = 1:length(Imaging_Fields)
        Field = Fields{curr_field};
        Field_Dates = Imaging_Fields{curr_field}.Date;
        for curr_date = 1:length(Field_Dates)
            Date = Field_Dates{curr_date};
            df_f_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'Pharm',SessionType,'df_f',Field);
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
                CaEvents_MovOnset_Alinged{curr_field,curr_date} = [];
                ZScore_MovOnset_Alinged{curr_field,curr_date} = [];
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date} = [];
                CaEvents_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date} = [];
                ZScore_MovOnset_Alinged_EachROI_MeanAXTrial = [];
                continue
            end
            load([df_f_Path filesep Animal '_' Date '_ROI_Traces.mat'],'roi_trace_df_2','CaEvents_2','ZScore_2','truncatePoint','-mat');

            for curr_session = 1:length(FrameOnset_index)
                disp([Animal ' ' Field '_' Date ' Session' num2str(curr_session)]);
                end_control = FrameOnset_index{curr_session}(end,:) <= round(truncatePoint{curr_session}(2));
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session}(:,end_control);
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session} - round(truncatePoint{curr_session}(1))+1;
                start_control = FrameOnset_index{curr_session}(1,:) > 0;
                FrameOnset_index{curr_session} = FrameOnset_index{curr_session}(:,start_control);
                for roi = 1:size(roi_trace_df_2{curr_session},1)
                    temp_trace = roi_trace_df_2{curr_session}(roi,:);
                    temp_CaEvents = CaEvents_2{curr_session}(roi,:);
                    temp_ZScore = ZScore_2{curr_session}(roi,:);
                    for trial = 1:size(FrameOnset_index{curr_session},2)
                        df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}{1,curr_session}(:,trial) = temp_trace(FrameOnset_index{curr_session}(:,trial));
                        CaEvents_MovOnset_Alinged{curr_field,curr_date}{roi,1}{1,curr_session}(:,trial) = temp_CaEvents(FrameOnset_index{curr_session}(:,trial));
                        ZScore_MovOnset_Alinged{curr_field,curr_date}{roi,1}{1,curr_session}(:,trial) = temp_ZScore(FrameOnset_index{curr_session}(:,trial));
                    end
                    clear temp_trace temp_CaEvents temp_ZScore
                end
            end
            for roi = 1:size(roi_trace_df_2{curr_session},1)
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1} = cell2mat(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1});
                CaEvents_MovOnset_Alinged{curr_field,curr_date}{roi,1} = cell2mat(CaEvents_MovOnset_Alinged{curr_field,curr_date}{roi,1});
                ZScore_MovOnset_Alinged{curr_field,curr_date}{roi,1} = cell2mat(ZScore_MovOnset_Alinged{curr_field,curr_date}{roi,1});
                % Get rid of noise, threshold = 20;
                temp_index = abs(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1})>=20;
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}(temp_index) = nan;
                df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1}(:,sum(temp_index)>=(size(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},1)/2)) = nan;
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);
                CaEvents_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(CaEvents_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);
                ZScore_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(ZScore_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);                
            end
        end        
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'Pharm',SessionType,'df_f');
    SessionTage = Imaging_Fields{1}.SessionType;
    save([TargetPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','df_f_MovOnset_Alinged_EachROI_MeanAXTrial',...
        'CaEvents_MovOnset_Alinged','CaEvents_MovOnset_Alinged_EachROI_MeanAXTrial',...
        'ZScore_MovOnset_Alinged','ZScore_MovOnset_Alinged_EachROI_MeanAXTrial','-v7.3');
end

%% Averaged activity after movement onset, combine trials across sessions first then average
clear all;
close all;
clc;

IN = 'VIP';
% SessionType = 'Nai_nAnt';
SessionType = 'Nai_mAnt';
switch SessionType
    case 'Nai_nAnt'
        Animals = {'CR_4259757-R','CR_4302983-LL','CR_4302983-R','CR_4302984-O','CR_4302984-R','CR_4302984-LR'};
        Drug = 'nAnt';
        Post_MovOnset_Aligend_df_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_Ca_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_ZS_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_df_sub_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_Ca_sub_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_ZS_sub_Pharm_Mean = nan(length(Animals),2);
    case 'Nai_mAnt'
        Animals = {'CR_4303048-L','CR_4303048-R','CR_4303048-LR','CR_4365784-O','CR_4365784-L','CR_4365807-O','CR_4365807-L'};
        Drug = 'mAnt';
        Post_MovOnset_Aligend_df_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_Ca_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_ZS_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_df_sub_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_Ca_sub_Pharm_Mean = nan(length(Animals),2);
        Post_MovOnset_Aligend_ZS_sub_Pharm_Mean = nan(length(Animals),2);
end
Fields = {'Field_1','Field_2'};

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    SessionType_2 = SessionType;
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'Pharm',SessionType_2,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','-mat');
    load([DataPath filesep Animal '_ImagingInfo.mat'],'Imaging_Fields','-mat');    
    for curr_field = 1:size(df_f_MovOnset_Alinged,1)
        nan_cells{curr_animal,1}{curr_field,1} = [];
        for curr_session = 1:min(5,size(df_f_MovOnset_Alinged(curr_field,:),2))
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
        end

        if strcmp(SessionType,'Nai_nAnt') && strcmp(Imaging_Fields{1}.SessionType{1},'Nai_nAnt')
            disp([Animal 'SWARP']);
            temp = Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,:);
            Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,1) = temp(:,2);
            Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,2) = temp(:,1);
            
            temp = Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,:);
            Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,1) = temp(:,2);
            Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,2) = temp(:,1);
        end            
        
        if strcmp(SessionType,'Nai_mAnt') && strcmp(Imaging_Fields{1}.SessionType{1},'Nai_mAnt')
            disp([Animal 'SWARP']);
            temp = Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,:);
            Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,1) = temp(:,2);
            Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,2) = temp(:,1);
            
            temp = Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,:);
            Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,1) = temp(:,2);
            Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,2) = temp(:,1);
        end      
        
        % Get mean then combine fields
        for ii = 1:length(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,:))
            Post_MovOnset_Aligend_df_field_sub_Pharm{curr_animal}{curr_field,1}(:,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field_sub{curr_animal}(curr_field,ii)),2);            
            Post_MovOnset_Aligend_df_field_Pharm{curr_animal}{curr_field,1}(:,ii) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,ii)),2);
        end
        
                        
        [nan_index,~] = find(isnan(sum(Post_MovOnset_Aligend_df_field_Pharm{curr_animal}{curr_field,1},2)));
        nan_cells{curr_animal,1}{curr_field,1} = [nan_cells{curr_animal,1}{curr_field,1}, nan_index];
        nan_cells{curr_animal,1}{curr_field,1} = unique(nan_cells{curr_animal,1}{curr_field,1});
        
        Post_MovOnset_Aligend_df_field_sub_Pharm{curr_animal}{curr_field,1}(nan_cells{curr_animal,1}{curr_field,1},:) = [];
        Post_MovOnset_Aligend_df_field_Pharm{curr_animal}{curr_field,1}(nan_cells{curr_animal,1}{curr_field,1},:) = [];
        
        NeuronNum{curr_animal}(curr_field) = size(Post_MovOnset_Aligend_df_field_Pharm{curr_animal}{curr_field,1},1);        
    end
    
    % Pool fields
    Post_MovOnset_Aligend_df_sub_Pharm{curr_animal} = cell2mat(Post_MovOnset_Aligend_df_field_sub_Pharm{curr_animal});
    Post_MovOnset_Aligend_df_Pharm{curr_animal} = cell2mat(Post_MovOnset_Aligend_df_field_Pharm{curr_animal});
    % Average
    Post_MovOnset_Aligend_df_sub_Pharm_Mean(curr_animal,1:size(Post_MovOnset_Aligend_df_Pharm{curr_animal},2)) = nanmean(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal});
    Post_MovOnset_Aligend_df_Pharm_Mean(curr_animal,1:size(Post_MovOnset_Aligend_df_Pharm{curr_animal},2)) = nanmean(Post_MovOnset_Aligend_df_Pharm{curr_animal});
    clear df_f_MovOnset_Alinged CaEvents_MovOnset_Alinged ZScore_MovOnset_Alinged SessionTage Imaging_Fields
end

%%
clear all
close all
clc

IN = 'VIP';
SessionType = 'Nai_mAnt';
Drug = 'mAnt';

load(['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType '\Figures_SubPre\VIP_ActivityAnalysis_SubPre.mat'],'Animals','Post_MovOnset_Aligend_df_sub_Pharm','-mat');
% Post mvmOnset activity avreage across neurons
Post_MovOnset_Aligend_df_sub_Pharm = Post_MovOnset_Aligend_df_sub_Pharm;
for curr_animal  = 1:length(Post_MovOnset_Aligend_df_sub_Pharm)
    neuronum(curr_animal) = size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal},1);
end
Animal_tag = [];
for ii = 1:length(neuronum)
    Animal_tag = [Animal_tag;ones(neuronum(ii),1)*ii];
end
for ii = 1:length(Post_MovOnset_Aligend_df_sub_Pharm)
    Post_MovOnset_Aligend_df_sub_Pharm{ii} = Post_MovOnset_Aligend_df_sub_Pharm{ii}(:,1:2);
end
All_cell_Pharm_conc = cell2mat(Post_MovOnset_Aligend_df_sub_Pharm');

[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
color_line = mean(color_value(2:3,:));
Pharms = {'Ctrl',Drug};

figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
temp_var = All_cell_Pharm_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
area([1.5,2.5],[1,1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
area([1.5,2.5],[-1,-1],'facecolor',[0.95,0.95,0.95],'edgecolor','none');
plot(temp_1,'color',color_value(1,:),'LineWidth',2);
for ii = 1:2
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_value(1,:),'LineWidth',2);
end
xlim([0.5,2.5]); xticks([1:2]); xticklabels(Pharms); ylabel('Mean df/f'); axis square
ylim([0,0.35]);
% yticks([0:0.1:0.3]);
% random effect modal
Task = repmat([0:1],size(temp_var,1),1);
Task = Task(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
Neurons = repmat([1:size(temp_var,1)],1,2);
Neurons = Neurons(:);
y = temp_var(:);
tbl = table(Animals_test,Task,y,Neurons);
tbl.Animals_test = nominal(tbl.Animals_test);
% tbl.Neurons = nominal(tbl.Neurons);
tbl.Task = nominal(tbl.Task);
% lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)');
lme = fitlme(tbl,'y ~ 1 + Task + (1|Animals_test)+(Task-1|Animals_test)');

[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue
text(1.5,0.3,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType filesep 'Figures_SubPre'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.pdf']); pause(1);

% Distribution
figure; set(gcf,'color','w','position',[200,200,200,200]); hold on;
for ii = 1:2
    Delta_df_against_C(:,ii) = All_cell_Pharm_conc(:,ii)-All_cell_Pharm_conc(:,1);
end
edges = [-2:0.1:2];
temp_n = histcounts(Delta_df_against_C(:,2),edges,'Normalization','probability');
% stairs(edges(1:end-1)+0.05,temp_n,'color',color_value(1,:),'linewidth',1);
histogram(Delta_df_against_C(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(1,:));
ylim([0 0.2]); xlim([edges(1) edges(end)]);
line([0 0],ylim,'color','k','linestyle',':'); yticks([0:0.1:0.2]); yticklabels({'0','10','20'});
xlabel('delta df/f'); ylabel('Fraction'); axis square;
temp_var = All_cell_Pharm_conc(:,[1,ii]);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_sub_Pharm{curr_animal},1),2)*curr_animal];
end
Task_test = repmat([0,1],size(temp_var,1),1);
y = temp_var(:);
Animals_test = Animals_test(:);
Task_test = Task_test(:);
tbl = table(Animals_test,Task_test,y);
tbl.Animals_test = nominal(tbl.Animals_test);
tbl.Task_test = nominal(tbl.Task_test);
lme = fitlme(tbl,'y ~ 1 + Task_test + (1|Animals_test) + (Task_test-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue
text(0,0.25,['p =' num2str(pValue(2),'%.4f')],'color','k','HorizontalAlignment','center')
axis square;

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.pdf']); pause(1);

% Temporal traces
load(['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType '\Figures_SubPre\' IN '_ActivityAnalysis_SubPre.mat'],'df_f_MovOnset_Alinged_EachROI_MeanAXTrial_sub_Pharm','-mat');
for ii = 1:2
    all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_sub_Pharm(:,ii));
end
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
warning off
color_map = cbrewer('div','RdBu',64);
warning on
color_map = flipud(color_map);
for ii = 1:2
    subplot(1,2,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
%     imagesc(temp,[-2 2]);
%     colormap(color_map);
    imagesc(temp,[-1 2]);
    colormap(mycmap);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end
% mycmap = colormap(gca);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_sort.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_sort.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_sort.pdf']); pause(1);

% sort based on naive
[~,I_ctrl] = sort(All_cell_Pharm_conc(:,1),'descend');
figure; hold on; set(gcf,'color','w','position',[50 50 400 500]);
for ii = 1:2
    subplot(1,2,ii);
    imagesc(all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii}(I_ctrl,:),[-0.5 2.2]);
    line([8,8],ylim,'color','w','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Pharms{ii});
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.pdf']); pause(1);

for ii = 1:length(Animals)
    temp_tag(:,ii) = Animal_tag(I_ctrl)==ii;
end
window_size = 50;
for ii = 1:length(Animals)
    Animal_tag_prob_sort(:,ii) = movsum(temp_tag(:,ii),window_size)/sum(temp_tag(:,ii));
%     Animal_tag_prob_sort(:,ii) = movsum(temp_tag(:,ii),window_size)/window_size;
end
figure; hold on; set(gcf,'color','w','position',[50 50 100 500]);
Animals_color = colormap;
Animals_color = Animals_color(4:6:64,:);
for ii = 1:length(Animals)
    plot(Animal_tag_prob_sort(:,ii),[1:length(Animal_tag_prob_sort(:,ii))],'color',Animals_color(ii,:));
end
xlim([0 0.25]);ylim([1 length(Animal_tag)]);set(gca,'ydir','reverse');
xlabel('Prob.');title('Animal ID');
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_AnimalID.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_AnimalID.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc_AnimalID.pdf']); pause(1);

figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 1;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
% ylim([-0.03 0.3]);
ylim([-0.05 0.45]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square;
for ii = 1:36
    temp_data = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{2}(:,ii)-all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_Pharm{1}(:,ii);
    [bootstat,~] = bootstrp(1000,@mean,temp_data);
    if prctile(bootstat,2.5)>0 || prctile(bootstat,97.5)<0
%         line([ii-0.5,ii+0.5],[0.36 0.36],'color','k','linewidth',1);
        line([ii-0.5,ii+0.5],[0.19 0.19],'color','k','linewidth',1);
    end
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.pdf']); pause(1);

% significantly changed fraction
load(['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType '\Figures_SubPre\' IN '_ActivityAnalysis_SubPre.mat'], 'Delta_df_L_N_CellTag')
All_Delta_df_against_C_CellTag = [];
for ii = 1:length(Animals)
    All_Delta_df_against_C_CellTag = [All_Delta_df_against_C_CellTag;cell2mat(Delta_df_L_N_CellTag{ii})];
end
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 2;
pie([sum(All_Delta_df_against_C_CellTag(:,ii)==1),sum(All_Delta_df_against_C_CellTag(:,ii)==-1),sum(All_Delta_df_against_C_CellTag(:,ii)==0)]);
axis square; axis off;

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.pdf']); pause(1);

animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_matrix = [];
    for jj = 1:length(curr_combo)
        temp_matrix = [temp_matrix;cell2mat(Delta_df_L_N_CellTag{curr_combo(jj)})];
    end
    for ii = 2
        Fraction_CellTag_combo{ii}(ii_combo,1) = sum(temp_matrix(:,ii)==1)/length(temp_matrix(:,ii));
        Fraction_CellTag_combo{ii}(ii_combo,2) = sum(temp_matrix(:,ii)==-1)/length(temp_matrix(:,ii));
        Fraction_CellTag_combo{ii}(ii_combo,3) = sum(temp_matrix(:,ii)==0)/length(temp_matrix(:,ii)); 
    end
end
figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
plot(sum(All_Delta_df_against_C_CellTag(:,2)==1)/length(All_Delta_df_against_C_CellTag)*100,sum(All_Delta_df_against_C_CellTag(:,2)==-1)/length(All_Delta_df_against_C_CellTag)*100,'marker','x','color','k');
xlim([0 40]);ylim([0 40]);
plot([0,40],[0,40],'linestyle',':','color','k');
xlabel('Inc. fraction');ylabel('Dec. fraction'); axis square;
xticks([0:10:40]);% xticklabels({'inc.','dec.','no'});
yticks([0:10:40]);    
axis square;

saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.png']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.pdf']); pause(1);        

figure; hold on; set(gcf,'color','w','position',[50 50 200 200]);
ii = 2;
temp_var = Fraction_CellTag_combo{ii};
temp_var_1 = Fraction_CellTag_combo{ii}(:,1)-Fraction_CellTag_combo{ii}(:,2);
shake_x = (rand(size(temp_var_1))-0.5)*0.3;
plot(shake_x,temp_var_1,'color',[0.5 0.5 0.5],'Linestyle','none','marker','.','markersize',8);
xlim([-0.5 0.5]);xticks([0]);xticklabels({'inc.-dec.'});
ylim([-0.35 0]);
% ylim([-0.15 0.05]);
line(xlim,[0 0],'color','k','linestyle',':');
yticks([-0.15:0.05:0.05]);yticklabels({'-15','-10','-5','0','5'});
ylabel('Fraction');
title(Pharms{ii});
axis square;

saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.png']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.pdf']); pause(1);        

%% Fraction of MVM modulated neurons
load(['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType '\Figures_SubPre\' IN '_ActivityAnalysis_SubPre.mat'],'movement_activated_cells_index','movement_suppressed_cells_index','-mat');

All_movement_activated_cells_index = cell2mat(movement_activated_cells_index');
All_movement_suppressed_cells_index = cell2mat(movement_suppressed_cells_index');
figure; hold on; set(gcf,'color','w','position',[200 200 400 200]);
for ii = 1:2
    subplot(1,2,ii);
    temp = All_movement_activated_cells_index(:,ii)+All_movement_suppressed_cells_index(:,ii);
    pie([sum(temp==1),sum(temp==-1),sum(temp==0)]);
%     clear temp;
    axis square;
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.pdf']); pause(1);

All_movement_modulated_cells_index_stages = All_movement_activated_cells_index+All_movement_suppressed_cells_index;
for ii = 1:2
    temp_bar(ii,:) = [sum(All_movement_modulated_cells_index_stages(:,ii)==1),sum(All_movement_modulated_cells_index_stages(:,ii)==-1),...
        sum(All_movement_modulated_cells_index_stages(:,ii)==0)];
end

figure; hold on; set(gcf,'color','w','position',[400 400 140 200]);
bar(temp_bar(1:2,:)/sum(temp_bar(1,:))*100,'stacked','barwidth',0.7);
xlim([0.5,2.5]);ylim([0 100]); % [0.75,0,0;0,0.34,0.42;0.8,0.8,0.8]
xticks([1:2]);xticklabels({'Ctrl',Drug});ylabel('Fraction');

saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_stackbar.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_stackbar.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_stackbar.pdf']); pause(1);


clear pval
x1 = [repmat('N',sum(temp_bar(1,:)),1);repmat('E',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,1);All_movement_modulated_cells_index_stages(:,2)];

% compare between stages
x1 = [repmat('C',size(All_movement_activated_cells_index,1),1);
    repmat('D',size(All_movement_activated_cells_index,1),1);];
x2 = All_movement_activated_cells_index+All_movement_suppressed_cells_index;
x2 = x2(:);
[tbl,chi2stat,pval] = crosstab(x1,x2);

animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_matrix = [];
    for jj = 1:length(curr_combo)
        temp_matrix = [temp_matrix;movement_activated_cells_index{curr_combo(jj)}+movement_suppressed_cells_index{curr_combo(jj)}];
    end
    for ii = 1:2
        Fraction_MvmMod_combo{ii}(ii_combo,1) = sum(temp_matrix(:,ii)==1)/length(temp_matrix(:,ii));
        Fraction_MvmMod_combo{ii}(ii_combo,2) = sum(temp_matrix(:,ii)==-1)/length(temp_matrix(:,ii));
        Fraction_MvmMod_combo{ii}(ii_combo,3) = sum(temp_matrix(:,ii)==0)/length(temp_matrix(:,ii)); 
    end
end
figure; hold on; set(gcf,'color','w','position',[50 50 400 200]);
for ii = 1:2
    subplot(1,2,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii};
    plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
    plot(temp_bar(ii,1)/length(All_movement_modulated_cells_index_stages)*100,temp_bar(ii,2)/length(All_movement_modulated_cells_index_stages)*100,'marker','x','color','k');
    xlim([0 60]);ylim([0 60]);
    plot([0,60],[0,60],'linestyle',':','color','k');
    xlabel('Act. fraction');ylabel('Sup. fraction'); axis square;
    xticks([0:20:60]);
    yticks([0:20:60]);
%     plot(temp_var','color',[0.5 0.5 0.5]);
%     temp_1 = nanmean(temp_var);
%     temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
%     plot(temp_1,'color','k','LineWidth',1);
%     for jj = 1:3
%         line([jj,jj],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color','k','LineWidth',1);
%     end
%     xlim([0.7 3.3]);xticks([1:3]);xticklabels({'act.','sup.','no'});
%     ylabel('Fraction')
    title(Pharms{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.png']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.pdf']); pause(1);        

figure; hold on; set(gcf,'color','w','position',[50 50 400 200]);
for ii = 1:2
    subplot(1,2,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii};
    temp_var_1 = Fraction_MvmMod_combo{ii}(:,1)-Fraction_MvmMod_combo{ii}(:,2);
    shake_x = (rand(size(temp_var_1))-0.5)*0.3;
    plot(shake_x,temp_var_1,'color',[0.5 0.5 0.5],'Linestyle','none','marker','.','markersize',8);
    xlim([-0.5 0.5]);xticks([0]);xticklabels({'act.-sup.'});
    line(xlim,[0 0],'color','k','linestyle',':');
    ylabel('Fraction');
    title(Pharms{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.pdf'],'-dpdf','-bestfit'); pause(1);        

figure; hold on; set(gcf,'color','w','position',[50 50 400 200]);
for ii = 2
    subplot(1,2,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii}-Fraction_MvmMod_combo{1};
    temp_var_1 = temp_var(:,1);
    temp_var_2 = temp_var(:,2);
    shake_x = (rand(size(temp_var_1))-0.5)*0.3;
    plot(shake_x-0.5,temp_var_1,'color',[191,0,0]/255,'Linestyle','none','marker','.','markersize',8);
    plot(shake_x+0.5,temp_var_2,'color',[0,87,107]/255,'Linestyle','none','marker','.','markersize',8);
    xlim([-1 1]);xticks([-0.5 0.5]);xticklabels({'act.','sup.'});
    line(xlim,[0 0],'color','k','linestyle',':');
    ylabel('Fraction');
    title(Pharms{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.pdf'],'-dpdf','-bestfit'); pause(1);        


% MVM modulated neurons activity onset
load(['Z:\People\Chi\TwoP_IN\' IN filesep 'Pharm' filesep SessionType '\Figures_SubPre\' IN '_ActivityAnalysis_SubPre.mat'],'movement_activated_cells_onset','movement_suppressed_cells_onset','-mat');
for ii = 1:2
    temp = cell2mat(movement_activated_cells_onset');
    temp = temp(:,ii);
    temp = temp(~isnan(temp));
    All_movement_activated_cells_onset{ii} = temp;
    clear temp;
    temp = cell2mat(movement_suppressed_cells_onset');
    temp = temp(:,ii);
    temp = temp(~isnan(temp));
    All_movement_suppressed_cells_onset{ii} = temp;
    clear temp;
end

figure; hold on; set(gcf,'color','w','position',[50 50 200 400]);
subplot(2,1,1); hold on;
ii = 1;
cellnum = length(All_movement_activated_cells_onset{ii});
plot(sort(All_movement_activated_cells_onset{ii},'ascend'),[1:cellnum]/cellnum,'color',[0.5,0.5,0.5],'linewidth',1);
ii = 2;
cellnum = length(All_movement_activated_cells_onset{ii});
plot(sort(All_movement_activated_cells_onset{ii},'ascend'),[1:cellnum]/cellnum,'color',color_value(1,:),'linewidth',1);
line([0 0],ylim,'color','k','linestyle',':');
xlim([-0.5 2]);
xlabel('Activity onset (s)'); ylabel('cum. prob.');
title('Activated'); axis square;
subplot(2,1,2); hold on;
ii = 1;
cellnum = length(All_movement_suppressed_cells_onset{ii});
plot(sort(All_movement_suppressed_cells_onset{ii},'ascend'),[1:cellnum]/cellnum,'color',[0.5,0.5,0.5],'linewidth',1);
ii = 2;
cellnum = length(All_movement_suppressed_cells_onset{ii});
plot(sort(All_movement_suppressed_cells_onset{ii},'ascend'),[1:cellnum]/cellnum,'color',color_value(1,:),'linewidth',1);
line([0 0],ylim,'color','k','linestyle',':');
xlim([-0.5 2]);
xlabel('Activity onset (s)'); ylabel('cum. prob.');
title('Suppressed'); axis square;
   
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.png']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.pdf']); pause(1);        

% Check trace
for ii = 1:2
    movement_activated_cells_trace_pool{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_sub_Pharm(:,ii));
    movement_activated_cells_trace_pool{ii} = movement_activated_cells_trace_pool{ii}(logical(All_movement_activated_cells_index(:,ii)),:);
end
for ii = 1:2
    movement_suppressed_cells_trace_pool{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_sub_Pharm(:,ii));
    movement_suppressed_cells_trace_pool{ii} = movement_suppressed_cells_trace_pool{ii}(logical(All_movement_suppressed_cells_index(:,ii)),:);
end
% across neurons
figure; hold on; set(gcf,'color','w','position',[200 200 200 400]);
subplot(2,1,1); hold on;
ii = 1;
temp = movement_activated_cells_trace_pool{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = movement_activated_cells_trace_pool{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1); ylim([-0.05, 0.8]);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square; title('Activated');
subplot(2,1,2); hold on;
ii = 1;
temp = movement_suppressed_cells_trace_pool{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.3);
plot(temp_1,'color',[0.5,0.5,0.5],'linewidth',2);
ii = 2;
temp = movement_suppressed_cells_trace_pool{ii};
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',color_value(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',color_value(1,:),'linewidth',2);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1); ylim([-0.7, 0.1]);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square; title('Suppressed');
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.pdf']); pause(1);


%% Check relationship between significalty changed activity and mvm modulated

% mvm in celltag
temp_1 = All_Delta_df_against_C_CellTag;
temp_2 = All_movement_activated_cells_index;
temp_3 = All_movement_suppressed_cells_index;
temp_4 = ones(size(temp_1));
temp_4 = temp_4-logical(temp_2+temp_3);
for ii = 2
    % celltage inc:
    % act to act, act to sup, act to no
    index_celltag = temp_1(:,ii)==1;
    index_naive = logical(temp_2(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_inc(ii,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % sup to act, sup to sup, sup to no
    index_naive = logical(temp_3(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_inc(ii,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % no to act, no to sup, no to no
    index_naive = logical(temp_4(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_inc(ii,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_inc(ii,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

    Fraction_mvm_in_celltage_inc(isnan(Fraction_mvm_in_celltage_inc))=0;

    % celltage dec:
    % act to act, act to sup, act to no
    index_celltag = temp_1(:,ii)==-1;
    index_naive = logical(temp_2(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_dec(ii,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % sup to act, sup to sup, sup to no
    index_naive = logical(temp_3(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_dec(ii,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % no to act, no to sup, no to no
    index_naive = logical(temp_4(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_dec(ii,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_dec(ii,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

    Fraction_mvm_in_celltage_dec(isnan(Fraction_mvm_in_celltage_inc))=0;

    % celltage no:
    % act to act, act to sup, act to no
    index_celltag = temp_1(:,ii)==0;
    index_naive = logical(temp_2(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_no(ii,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % sup to act, sup to sup, sup to no
    index_naive = logical(temp_3(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_no(ii,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
    % no to act, no to sup, no to no
    index_naive = logical(temp_4(:,1));
    index_act = logical(temp_2(:,ii));
    index_sup = logical(temp_3(:,ii));
    index_no = logical(temp_4(:,ii));
    Fraction_mvm_in_celltage_no(ii,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
    Fraction_mvm_in_celltage_no(ii,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

    Fraction_mvm_in_celltage_no(isnan(Fraction_mvm_in_celltage_inc))=0;
end       

figure; set(gcf,'color','w','position',[2000 100 600 200]); hold on;
ii = 2;
subplot(1,3,1); hold on
temp = Fraction_mvm_in_celltage_inc(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');
subplot(1,3,2); hold on
temp = Fraction_mvm_in_celltage_dec(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');
subplot(1,3,3); hold on
temp = Fraction_mvm_in_celltage_no(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');

saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.pdf']); pause(1);   

animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_1 = [];
    for jj = 1:length(curr_combo)
        temp_1 = [temp_1;cell2mat(Delta_df_L_N_CellTag{curr_combo(jj)})];
    end
    temp_2 = [];
    for jj = 1:length(curr_combo)
        temp_2 = [temp_2;movement_activated_cells_index{curr_combo(jj)}];
    end
    temp_3 = [];
    for jj = 1:length(curr_combo)
        temp_3 = [temp_3;movement_suppressed_cells_index{curr_combo(jj)}];
    end
    temp_4 = ones(size(temp_1));
    temp_4 = temp_4-logical(temp_2+temp_3);
    for ii = 2
        % celltage inc:
        % act to act, act to sup, act to no
        index_celltag = temp_1(:,ii)==1;
        index_naive = logical(temp_2(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % sup to act, sup to sup, sup to no
        index_naive = logical(temp_3(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % no to act, no to sup, no to no
        index_naive = logical(temp_4(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_inc{ii}(ii_combo,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

        combo_Fraction_mvm_in_celltage_inc{ii}(isnan(combo_Fraction_mvm_in_celltage_inc{ii}))=0;

        % celltage dec:
        % act to act, act to sup, act to no
        index_celltag = temp_1(:,ii)==-1;
        index_naive = logical(temp_2(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % sup to act, sup to sup, sup to no
        index_naive = logical(temp_3(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % no to act, no to sup, no to no
        index_naive = logical(temp_4(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_dec{ii}(ii_combo,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

        combo_Fraction_mvm_in_celltage_dec{ii}(isnan(combo_Fraction_mvm_in_celltage_inc{ii}))=0;

        % celltage no:
        % act to act, act to sup, act to no
        index_celltag = temp_1(:,ii)==0;
        index_naive = logical(temp_2(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,1) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,2) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,3) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % sup to act, sup to sup, sup to no
        index_naive = logical(temp_3(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,4) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,5) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,6) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);
        % no to act, no to sup, no to no
        index_naive = logical(temp_4(:,1));
        index_act = logical(temp_2(:,ii));
        index_sup = logical(temp_3(:,ii));
        index_no = logical(temp_4(:,ii));
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,7) = sum(index_naive.*index_act.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,8) = sum(index_naive.*index_sup.*index_celltag)/sum(index_celltag);
        combo_Fraction_mvm_in_celltage_no{ii}(ii_combo,9) = sum(index_naive.*index_no.*index_celltag)/sum(index_celltag);

        combo_Fraction_mvm_in_celltage_no{ii}(isnan(combo_Fraction_mvm_in_celltage_inc{ii}))=0;
    end       
end
figure; set(gcf,'color','w','position',[2000 100 600 200]); hold on
ii = 2;
subplot(1,3,1); hold on
temp = combo_Fraction_mvm_in_celltage_inc{ii};
plot(temp','color',[0.8,0.8,0.8]);
bar([1:9],nanmean(temp),'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');
subplot(1,3,2); hold on
temp = combo_Fraction_mvm_in_celltage_dec{ii};
plot(temp','color',[0.8,0.8,0.8]);
bar([1:9],nanmean(temp),'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');
subplot(1,3,3); hold on
temp = combo_Fraction_mvm_in_celltage_no{ii};
plot(temp','color',[0.8,0.8,0.8]);
bar([1:9],nanmean(temp),'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');

saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.png']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.pdf']); pause(1);   

% decrease only    
figure; set(gcf,'color','w','position',[2000 100 600 200]); hold on
ii = 2;
subplot(1,3,2); hold on
temp = combo_Fraction_mvm_in_celltage_dec{ii};
plot(temp','color',[0.8,0.8,0.8]);
temp = Fraction_mvm_in_celltage_dec(ii,:);
bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
xlim([0 10]);ylim([0 1]);
line([3.5 3.5],ylim,'color','k','linestyle',':');
line([6.5 6.5],ylim,'color','k','linestyle',':');
text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
ylabel('Fraction');

saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.png']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.pdf']); pause(1);   
    
    