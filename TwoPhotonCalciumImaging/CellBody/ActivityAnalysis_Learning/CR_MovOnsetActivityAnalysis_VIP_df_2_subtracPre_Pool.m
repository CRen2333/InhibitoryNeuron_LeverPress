%% 2P Data VIP
%% Get activity aligned to movement onset
clear all;
close all;
clc;

IN = 'VIP';
switch IN
    case 'VIP'
        Animals = {'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','WL_3526641-R','WL_3526642-L','WL_3526642-R'};
end
Fields = {'Field_1','Field_2','Field_3','Field_4'};

for curr_animal = 1:length(Animals)
    clearvars -except Animals IN curr_animal Fields
    Animal = Animals{curr_animal};
    LeverTrace_Path = fullfile('Z:\People\Chi\TwoP_IN',IN,Animal,'LeverTrace');
    Dates = dir(LeverTrace_Path);
    Dates = {Dates.name};
    Dates = sort(Dates(cellfun(@(x) contains(x, '17')||contains(x, '18')||contains(x, '19'), Dates)))';
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
                df_f_MovOnset_Alinged_EachROI_MeanAXTrial{curr_field,curr_date}(roi,:) = nanmean(df_f_MovOnset_Alinged{curr_field,curr_date}{roi,1},2);
            end
        end        
    end
    TargetPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    save([TargetPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','df_f_MovOnset_Alinged_EachROI_MeanAXTrial','-v7.3');
end

%% Averaged activity after movement onset, combine trials across sessions first then average
clear all;
close all;
clc;

IN = 'VIP';
Animals = {'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','WL_3526641-R','WL_3526642-L','WL_3526642-R'};

Fields = {'Field_1','Field_2','Field_3','Field_4'};

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    DataPath = fullfile('Z:\People\Chi\TwoP_IN\',IN,Animal,'df_f');
    load([DataPath filesep Animal '_MovOnsetAligendTraces.mat'],'df_f_MovOnset_Alinged','-mat');
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
        % Get 4 stages first
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,1) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,1)),2);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,2) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,2:4)),2);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,3) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,5:8)),2);
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,4) = nanmean(cell2mat(Post_MovOnset_Aligend_df_field{curr_animal}(curr_field,9:end)),2);
        [nan_index,~] = find(isnan(Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}));
        nan_cells{curr_animal,1}{curr_field,1} = [nan_cells{curr_animal,1}{curr_field,1}, nan_index];
        nan_cells{curr_animal,1}{curr_field,1} = unique(nan_cells{curr_animal,1}{curr_field,1});
        Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(nan_cells{curr_animal,1}{curr_field,1},:) = [];
        
        Pre_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,1) = nanmean(cell2mat(Pre_MovOnset_Aligend_df_field{curr_animal}(curr_field,1)),2);
        Pre_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,2) = nanmean(cell2mat(Pre_MovOnset_Aligend_df_field{curr_animal}(curr_field,2:4)),2);
        Pre_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,3) = nanmean(cell2mat(Pre_MovOnset_Aligend_df_field{curr_animal}(curr_field,5:8)),2);
        Pre_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(:,4) = nanmean(cell2mat(Pre_MovOnset_Aligend_df_field{curr_animal}(curr_field,9:end)),2);
        
        Pre_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1}(nan_cells{curr_animal,1}{curr_field,1},:) = [];
        
        NeuronNum{curr_animal}(curr_field) = size(Post_MovOnset_Aligend_df_field_stage{curr_animal}{curr_field,1},1);        
    end
    % Pool fields
    Post_MovOnset_Aligend_df_stage{curr_animal} = cell2mat(Post_MovOnset_Aligend_df_field_stage{curr_animal});
    Pre_MovOnset_Aligend_df_stage{curr_animal} = cell2mat(Pre_MovOnset_Aligend_df_field_stage{curr_animal});
    % Average
    Post_MovOnset_Aligend_df_stage_Mean(curr_animal,:) = nanmean(Post_MovOnset_Aligend_df_stage{curr_animal});
    Pre_MovOnset_Aligend_df_stage_Mean(curr_animal,:) = nanmean(Pre_MovOnset_Aligend_df_stage{curr_animal});
    clear df_f_MovOnset_Alinged
end

%% Pool neurons across animals
clear all
close all
clc

IN = 'VIP';
load(['Z:\People\Chi\TwoP_IN\' IN '\Figures_SubPre\VIP_ActivityAnalysis_SubPre.mat'],'Animals','Post_MovOnset_Aligend_df_stage','nan_cells','-mat');
% Post mvmOnset activity avreage across neurons
for curr_animal  = 1:length(Post_MovOnset_Aligend_df_stage)
    neuronum(curr_animal) = size(Post_MovOnset_Aligend_df_stage{curr_animal},1);
end
Animal_tag = [];
for ii = 1:length(neuronum)
    Animal_tag = [Animal_tag;ones(neuronum(ii),1)*ii];
end
All_cell_stage_conc = cell2mat(Post_MovOnset_Aligend_df_stage');

[CX] = cbrewer('div','PRGn',10);
color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
color_line = mean(color_value(2:3,:));
stages = {'Naive','Early','Middle','Late'};

figure; set(gcf,'position',[200,200,200,200]); hold on;
temp_var = All_cell_stage_conc;
temp_1 = nanmean(temp_var);
temp_2 = nanstd(temp_var)./sqrt(sum(~isnan((temp_var))));
plot(temp_1,'color',color_line,'LineWidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',color_line,'LineWidth',2);
end
xlim([0.5,4.5]); ylim([0,0.35]); xticks([1:4]); xticklabels(stages); ylabel('Mean df/f');
yticks([0:0.1:0.3]);
% random effect modal
Sessions = repmat([1:4],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = [];
for curr_animal = 1:length(Animals)
    Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal}))*curr_animal];
end
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
text(3,0.25,'p < 0.0001','color','k','HorizontalAlignment','center')
axis square

% Compare to naive stage
tbl.Sessions = nominal(tbl.Sessions);
glme = fitglme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[gbeta,gbetanames,gstats] = fixedEffects(glme);
PValues = gstats.pValue(2:end);
[FDR] = mafdr(PValues,'BHFDR', true);

FigurePath = ['Z:\People\Chi\TwoP_IN\' IN filesep 'Figures_SubPre'];
if ~exist(FigurePath)
    mkdir(FigurePath);
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_Post.pdf']); pause(1);

% Distribution
figure; set(gcf,'color','w','position',[100,100,600,600]); hold on;
for ii = 1:4
    Delta_df_against_N(:,ii) = All_cell_stage_conc(:,ii)-All_cell_stage_conc(:,1);
end
for ii = 1:4
    Delta_df_against_E(:,ii) = All_cell_stage_conc(:,ii)-All_cell_stage_conc(:,2);
end
for ii = 1:4
    Delta_df_against_M(:,ii) = All_cell_stage_conc(:,ii)-All_cell_stage_conc(:,3);
end
edges = [-2:0.1:2];
for ii = 2:4
    subplot(4,4,ii);
    temp_n = histcounts(Delta_df_against_N(:,ii),edges,'Normalization','probability');
    histogram(Delta_df_against_N(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(ii,:));
    ylim([0 0.25]); xlim([edges(1) edges(end)]);
    line([0 0],ylim,'color','k','linestyle',':');
    xlabel('delta df/f'); ylabel('Prob.'); axis square;
    temp_var = All_cell_stage_conc(:,[1,ii]);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal},1),2)*curr_animal];
    end
    Stages_test = repmat([0,1],size(temp_var,1),1);
    y = temp_var(:);
    Animals_test = Animals_test(:);
    Stages_test = Stages_test(:);
    tbl = table(Animals_test,Stages_test,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Stages_test = nominal(tbl.Stages_test);
    lme = fitlme(tbl,'y ~ 1 + Stages_test + (1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
    [FDR] = mafdr(PValues,'BHFDR', true);
end

for ii = 3:4
    subplot(4,4,ii+4);
    temp_n = histcounts(Delta_df_against_E(:,ii),edges,'Normalization','probability');
    histogram(Delta_df_against_N(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(ii,:));    
    ylim([0 0.25]); xlim([edges(1) edges(end)]);
    line([0 0],ylim,'color','k','linestyle',':');
    xlabel('delta df/f'); ylabel('Prob.'); axis square;
    temp_var = All_cell_stage_conc(:,[2,ii]);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal},1),2)*curr_animal];
    end
    Stages_test = repmat([0,1],size(temp_var,1),1);
    y = temp_var(:);
    Animals_test = Animals_test(:);
    Stages_test = Stages_test(:);
    tbl = table(Animals_test,Stages_test,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Stages_test = nominal(tbl.Stages_test);
    lme = fitlme(tbl,'y ~ 1 + Stages_test + (1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
end

for ii = 4
    subplot(4,4,ii+8);
    temp_n = histcounts(Delta_df_against_M(:,ii),edges,'Normalization','probability');
    histogram(Delta_df_against_N(:,ii),edges,'Normalization','probability','edgecolor','none','facecolor',color_value(ii,:));    
    ylim([0 0.25]); xlim([edges(1) edges(end)]);
    line([0 0],ylim,'color','k','linestyle',':');
    xlabel('delta df/f'); ylabel('Prob.'); axis square;
    temp_var = All_cell_stage_conc(:,[3,ii]);
    Animals_test = [];
    for curr_animal = 1:length(Animals)
        Animals_test = [Animals_test;ones(size(Post_MovOnset_Aligend_df_stage{curr_animal},1),2)*curr_animal];
    end
    Stages_test = repmat([0,1],size(temp_var,1),1);
    y = temp_var(:);
    Animals_test = Animals_test(:);
    Stages_test = Stages_test(:);
    tbl = table(Animals_test,Stages_test,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    tbl.Stages_test = nominal(tbl.Stages_test);
    lme = fitlme(tbl,'y ~ 1 + Stages_test + (1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    pValue = stats.pValue;
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Delta_df_distribution.pdf']); pause(1);

% Temporal traces
load(['Z:\People\Chi\TwoP_IN\' IN '\Figures_SubPre\' IN '_ActivityAnalysis_SubPre.mat'],'df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages','mycmap','MovOnset_Frame','-mat');
for ii = 1:4
    all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages(:,ii));
end
% sort based on naive
[~,I_naive] = sort(All_cell_stage_conc(:,1),'descend');
% sort within each stage
figure; hold on; set(gcf,'color','w','position',[1500 50 800 500]);
Stages = {'Naive','Early','Middle','Late'};
warning off
color_map = cbrewer('div','RdBu',64);
warning on
color_map = flipud(color_map);
for ii = 1:4
    subplot(1,4,ii);
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    [~,I] = sort(sum(temp(:,MovOnset_Frame:end),2),'descend');
    temp = temp(I,:);
%     imagesc(temp,[-2.6 2.6]); colormap(color_map);
    imagesc(temp,[-0.7 2.6]); colormap(mycmap);
%     imagesc(all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii}(I_naive,:),[-0.5 2.5]);
    line([8,8],ylim,'color','k','linestyle',':');
    ylabel('# of cells');
    xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (s)');
    title(Stages{ii});
end
% mycmap = colormap(gca);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MeanXTrial_dfof_conc.pdf']); pause(1);

figure; hold on; set(gcf,'color','w','position',[200 200 200 200]);
for ii = 1:4
    temp = all_df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
ylim([-0.1 0.6]);
line([8,8],ylim,'color','k','linestyle',':','linewidth',1);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Mean_dfof_temporal.pdf']); pause(1);

% significantly changed fraction
load('Z:\People\Chi\TwoP_IN\VIP\Figures_SubPre\VIP_ActivityAnalysis_SubPre.mat', 'Delta_df_against_N_CellTag')
All_Delta_df_against_N_CellTag = [];
for ii = 1:length(Animals)
    All_Delta_df_against_N_CellTag = [All_Delta_df_against_N_CellTag;cell2mat(Delta_df_against_N_CellTag{ii})];
end
figure; hold on; set(gcf,'color','w','position',[50 50 800 200]);
for ii = 2:4
    subplot(1,4,ii);
    pie([sum(All_Delta_df_against_N_CellTag(:,ii)==1),sum(All_Delta_df_against_N_CellTag(:,ii)==-1),sum(All_Delta_df_against_N_CellTag(:,ii)==0)]);
    axis square;
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.png']); pause(1);
print([FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Pie.pdf'],'-dpdf','-bestfit'); pause(1);

% stack bar
for ii = 1:4
    temp_bar(ii,:) = [sum(All_Delta_df_against_N_CellTag(:,ii)==1),sum(All_Delta_df_against_N_CellTag(:,ii)==-1),sum(All_Delta_df_against_N_CellTag(:,ii)==0)];
end
figure; hold on; set(gcf,'color','w','position',[400 400 200 200]);
bar(temp_bar(2:4,:)/sum(temp_bar(1,:))*100,'stacked','barwidth',0.7);
xlim([0.5,3.5]);ylim([0 100]); % [1,0.7,0.1;0.3,0.75,0.93;0.8,0.8,0.8]
xticks([1:3]);xticklabels(Stages(2:4));ylabel('Fraction');
axis square;
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Bar.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Bar.png']); pause(1);
print([FigurePath filesep 'Pool_' IN '_Deltadf_L_N_Bar.pdf'],'-dpdf','-bestfit'); pause(1);

clear pval
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('M',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,3)];
[tbl,~,pval(1),~] = crosstab(x1,x2);
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,2);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(2),~] = crosstab(x1,x2);
x1 = [repmat('M',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_Delta_df_against_N_CellTag(:,3);All_Delta_df_against_N_CellTag(:,4)];
[tbl,~,pval(3),~] = crosstab(x1,x2);

% randomly select 50% animals
animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_matrix = [];
    for jj = 1:length(curr_combo)
        temp_matrix = [temp_matrix;cell2mat(Delta_df_against_N_CellTag{curr_combo(jj)})];
    end
    for ii = 2:4
        Fraction_CellTag_combo{ii}(ii_combo,1) = sum(temp_matrix(:,ii)==1)/length(temp_matrix(:,ii));
        Fraction_CellTag_combo{ii}(ii_combo,2) = sum(temp_matrix(:,ii)==-1)/length(temp_matrix(:,ii));
        Fraction_CellTag_combo{ii}(ii_combo,3) = sum(temp_matrix(:,ii)==0)/length(temp_matrix(:,ii)); 
    end
end
figure; hold on; set(gcf,'color','w','position',[50 50 800 200]);
for ii = 2:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_CellTag_combo{ii};
    plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
    plot(temp_bar(ii,1)/length(All_Delta_df_against_N_CellTag)*100,temp_bar(ii,2)/length(All_Delta_df_against_N_CellTag)*100,'marker','x','color','k');
    xlim([0 60]);ylim([0 80]);
    plot([0,60],[0,60],'linestyle',':','color','k');
    xlabel('Inc. fraction');ylabel('Dec. fraction'); axis square;
    xlim([0 60]);xticks([0:20:60]);
    title(Stages{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie.pdf'],'-dpdf','-bestfit'); pause(1);        

% increased fraction - decreased fraction 
figure; hold on; set(gcf,'color','w','position',[50 50 800 200]);
for ii = 2:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_CellTag_combo{ii};
    temp_var_1 = Fraction_CellTag_combo{ii}(:,1)-Fraction_CellTag_combo{ii}(:,2);
    shake_x = (rand(size(temp_var_1))-0.5)*0.3;
    plot(shake_x,temp_var_1,'color',[0.5 0.5 0.5],'Linestyle','none','marker','.','markersize',8);
    xlim([-0.5 0.5]);xticks([0]);xticklabels({'inc.-dec.'});
    ylim([-0.6 0]);
    line(xlim,[0 0],'color','k','linestyle',':');
    ylabel('Fraction');
    title(stages{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Deltadf_L_N_Pie_incdec.pdf'],'-dpdf','-bestfit'); pause(1);        

%% Fraction of MVM modulated neurons
% based on stages, similar as loading from SOM_ActivityAnalysis_SubPre
clear movement_activated_cells_index_stages movement_suppressed_cells_index_stages
for ii = 1:length(Animals)
    Animal = Animals{ii};
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal '\df_f' filesep Animal '_Activity_bootstats_stage.mat']);
    for jj = 1:length(Label.MvmModuLabel_Act)
        Label.MvmModuLabel_Act{jj}(nan_cells{ii}{jj},:) = [];
        Label.MvmModuLabel_Sup{jj}(nan_cells{ii}{jj},:) = [];
    end
    movement_activated_cells_index_stages{ii,1} = logical(cell2mat(Label.MvmModuLabel_Act));
    movement_suppressed_cells_index_stages{ii,1} = logical(cell2mat(Label.MvmModuLabel_Sup));
end
All_movement_activated_cells_index_stages = cell2mat(movement_activated_cells_index_stages);
All_movement_suppressed_cells_index_stages = cell2mat(movement_suppressed_cells_index_stages);
All_movement_modulated_cells_index_stages = All_movement_activated_cells_index_stages-All_movement_suppressed_cells_index_stages;

for ii = 1:4
    temp_bar(ii,:) = [sum(All_movement_modulated_cells_index_stages(:,ii)==1),sum(All_movement_modulated_cells_index_stages(:,ii)==-1),...
        sum(All_movement_modulated_cells_index_stages(:,ii)==0)];
end
figure; hold on; set(gcf,'color','w','position',[400 400 200 200]);
bar(temp_bar(1:4,:)/sum(temp_bar(1,:))*100,'stacked','barwidth',0.7);
xlim([0.5,4.5]);ylim([0 100]); % [0.75,0,0;0,0.34,0.42;0.8,0.8,0.8]
xticks([1:4]);xticklabels(Stages);ylabel('Fraction');
axis square;

clear pval
x1 = [repmat('N',sum(temp_bar(1,:)),1);repmat('E',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,1);All_movement_modulated_cells_index_stages(:,2)];
[tbl,~,pval(1),~] = crosstab(x1,x2);
x1 = [repmat('N',sum(temp_bar(1,:)),1);repmat('M',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,1);All_movement_modulated_cells_index_stages(:,3)];
[tbl,~,pval(2),~] = crosstab(x1,x2);
x1 = [repmat('N',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,1);All_movement_modulated_cells_index_stages(:,4)];
[tbl,~,pval(3),~] = crosstab(x1,x2);
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('M',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,2);All_movement_modulated_cells_index_stages(:,3)];
[tbl,~,pval(4),~] = crosstab(x1,x2);
x1 = [repmat('E',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,2);All_movement_modulated_cells_index_stages(:,4)];
[tbl,~,pval(5),~] = crosstab(x1,x2);
x1 = [repmat('M',sum(temp_bar(1,:)),1);repmat('L',sum(temp_bar(1,:)),1)];
x2 = [All_movement_modulated_cells_index_stages(:,3);All_movement_modulated_cells_index_stages(:,4)];
[tbl,~,pval(6),~] = crosstab(x1,x2);
[~,I] = sort(pval,'ascend');
for ii = 1:6
    pval_corrected_threshold(ii) = ii*0.05/6;
end
pval_corrected_threshold = pval_corrected_threshold(I);
pval < pval_corrected_threshold
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_Bar.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_Bar.png']); pause(1);
print([FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup_Bar.pdf'],'-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','position',[500 500 800 200]);
for ii = 1:4
    subplot(1,4,ii);
    temp = All_movement_activated_cells_index_stages(:,ii)-All_movement_suppressed_cells_index_stages(:,ii);
    pie([sum(temp==1),sum(temp==-1),sum(temp==0)]);
%     clear temp;
    axis square;
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.png']); pause(1);
print([FigurePath filesep 'Pool_' IN '_Fraction_MvmActSup.pdf'],'-dpdf','-bestfit'); pause(1);
% compare between stages
x1 = [repmat('E',size(All_movement_activated_cells_index_stages,1),1);
    repmat('M',size(All_movement_activated_cells_index_stages,1),1);
    repmat('L',size(All_movement_activated_cells_index_stages,1),1);];
x2 = All_movement_activated_cells_index_stages(:,2:4)-All_movement_suppressed_cells_index_stages(:,2:4);
x2 = x2(:);
[tbl,chi2stat,pval] = crosstab(x1,x2);

animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_matrix = [];
    for jj = 1:length(curr_combo)
        temp_matrix = [temp_matrix;movement_activated_cells_index_stages{curr_combo(jj)}-movement_suppressed_cells_index_stages{curr_combo(jj)}];
    end
    for ii = 1:4
        Fraction_MvmMod_combo{ii}(ii_combo,1) = sum(temp_matrix(:,ii)==1)/length(temp_matrix(:,ii));
        Fraction_MvmMod_combo{ii}(ii_combo,2) = sum(temp_matrix(:,ii)==-1)/length(temp_matrix(:,ii));
        Fraction_MvmMod_combo{ii}(ii_combo,3) = sum(temp_matrix(:,ii)==0)/length(temp_matrix(:,ii)); 
    end
end
figure; hold on; set(gcf,'color','w','position',[500 500 800 200]);
for ii = 1:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii};
    plot(temp_var(:,1)*100,temp_var(:,2)*100,'linestyle','none','marker','.','color',[0.5 0.5 0.5]);
    plot(temp_bar(ii,1)/length(All_movement_modulated_cells_index_stages)*100,temp_bar(ii,2)/length(All_movement_modulated_cells_index_stages)*100,'marker','x','color','k');
    xlim([0 80]);ylim([0 80]);
    plot([0,80],[0,80],'linestyle',':','color','k');
    xlabel('Act. fraction');ylabel('Sup. fraction'); axis square;
    xticks([0:20:80]);
    title(Stages{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmActSup.pdf'],'-dpdf','-bestfit'); pause(1);        

figure; hold on; set(gcf,'color','w','position',[500 500 800 200]);
for ii = 1:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii};
    temp_var_1 = Fraction_MvmMod_combo{ii}(:,1)-Fraction_MvmMod_combo{ii}(:,2);
    shake_x = (rand(size(temp_var_1))-0.5)*0.3;
    plot(shake_x,temp_var_1,'color',[0.5 0.5 0.5],'Linestyle','none','marker','.','markersize',8);
    xlim([-0.5 0.5]);xticks([0]);xticklabels({'act.-sup.'});
    line(xlim,[0 0],'color','k','linestyle',':');
    ylabel('Fraction');
    title(stages{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmAct_minus_Sup.pdf'],'-dpdf','-bestfit'); pause(1);        

figure; hold on; set(gcf,'color','w','position',[500 500 800 200]);
for ii = 2:4
    subplot(1,4,ii); hold on;
    temp_var = Fraction_MvmMod_combo{ii}-Fraction_MvmMod_combo{1};
    temp_var_1 = temp_var(:,1);
    temp_var_2 = temp_var(:,2);
    shake_x = (rand(size(temp_var_1))-0.5)*0.3;
    plot(shake_x-0.5,temp_var_1,'color',[191,0,0]/255,'Linestyle','none','marker','.','markersize',8);
    plot(shake_x+0.5,temp_var_2,'color',[0,87,107]/255,'Linestyle','none','marker','.','markersize',8);
    xlim([-1 1]);xticks([-0.5 0.5]);xticklabels({'act.','sup.'});
    line(xlim,[0 0],'color','k','linestyle',':');
    ylabel('Fraction');
    title(stages{ii});
    axis square;
end
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_Fraction_MvmMod_minus_Naive.pdf'],'-dpdf','-bestfit'); pause(1);        

% MVM modulated neurons activity onset
load('Z:\People\Chi\TwoP_IN\VIP\Figures_SubPre\VIP_ActivityAnalysis_SubPre.mat','movement_activated_cells_onset_stage','movement_suppressed_cells_onset_stage','-mat');
for ii = 1:4
    All_movement_activated_cells_onset_stage{ii} = cell2mat(movement_activated_cells_onset_stage(:,ii));
    All_movement_suppressed_cells_onset_stage{ii} = cell2mat(movement_suppressed_cells_onset_stage(:,ii));
end

figure; hold on; set(gcf,'color','w','position',[50 50 200 400]);
subplot(2,1,1); hold on;
for ii = 1:4
    cellnum = length(All_movement_activated_cells_onset_stage{ii});
    plot(sort(All_movement_activated_cells_onset_stage{ii},'ascend'),[1:cellnum]/cellnum,'color',color_value(ii,:),'linewidth',1);
end
line([0 0],ylim,'color','k','linestyle',':');
xlim([-0.5 2]);
xlabel('Activity onset (s)'); ylabel('cum. prob.');
title('Activated'); axis square;
subplot(2,1,2); hold on;
for ii = 1:4
    cellnum = length(All_movement_suppressed_cells_onset_stage{ii});
    plot(sort(All_movement_suppressed_cells_onset_stage{ii},'ascend'),[1:cellnum]/cellnum,'color',color_value(ii,:),'linewidth',1);
end
line([0 0],ylim,'color','k','linestyle',':');
xlim([-0.5 2]);
xlabel('Activity onset (s)'); ylabel('cum. prob.');
title('Suppressed'); axis square;
   
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.png']); pause(1);
print([FigurePath filesep 'comb_Pool_' IN '_MvmMod_ActivityOnset.pdf'],'-dpdf','-bestfit'); pause(1);        

% Check trace
for ii = 1:4
    movement_activated_cells_trace_pool{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages(:,ii));
    movement_activated_cells_trace_pool{ii} = movement_activated_cells_trace_pool{ii}(logical(All_movement_activated_cells_index_stages(:,ii)),:);
end
% temp_index = cell2mat(movement_suppressed_cells_index_stages');
for ii = 1:4
    movement_suppressed_cells_trace_pool{ii} = cell2mat(df_f_MovOnset_Alinged_EachROI_MeanAXTrial_stages(:,ii));
    movement_suppressed_cells_trace_pool{ii} = movement_suppressed_cells_trace_pool{ii}(logical(All_movement_suppressed_cells_index_stages(:,ii)),:);
end
% across neurons
figure; hold on; set(gcf,'color','w','position',[200 200 200 400]);
subplot(2,1,1); hold on;
for ii = 1:4
    temp = movement_activated_cells_trace_pool{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
line([8,8],ylim,'color','k','linestyle',':','linewidth',1); ylim([-0.1, 0.8]);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square; title('Activated');
subplot(2,1,2); hold on;
for ii = 1:4
    temp = movement_suppressed_cells_trace_pool{ii};
    temp_1 = nanmean(temp);
    temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
    h = area([(temp_1-temp_2)',(2*temp_2)']);
    set(h(1),'EdgeColor','none','FaceColor','none');
    set(h(2),'EdgeColor','none','FaceColor',color_value(ii,:),'FaceAlpha',0.3);
    plot(temp_1,'color',color_value(ii,:),'linewidth',2);
end
line([8,8],ylim,'color','k','linestyle',':','linewidth',1); ylim([-0.6, 0.1]);
xlim([1,36]); xticks([8,22,36]); xticklabels({'0','1','2'}); xlabel('Time (sec)');
ylabel('Mean df/f'); axis square; title('Suppressed');
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_MvmActSup_trace.pdf']); pause(1);

%% Check relationship between significalty changed activity and mvm modulated
% mvm in celltag
temp_1 = All_Delta_df_against_N_CellTag;
temp_2 = logical(All_movement_activated_cells_index_stages);
temp_3 = logical(All_movement_suppressed_cells_index_stages);
temp_4 = ones(size(temp_1));
temp_4 = temp_4-logical(temp_2+temp_3);
for ii = 2:4
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

figure; set(gcf,'color','w','position',[2000 100 600 400]); hold on
for ii = 1:3
    subplot(3,3,ii); hold on
    temp = Fraction_mvm_in_celltage_inc(ii+1,:);
    bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
    xlim([0 10]);ylim([0 1]);
    line([3.5 3.5],ylim,'color','k','linestyle',':');
    line([6.5 6.5],ylim,'color','k','linestyle',':');
    text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
    xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
    ylabel('Fraction');
end
for ii = 1:3
    subplot(3,3,ii+3); hold on
    temp = Fraction_mvm_in_celltage_dec(ii+1,:);
    bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
    xlim([0 10]);ylim([0 1]);
    line([3.5 3.5],ylim,'color','k','linestyle',':');
    line([6.5 6.5],ylim,'color','k','linestyle',':');
    text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
    xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
    ylabel('Fraction');
end    
for ii = 1:3
    subplot(3,3,ii+6); hold on
    temp = Fraction_mvm_in_celltage_no(ii+1,:);
    bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
    xlim([0 10]);ylim([0 1]);
    line([3.5 3.5],ylim,'color','k','linestyle',':');
    line([6.5 6.5],ylim,'color','k','linestyle',':');
    text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
    xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
    ylabel('Fraction');
end
saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.png']); pause(1);
saveas(gcf, [FigurePath filesep 'Pool_' IN '_fraction_mvm_in_celltag.pdf']); pause(1);   

animal_combo = nchoosek([1:length(Animals)],round(length(Animals)/2));
for ii_combo = 1:size(animal_combo,1)
    curr_combo = animal_combo(ii_combo,:);
    temp_1 = [];
    for jj = 1:length(curr_combo)
        temp_1 = [temp_1;cell2mat(Delta_df_against_N_CellTag{curr_combo(jj)})];
    end
    temp_2 = [];
    for jj = 1:length(curr_combo)
        temp_2 = [temp_2;movement_activated_cells_index_stages{curr_combo(jj)}];
    end
    temp_3 = [];
    for jj = 1:length(curr_combo)
        temp_3 = [temp_3;movement_suppressed_cells_index_stages{curr_combo(jj)}];
    end
    temp_2 = logical(temp_2);
    temp_3 = logical(temp_3);
    temp_4 = ones(size(temp_1));
    temp_4 = temp_4-logical(temp_2+temp_3);
    for ii = 2:4
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

figure; set(gcf,'color','w','position',[2000 100 600 400]); hold on
for ii = 1:3
    subplot(3,3,ii); hold on
    temp = combo_Fraction_mvm_in_celltage_inc{ii+1};
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
end
for ii = 1:3
    subplot(3,3,ii+3); hold on
    temp = combo_Fraction_mvm_in_celltage_dec{ii+1};
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
end    
for ii = 1:3
    subplot(3,3,ii+6); hold on
    temp = combo_Fraction_mvm_in_celltage_no{ii+1};
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
end
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.png']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag.pdf']); pause(1);   

% decrease only    
figure; set(gcf,'color','w','position',[2000 100 600 400]); hold on
for ii = 1:3
    subplot(3,3,ii+3); hold on
    temp = combo_Fraction_mvm_in_celltage_dec{ii+1};
    plot(temp','color',[0.8,0.8,0.8]);
    temp = Fraction_mvm_in_celltage_dec(ii+1,:);
    bar([1:9],temp,'facecolor','none','edgecolor','k','barwidth',0.6);
    xlim([0 10]);ylim([0 1]);
    line([3.5 3.5],ylim,'color','k','linestyle',':');
    line([6.5 6.5],ylim,'color','k','linestyle',':');
    text(2,0.9,'Act.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(5,0.9,'Sup.','color','k','fontsize',8,'HorizontalAlignment','center');
    text(8,0.9,'no','color','k','fontsize',8,'HorizontalAlignment','center');
    xticks([1:9]);xticklabels({'a.','s.','n.','a.','s.','n.','a.','s.','n.'});
    ylabel('Fraction');
end 

saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.fig']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.png']); pause(1);
saveas(gcf, [FigurePath filesep 'combo_Pool_' IN '_fraction_mvm_in_celltag_dec.pdf']); pause(1);   
