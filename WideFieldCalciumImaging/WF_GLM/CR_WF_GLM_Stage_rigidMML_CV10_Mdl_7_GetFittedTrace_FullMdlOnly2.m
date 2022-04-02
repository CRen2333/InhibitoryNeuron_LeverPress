% Reconstruct activity based on beta from full model, cv10
% Calculate correlation for cued rewarded movement epoch
close all
clc

clearvars -except VIP SOM PV Stats ColorMap FullOnly

Model = 'cv10_Model_7';
cd(['Z:\People\Chi\WFLP_IN\GLM\' Model]);

INs = {'VIP','SOM','PV'};
cv_num = 10;
FullOnly = true;

for curr_IN = 1:length(INs)
    
    clearvars -except VIP SOM PV INs curr_IN Model Stats ColorMap cv_num FullOnly
    IN = INs{curr_IN};
    Initial = 'CR';
    
    switch IN
        case 'VIP'
            Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};
        case 'SOM'
            Animals = {'3183958-L','3183958-R','3183959-LL','3218181-L','3218181-O','3218181-R','3438521-O','3438521-L','3438521-R','3438521-LR','3453262-O','3453262-L'};
        case 'PV'
            Animals = {'3161016-O','3161018-R','3233232-L','3233232-O','3233232-R','3233233-O','3233233-L','3233233-R','3491479-L','3491479-LR','3547207-LR'};
    end


    ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
    Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
    for ii = 1:length(ROI)
        temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
        reordered_module(ii) = find(temp==1);
        clear temp
    end

    for curr_animal = 1:length(Animals)
        tic
        Animal = Animals{curr_animal};
        disp([IN ' ' Animal]);
        data_folder = ['Z:\People\Chi\WFLP_IN\GLM\' Model];
        if FullOnly
            file_name = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_rigidMML_CV10_FullOnly.mat'];
        else
            file_name = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_rigidMML_CV10.mat'];
        end
        load(file_name,'CV10','shuffle_CV10','predictors_tag','-mat');
        load([data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_rigidMML_CV10.mat'],'predictors_all');
        
        % Organize MvmOnset to stage
        General_path = fullfile('Z:\People\Chi\WFLP_IN',IN,[Initial '_' Animal]);
        load([General_path filesep 'EventAligned_Gap500' filesep Initial '_' Animal '_FrameIndexForGLM.mat'],'FrameIndex');
        load(fullfile(General_path,'MovAnalysis',[Initial '_' Animal '_FrameIndex_Info.mat']),'Imaging_Day','IndexInfo');
        Imaging_Day = Imaging_Day(1:11);
        Image_Range_MovOnset = [-15:60];
        for curr_session = 1:length(Imaging_Day)
            RwdMvmOnset{curr_session,1} = [];
            if isempty(IndexInfo{curr_session})
                RwdMvmOnset{curr_session,1}{curr_block,1} = [];
                continue
            end
            for curr_block = 1:length(IndexInfo{curr_session})
                Temp = IndexInfo{curr_session}{curr_block};
                if isempty(Temp)
                    RwdMvmOnset{curr_session,1}{curr_block,1} = [];
                    continue
                end
                Temp = Temp(~isnan(Temp(:,6)),:); % rewarded
                Temp = Temp(~isnan(Temp(:,3)),:); % cued
                Temp = round(Temp(:,7)*29.98);
                Temp(Temp+Image_Range_MovOnset(end)>9000,:) = [];
                Temp(Temp-Image_Range_MovOnset(1)<0,:) = [];
                RwdMvmOnset{curr_session,1}{curr_block,1} = Temp;
            end
        end
        
        Stage_sessions = {[1],[2,3,4],[5,6,7,8],[9,10,11]};
        Stages = {'Naive','Early','Middle','Late'};
        for ii = 1:4
            block_count = 0;
            RwdMvmOnset_Stage{ii} = [];
            for jj = 1:length(Stage_sessions{ii})
                curr_session = Stage_sessions{ii}(jj);
                for kk = 1:length(RwdMvmOnset{curr_session,1})
                    if isempty(FrameIndex.MvmOn{curr_session}{kk})
                        continue
                    end
                    if ~isempty(RwdMvmOnset{curr_session,1}{kk,1})
                        RwdMvmOnset_Stage{ii} = [RwdMvmOnset_Stage{ii};RwdMvmOnset{curr_session,1}{kk,1}+block_count*9000];
                        block_count = block_count+1;
                    else
                        block_count = block_count+1;
                    end
                end
            end
        end
        
        for ii = 1:4
            PredictorsStage{ii} = cell2mat(predictors_all(Stage_sessions{ii}));
        end
                
        % Full and true first      
        field = 'Full';
        for ii = 1:4
            for cv_ii = 1:10
                training_index = true(1,cv_num*size(CV10.FitTrace_test{ii},1));
                test_size = length(training_index)/cv_num;
                CV10_FitTest.(field){ii}([(cv_ii-1)*test_size+1:cv_ii*test_size],:) = CV10.FitTrace_test{ii}(:,:,cv_ii);
            end
        end
        for ii = 1:4
            for cv_ii = 1:10
                training_index = true(1,cv_num*size(CV10.TrueTrace_test{ii},1));
                test_size = length(training_index)/cv_num;
                CV10_FitTest.True{ii}([(cv_ii-1)*test_size+1:cv_ii*test_size],:) = CV10.TrueTrace_test{ii}(:,:,cv_ii);
            end
        end
        
        % each predictor
        fields = {'Mvm','Lick','Cue','Reward'};
        for field_ii = 1:length(fields)
            field = fields{field_ii};
            predictor_index = contains(predictors_tag,field);
            for ii = 1:4
                curr_predictor = PredictorsStage{ii}(predictor_index,:);
                for cv_ii = 1:cv_num 
                    test_size = size(curr_predictor,2)/cv_num;
                    curr_beta = CV10.Betas{ii}(predictor_index,:,cv_ii);
                    curr_bias = CV10.Bias{ii}(:,:,cv_ii);
                    CV10_FitTest.(field){ii}([(cv_ii-1)*test_size+1:cv_ii*test_size],:) = ...
                        (curr_beta'*curr_predictor(:,[(cv_ii-1)*test_size+1:cv_ii*test_size])+curr_bias')';
                end
            end
        end
                    
        fields = fieldnames(CV10_FitTest);
        for field_ii = 1:length(fields)
            field = fields{field_ii};
            for roi = 1:16
                for ii = 1:4
                    for mvmonset_ii = 1:length(RwdMvmOnset_Stage{ii})
                        curr_onset = RwdMvmOnset_Stage{ii}(mvmonset_ii)+Image_Range_MovOnset;
                        CV10_FitTest_RMD.(field){roi}{ii}(mvmonset_ii,:) = CV10_FitTest.(field){ii}(curr_onset,roi);
                        CV10_FitTest_RMD_sub.(field){roi}{ii}(mvmonset_ii,:) = CV10_FitTest_RMD.(field){roi}{ii}(mvmonset_ii,:)-nanmean(CV10_FitTest_RMD.(field){roi}{ii}(mvmonset_ii,5:9));
                    end
                    CV10_FitTest_RMD_sub_mean.(field){roi}{ii}(curr_animal,:) = nanmean(CV10_FitTest_RMD_sub.(field){roi}{ii},1);
                end
            end
        end
                
        % Calculate corr only during cued rewarded movement        
        fields = fieldnames(CV10_FitTest);
        for field_ii = 1:length(fields)
            field = fields{field_ii};
            for roi = 1:16
                for ii = 1:4
                    temp_1 = CV10_FitTest_RMD.(field){roi}{ii}';
                    temp_1 = temp_1(:);
                    temp_2 = CV10_FitTest_RMD.True{roi}{ii}';
                    temp_2 = temp_2(:);
                    temp_corr = corrcoef(temp_1,temp_2);
                    CV10_FitTest_RMD_corr.(field)(roi,ii) = temp_corr(1,2);
                    CV10_FitTest_RMD_corr_all.(field){roi}(curr_animal,ii) = temp_corr(1,2);
                end
            end
        end
        savefile_name = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_rigidMML_CV10_FullBetaOnly2.mat'];
        save(savefile_name,'CV10_FitTest','CV10_FitTest_RMD','CV10_FitTest_RMD_sub','CV10_FitTest_RMD_corr','-v7.3');
        clear CV10 CV10_FitTest CV10_FitTest_RMD CV10_FitTest_RMD_sub shuffle_CV10 shuffle_CV10_FitTest shuffle_CV10_FitTest_RMD shuffle_CV10_FitTest_RMD_sub
        toc
    end
    switch IN
        case 'VIP'
            VIP.CV10_FitTest_RMD_sub_mean = CV10_FitTest_RMD_sub_mean;
            VIP.CV10_FitTest_RMD_corr_all = CV10_FitTest_RMD_corr_all;
        case 'SOM'
            SOM.CV10_FitTest_RMD_sub_mean = CV10_FitTest_RMD_sub_mean;
            SOM.CV10_FitTest_RMD_corr_all = CV10_FitTest_RMD_corr_all;
        case 'PV'
            PV.CV10_FitTest_RMD_sub_mean = CV10_FitTest_RMD_sub_mean;
            PV.CV10_FitTest_RMD_corr_all = CV10_FitTest_RMD_corr_all;
    end
    clear CV10_FitTest_RMD_sub_mean CV10_FitTest_RMD_corr_all shuffle_CV10_FitTest_RMD_sub_mean shuffle_CV10_FitTest_RMD_corr_all     
end

save(['GLM_' Model '_FitTraceRwdMVM_FullBetaOnly2.mat'],'VIP','SOM','PV','-v7.3');


%% load true trace as comparison
load('Z:\People\Chi\WFLP_IN\VIP\Craniotomy\GAP500LEEWAY150_THY1MASK\VIP_Activity_THY1MASK_OPExclude.mat', 'Cued_subtract_postActivity_all','Cued_subtract_prepActivity_all')
Stages = {'Naive','Early','Middle','Late'};
for curr_ROI = 1:16
    for curr_field = 1:length(Stages)
        field = Stages{curr_field};
        VIP.CV10_FitTest_RMD_sub_mean.TrueTrue{curr_ROI}{curr_field} = [Cued_subtract_prepActivity_all.(field){curr_ROI},Cued_subtract_postActivity_all.(field){curr_ROI}];
    end
end

load('Z:\People\Chi\WFLP_IN\SOM\Craniotomy\GAP500LEEWAY150_THY1MASK\SOM_Activity_THY1MASK_OPExclude.mat', 'Cued_subtract_postActivity_all','Cued_subtract_prepActivity_all')
Stages = {'Naive','Early','Middle','Late'};
for curr_ROI = 1:16
    for curr_field = 1:length(Stages)
        field = Stages{curr_field};
        SOM.CV10_FitTest_RMD_sub_mean.TrueTrue{curr_ROI}{curr_field} = [Cued_subtract_prepActivity_all.(field){curr_ROI},Cued_subtract_postActivity_all.(field){curr_ROI}];
    end
end

load('Z:\People\Chi\WFLP_IN\PV\Craniotomy\GAP500LEEWAY150_THY1MASK\PV_Activity_THY1MASK_OPExclude.mat', 'Cued_subtract_postActivity_all','Cued_subtract_prepActivity_all')
Stages = {'Naive','Early','Middle','Late'};
for curr_ROI = 1:16
    for curr_field = 1:length(Stages)
        field = Stages{curr_field};
        PV.CV10_FitTest_RMD_sub_mean.TrueTrue{curr_ROI}{curr_field} = [Cued_subtract_prepActivity_all.(field){curr_ROI},Cued_subtract_postActivity_all.(field){curr_ROI}];
    end
end
% 
save(['GLM_' Model '_FitTraceRwdMVM_FullBetaOnly2.mat'],'VIP','SOM','PV','-append');

% Correlation 
aaa = [];
for ii = 1:16
    aaa(:,ii) = nanmean(PV.CV10_FitTest_RMD_corr_all.Full{ii},2);
end
aaa = nanmean(aaa,2);
nanmean(aaa)
nanstd(aaa)/sqrt(length(aaa))

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
fields = fieldnames(VIP.CV10_FitTest_RMD_sub_mean);
fields = {'Full','Mvm','Lick','Cue','Reward'}; % reorder
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 1:4
            temp_var = VIP.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 1:4
            temp_var = VIP.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.002 0.028]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'VIP_FitTrace_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'VIP_FitTrace_RwdMVM_FullBetaOnly2.png','png');

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 1:4
            temp_var = SOM.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 1:4
            temp_var = SOM.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.002 0.028]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'SOM_FitTrace_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'SOM_FitTrace_RwdMVM_FullBetaOnly2.png','png');

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 1:4
            temp_var = PV.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 1:4
            temp_var = PV.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.002 0.028]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'PV_FitTrace_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'PV_FitTrace_RwdMVM_FullBetaOnly2.png','png');

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
fields = fieldnames(VIP.CV10_FitTest_RMD_sub_mean);
fields = fields([1,6,7,4,5,3,8]); % reorder
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 2:4
            temp_var = VIP.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii}-...
                VIP.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 2:4
            temp_var = VIP.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii}-...
                VIP.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.007 0.007]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'VIP_FitTraceDelta_RwdMVM_FullBetaOnly.fig','fig');
saveas(gcf,'VIP_FitTraceDelta_RwdMVM_FullBetaOnly.png','png');

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 2:4
            temp_var = SOM.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii}-...
            SOM.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 2:4
            temp_var = SOM.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii}-...
                SOM.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.005 0.013]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'SOM_FitTraceDelta_RwdMVM_FullBetaOnly.fig','fig');
saveas(gcf,'SOM_FitTraceDelta_RwdMVM_FullBetaOnly.png','png');

figure; set(gcf,'color','w','position',[100 100 1500 500]); hold on;
colors_true = repmat([0.6;0.45;0.3;0.15],1,3);
colors_fit = cbrewer('seq','Blues',9);
colors_fit = colors_fit([5:8],:);
for roi = 1:16
    for curr_field = 1:length(fields)
        subplot(length(fields),16,(curr_field-1)*16+roi); hold on;
        field = fields{curr_field};
        for ii = 2:4
            temp_var = PV.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{ii}-...
                PV.CV10_FitTest_RMD_sub_mean.True{reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_true(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_true(ii,:));
        end
        for ii = 2:4
            temp_var = PV.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{ii}-...
                PV.CV10_FitTest_RMD_sub_mean.(field){reordered_module(roi)}{1};
            temp_1 = nanmean(temp_var);
            temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
%             h = area([(temp_1-temp_2)',(2*temp_2)']);
%             set(h(1),'EdgeColor','none','FaceColor','none');
%             set(h(2),'EdgeColor','none','FaceColor',colors_fit(ii,:),'FaceAlpha',0.3);                
            plot(temp_1,'color',colors_fit(ii,:));
        end
        ylim([-0.01 0.01]);xlim([1 76]);
        if curr_field == 1
            title(Ordered_ROI{roi}, 'fontsize',8);
        end
        if roi == 1
            ylabel(field);
        end
        xticks([]);
    end
end        
saveas(gcf,'PV_FitTraceDelta_RwdMVM_FullBetaOnly.fig','fig');
saveas(gcf,'PV_FitTraceDelta_RwdMVM_FullBetaOnly.png','png');

% plot example trace for r-M1
colors_field = [0,0,0;...
    63.75*1.5,63.75*1.5,63.75*1.5;...
    166,73,90;...
%     217,156,161;...
    60,94,115;...
%     124,158,166;...
    217,181,150;...
    166,120,93]/255;
Stages = {'Naive','Early','Middle','Late'};
fields = {'True','Full','Mvm','Lick','Cue','Reward'}; % reorder
figure; set(gcf,'color','w','position',[100 100 800 200]); hold on;
roi = 8;
for ii = 1:4
    subplot(1,4,ii); hold on;
    for field_ii = 1:length(fields)
        field = fields{field_ii};
        temp_var = VIP.CV10_FitTest_RMD_sub_mean.(field){roi}{ii};
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colors_field(field_ii,:),'FaceAlpha',0.3);                
        plot(temp_1,'color',colors_field(field_ii,:));
    end
    xlim([1 76]); xticks([16,46,76]); xticklabels([0,1,2]); xlabel('Time (s)');
    ylim([-0.002 0.015]); ylabel('df/f');
    line([16 16],ylim,'color','k','linestyle',':');
    title(Stages{ii});
    axis square;
end
saveas(gcf,'rM1_VIP_FitTrace_RwdMVM_FullBetaOnly2.fig','fig')
saveas(gcf,'rM1_VIP_FitTrace_RwdMVM_FullBetaOnly2.png','png')
saveas(gcf,'rM1_VIP_FitTrace_RwdMVM_FullBetaOnly2.pdf','pdf')

figure; set(gcf,'color','w','position',[100 100 800 200]); hold on;
roi = 8;
for ii = 1:4
    subplot(1,4,ii); hold on;
    for field_ii = 1:length(fields)
        field = fields{field_ii};
        temp_var = SOM.CV10_FitTest_RMD_sub_mean.(field){roi}{ii};
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colors_field(field_ii,:),'FaceAlpha',0.3);                
        plot(temp_1,'color',colors_field(field_ii,:));
    end
    xlim([1 76]); xticks([16,46,76]); xticklabels([0,1,2]); xlabel('Time (s)');
    ylim([-0.007 0.028]); ylabel('df/f');
    line([16 16],ylim,'color','k','linestyle',':');
    title(Stages{ii});
    axis square;
end
saveas(gcf,'rM1_SOM_FitTrace_RwdMVM_FullBetaOnly2.fig','fig')
saveas(gcf,'rM1_SOM_FitTrace_RwdMVM_FullBetaOnly2.png','png')
saveas(gcf,'rM1_SOM_FitTrace_RwdMVM_FullBetaOnly2.pdf','pdf')

figure; set(gcf,'color','w','position',[100 100 800 200]); hold on;
roi = 8;
for ii = 1:4
    subplot(1,4,ii); hold on;
    for field_ii = 1:length(fields)
        field = fields{field_ii};
        temp_var = PV.CV10_FitTest_RMD_sub_mean.(field){roi}{ii};
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var)));
        h = area([(temp_1-temp_2)',(2*temp_2)']);
        set(h(1),'EdgeColor','none','FaceColor','none');
        set(h(2),'EdgeColor','none','FaceColor',colors_field(field_ii,:),'FaceAlpha',0.3);                
        plot(temp_1,'color',colors_field(field_ii,:));
    end
    xlim([1 76]); xticks([16,46,76]); xticklabels([0,1,2]); xlabel('Time (s)');
    ylim([-0.003 0.022]); ylabel('df/f');
    line([16 16],ylim,'color','k','linestyle',':');
    title(Stages{ii});
    axis square;
end
saveas(gcf,'rM1_PV_FitTrace_RwdMVM_FullBetaOnly2.fig','fig')
saveas(gcf,'rM1_PV_FitTrace_RwdMVM_FullBetaOnly2.png','png')
saveas(gcf,'rM1_PV_FitTrace_RwdMVM_FullBetaOnly2.pdf','pdf')

% Get mean over 2 sec and norm to true trace
fields = fieldnames(VIP.CV10_FitTest_RMD_sub_mean);
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        for ii = 1:4
            VIP.CV10_PostActivity.(field){roi}(:,ii) = nanmean(VIP.CV10_FitTest_RMD_sub_mean.(field){roi}{ii}(:,16:end),2);
        end
        VIP.CV10_PostActivity_delta.(field){roi} = VIP.CV10_PostActivity.(field){roi}-repmat(VIP.CV10_PostActivity.(field){roi}(:,1),1,4);
    end
end
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        temp_matrix = nanmean(VIP.CV10_PostActivity_delta.True{roi});
        VIP.CV10_PostActivity_deltaNorm.(field){roi} = VIP.CV10_PostActivity_delta.(field){roi}./repmat(temp_matrix,size(VIP.CV10_PostActivity_delta.(field){roi},1),1);
    end
end
fields = fieldnames(SOM.CV10_FitTest_RMD_sub_mean);
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        for ii = 1:4
            SOM.CV10_PostActivity.(field){roi}(:,ii) = nanmean(SOM.CV10_FitTest_RMD_sub_mean.(field){roi}{ii}(:,16:end),2);
        end
        SOM.CV10_PostActivity_delta.(field){roi} = SOM.CV10_PostActivity.(field){roi}-repmat(SOM.CV10_PostActivity.(field){roi}(:,1),1,4);
    end
end
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        temp_matrix = nanmean(SOM.CV10_PostActivity_delta.True{roi});
        SOM.CV10_PostActivity_deltaNorm.(field){roi} = SOM.CV10_PostActivity_delta.(field){roi}./repmat(temp_matrix,size(SOM.CV10_PostActivity_delta.(field){roi},1),1);
    end
end
fields = fieldnames(PV.CV10_FitTest_RMD_sub_mean);
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        for ii = 1:4
            PV.CV10_PostActivity.(field){roi}(:,ii) = nanmean(PV.CV10_FitTest_RMD_sub_mean.(field){roi}{ii}(:,16:end),2);
        end
        PV.CV10_PostActivity_delta.(field){roi} = PV.CV10_PostActivity.(field){roi}-repmat(PV.CV10_PostActivity.(field){roi}(:,1),1,4);
    end
end
for field_ii = 1:length(fields)
    field = fields{field_ii};
    for roi = 1:16
        temp_matrix = nanmean(PV.CV10_PostActivity_delta.True{roi});
        PV.CV10_PostActivity_deltaNorm.(field){roi} = PV.CV10_PostActivity_delta.(field){roi}./repmat(temp_matrix,size(PV.CV10_PostActivity_delta.(field){roi},1),1);
    end
end

% colors_field = [0,0,0;...
%     127.5,127.5,127.5;...
%     166,73,90;...
%     217,156,161;...
%     60,94,115;...
%     124,158,166;...
%     217,181,150;...
%     166,120,93]/255;
figure; hold on;
subplot(3,1,1); hold on;
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = VIP.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*4)*(roi-1)+[1:4]+(ii_field-1)*4+4*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:4
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end

        for ii = 1:4
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*4)*(roi-1)+[1:4]+(4-1)*4+4*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.001 0.012]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('df/f'); title('VIP postMvmOnset activity');
subplot(3,1,2); hold on;
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = SOM.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*4)*(roi-1)+[1:4]+(ii_field-1)*4+4*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:4
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*4)*(roi-1)+[1:4]+(4-1)*4+4*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.001 0.023]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('df/f'); title('SOM postMvmOnset activity');        
subplot(3,1,3); hold on;
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = PV.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*4)*(roi-1)+[1:4]+(ii_field-1)*4+4*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:4
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*4)*(roi-1)+[1:4]+(4-1)*4+4*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.001 0.018]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('df/f'); title('PV. postMvmOnset activity');        
saveas(gcf,'PostMvmOnsetActivity_FitTrace_FullBetaOnly2_bar.fig','fig');
saveas(gcf,'PostMvmOnsetActivity_FitTrace_FullBetaOnly2_bar.png','png');
print('PostMvmOnsetActivity_FitTrace_FullBetaOnly2_bar.pdf','-dpdf','-bestfit'); pause(1);


figure; hold on;
subplot(3,1,1); hold on;
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = VIP.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*3)*(roi-1)+[1:3]+(ii_field-1)*3+3*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:3
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.006 0.003]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('delta df/f'); title('VIP. delta postMvmOnset activity (- naive)');
subplot(3,1,2); hold on;
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = SOM.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*3)*(roi-1)+[1:3]+(ii_field-1)*3+3*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:3
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.001 0.011]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('delta df/f'); title('SOM. delta postMvmOnset activity (- naive)');        
subplot(3,1,3); hold on;
for roi = 1:16
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = PV.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = (length(fields)*3)*(roi-1)+[1:3]+(ii_field-1)*3+3*(roi-1);
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:3
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
    x_tickpos(roi) = x_control_1(end)+0.5;
end
xlim([-2,x_control(end)+3]); ylim([-0.0025 0.005]);
xticks(x_tickpos);xticklabels(Ordered_ROI);
ylabel('delta df/f'); title('PV. delta postMvmOnset activity (- naive)');        
saveas(gcf,'PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2_bar.fig','fig');
saveas(gcf,'PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2_bar.png','png');
print('PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2_bar.pdf','-dpdf','-bestfit'); pause(1);

% use line plotting
figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = VIP.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:4],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:4
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 4.5]); xticks([]); xticklabels;
    ylim([-0.001 0.013]); yticks([0:0.005:0.01]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'VIP_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'VIP_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.png','png');
print('V_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = SOM.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:4],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:4
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 4.5]); xticks([]); xticklabels;
    ylim([-0.001 0.022]); yticks([0:0.005:0.02]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'SOM_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'SOM_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.png','png');
print('S_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = PV.CV10_PostActivity.(field){reordered_module(roi)}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:4],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:4
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 4.5]); xticks([]); xticklabels;
    ylim([-0.001 0.017]); yticks([0:0.005:0.02]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'PV_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'PV_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.png','png');
print('P_PostMvmOnsetActivity_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = VIP.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:3],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:3
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 3.5]); xticks([]); xticklabels;
    line(xlim,[0 0],'color','k','linestyle',':')
    ylim([-0.006 0.003]); yticks([-0.005:0.005:0.01]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'VIP_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'VIP_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.png','png');
print('V_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = SOM.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:3],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:3
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 3.5]); xticks([]); xticklabels;
    line(xlim,[0 0],'color','k','linestyle',':')
    ylim([-0.001 0.011]); yticks([-0.005:0.005:0.01]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'SOM_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'SOM_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.png','png');
print('S_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

figure; hold on; set(gcf,'color','w','pos',[100 100 1500 100]);
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
for roi = 1:16
    subplot(1,16,roi); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = PV.CV10_PostActivity_delta.(field){reordered_module(roi)}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        plot([1:1:3],temp_1,'color',colors_field(ii_field,:));
        for ii = 1:3
            line([ii ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0.5 3.5]); xticks([]); xticklabels;
    line(xlim,[0 0],'color','k','linestyle',':')
    ylim([-0.0025 0.005]); yticks([-0.005:0.005:0.01]);yticklabels([]);
    set(gca,'XColor','k','YColor','k');
    box on;
end
saveas(gcf,'PV_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.fig','fig');
saveas(gcf,'PV_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.png','png');
print('PV_PostMvmOnsetActivityChange_FitTrace_FullBetaOnly2.pdf','-dpdf','-bestfit'); pause(1);

save(['GLM_' Model '_FitTraceRwdMVM_FullBetaOnly2.mat'],'VIP','SOM','PV','-append');

%% M1, VIP and SOM, 4 bins after MvmOnset
roi = 8;
bins = [16:30;31:45;46:60;61:75];
fields = fieldnames(VIP.CV10_FitTest_RMD_sub_mean);
for bin_ii = 1:size(bins,1)
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        for ii = 1:4
            temp = VIP.CV10_FitTest_RMD_sub_mean.(field){roi}{ii}-VIP.CV10_FitTest_RMD_sub_mean.(field){roi}{1};
            rM1_BinnedActivityChange.VIP.(field){bin_ii}(:,ii) = nanmean(temp(:,bins(bin_ii,:)),2);
            temp = SOM.CV10_FitTest_RMD_sub_mean.(field){roi}{ii}-SOM.CV10_FitTest_RMD_sub_mean.(field){roi}{1};
            rM1_BinnedActivityChange.SOM.(field){bin_ii}(:,ii) = nanmean(temp(:,bins(bin_ii,:)),2);
            temp = SOM.CV10_FitTest_RMD_sub_mean.(field){roi}{ii};
            rM1_BinnedActivity.SOM.(field){bin_ii}(:,ii) = nanmean(temp(:,bins(bin_ii,:)),2);
        end
    end
end

IN = 'VIP';
fields = {'True','Full','Mvm','Lick','Cue','Reward'};
figure; set(gcf,'color','w','pos',[200,200,800,200]); hold on;   
for bin_ii = 1:4
    subplot(1,4,bin_ii); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = rM1_BinnedActivityChange.(IN).(field){bin_ii}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = [1:3]+(ii_field-1)*4;
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:3
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
        x_tickpos(ii_field) = x_control(2);
    end
    xlim([-2,x_control(end)+3]); ylim([-0.007 0.003]);
    xticks(x_tickpos); xticklabels([]);
    ylabel('delta df/f');
%     x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
    if bin_ii == 1
        title('0-0.5 s');
    elseif bin_ii == 2
        title('0.5-1 s')
    elseif bin_ii == 3
        title('1-1.5 s')
    elseif bin_ii == 4
        title('1.5-2 s')
    end
    axis square;
end
saveas(gcf,'rM1_VIP_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'rM1_VIP_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.png','png');
saveas(gcf,'rM1_VIP_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.pdf','pdf');

IN = 'SOM';
figure; set(gcf,'color','w','pos',[200,200,800,200]); hold on;   
for bin_ii = 1:4
    subplot(1,4,bin_ii); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = rM1_BinnedActivityChange.(IN).(field){bin_ii}(:,2:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = [1:3]+(ii_field-1)*4;
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:3
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
        x_tickpos(ii_field) = x_control(2);
    end
    xlim([-2,x_control(end)+3]); ylim([-0.001 0.015]);
    xticks(x_tickpos); xticklabels([]);
    ylabel('delta df/f');
%     x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
%     x_tickpos(roi) = x_control_1(end)+0.5;
    if bin_ii == 1
        title('0-0.5 s');
    elseif bin_ii == 2
        title('0.5-1 s')
    elseif bin_ii == 3
        title('1-1.5 s')
    elseif bin_ii == 4
        title('1.5-2 s')
    end
    axis square;
end
saveas(gcf,'rM1_SOM_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'rM1_SOM_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.png','png');
saveas(gcf,'rM1_SOM_Binned_FitTraceDelta_RwdMVM_FullBetaOnly2.pdf','pdf');

IN = 'SOM';
figure; set(gcf,'color','w','pos',[200,200,800,200]); hold on;   
for bin_ii = 1:4
    subplot(1,4,bin_ii); hold on;
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = rM1_BinnedActivity.(IN).(field){bin_ii}(:,1:end);
        temp_1 = nanmean(temp_var);
        temp_2 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
        x_control = [1:4]+(ii_field-1)*4;
        bar(x_control,temp_1,'Edgecolor','none','Facecolor',colors_field(ii_field,:));
        for ii = 1:4
            line([x_control(ii) x_control(ii)],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',colors_field(ii_field,:));
        end
        x_tickpos(ii_field) = (x_control(2)+x_control(3))/2;
    end
    xlim([-2,x_control(end)+3]); ylim([-0.005 0.026]);
    xticks(x_tickpos); xticklabels([]);
    ylabel('delta df/f');
%     x_control_1 = (length(fields)*3)*(roi-1)+[1:3]+(4-1)*3+3*(roi-1);
%     x_tickpos(roi) = x_control_1(end)+0.5;
    if bin_ii == 1
        title('0-0.5 s');
    elseif bin_ii == 2
        title('0.5-1 s')
    elseif bin_ii == 3
        title('1-1.5 s')
    elseif bin_ii == 4
        title('1.5-2 s')
    end
    axis square;
end
saveas(gcf,'rM1_SOM_Binned_FitTrace_RwdMVM_FullBetaOnly2.fig','fig');
saveas(gcf,'rM1_SOM_Binned_FitTrace_RwdMVM_FullBetaOnly2.png','png');
saveas(gcf,'rM1_SOM_Binned_FitTrace_RwdMVM_FullBetaOnly2.pdf','pdf');

save(['GLM_' Model '_FitTraceRwdMVM_FullBetaOnly.mat'],'rM1_BinnedActivity','rM1_BinnedActivityChange','Model','-append');

%% Plot example predictors for model
load('Z:\People\Chi\WFLP_IN\SOM\CR_3183958-L\EventAligned_Gap500\CR_3183958-L_FrameIndexForGLM.mat')
figure; hold on
ii = 1; jj = 2;
plot(Traces.LeverSpeed{1,1}{1,jj}*5)
plot(Traces.LickRate{1,1}{1,jj}/10-5)
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.MvmOn{1,1}{1,jj}) = 1;
plot(temp_trace-1.2)
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.LickBoutOn{1,1}{1,jj}) = 1;
plot(temp_trace-2.4)
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.Cue{1,1}{1,jj}) = 1;
plot(temp_trace+2)
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.RewardOn{1,1}{1,jj}) = 1;
plot(temp_trace+4)
xlim([4170,4400])
ylim([-6 6])

% Plot real trace for figure 1
figure; hold on
ii = 1; jj = 2;
plot(Traces.LeverPosition{1,ii}{1,jj})
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.Cue{1,ii}{1,jj}) = 1;
plot(temp_trace)
temp_trace = zeros(1,9000);
temp_trace(FrameIndex.RewardOn{1,ii}{1,jj}) = 1;
plot(temp_trace+1.5)
line(xlim,[1.75-2.5132,1.75-2.5132])
line(xlim,[1.5-2.5132,1.5-2.5132])
xlim([4150,4900])
ylim([-2 3])



