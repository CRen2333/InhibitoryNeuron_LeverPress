%% GLM
% Fit for each module, session based
% Specifying the number of coefficients that corresponds to the number of data points before/after each behavior event
% Pinto: PN: 0–4 s; RW: 0–2 s; stimulus: 0–0.2 s; prep cue: 0–0.8 s; licking onset: ~0.5–3 s; licking offset: ~0.5–8 s
% Musall: stimuli: start to trial end, mvm events: -0.5-2 s
% use ridgeMML (Musall, NN, 2019) to make it faster

clear all
close all
clc
INs = {'VIP','SOM','PV'};
Model = 'cv10_Model_7';
full_only = false;

cv_num = 10; % 10-fold

for curr_IN = 1:length(INs)
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

    for curr_animal = 1:length(Animals)
        clearvars -except IN INs curr_IN Initial Animals curr_animal Model full_only cv_num

        Animal = Animals{curr_animal};
        General_path = fullfile('Z:\People\Chi\WFLP_IN',IN,[Initial '_' Animal]);

        save_folder = ['Z:\People\Chi\WFLP_IN\GLM\' Model];

        % Shift for binary traces
        % Use analog trace for licking rate and lever position speed, no shift for analog traces
        framerate = 30;
        Shift.CueOn_start = 0*framerate; Shift.CueOn_end = 0*framerate;
        Shift.CueMiddle_start = 0*framerate; Shift.CueMiddle_end = 0.5*framerate;
        Shift.RewardOn_start = 0; Shift.RewardOn_end = 0.5*framerate;
        Shift.PunishOn_start = 0; Shift.PunishOn_end = 0.2*framerate;
        Shift.LickBoutOn_start = 0.5*framerate; Shift.LickBoutOn_end = 2*framerate;
        Shift.LickBoutOff_start = 0.5*framerate; Shift.LickBoutOff_end = 2*framerate;
        Shift.MvmOn_start = 0.5*framerate; Shift.MvmOn_end = 2*framerate;
        Shift.MvmOff_start = 0.5*framerate; Shift.MvmOff_end = 2*framerate;
        Shift.Analog_start =0*framerate; Shift.Analog_end = 0*framerate;

        predictors_tag = [repmat({'CueMiddle'},Shift.CueMiddle_start+Shift.CueMiddle_end+1,1);...
                repmat({'RewardOn'},Shift.RewardOn_start+Shift.RewardOn_end+1,1);...
                repmat({'MvmOn'},Shift.MvmOn_start+Shift.MvmOn_end+1,1);...
                repmat({'Mvm_LeverSpeed'},Shift.Analog_start+Shift.Analog_start+1,1);...
                repmat({'LickBoutOn'},Shift.LickBoutOn_start+Shift.LickBoutOn_end+1,1);...
                repmat({'LickRate'},Shift.Analog_start+Shift.Analog_start+1,1)];

        load([General_path filesep 'EventAligned_Gap500' filesep Initial '_' Animal '_FrameIndexForGLM.mat'],'FrameIndex','Traces');
        df_f_folder = ['Z:\People\Chi\WFLP_IN\WFIN_DFOF_TrainingGroup\' IN filesep Initial '_' Animal];
        Dates = dir(df_f_folder);
        Dates = {Dates.name};
        Dates = Dates(3:end);

        for curr_session = 1:11
            Date = Dates{curr_session};      
            load([df_f_folder filesep Date filesep Initial '_' Date '_' Animal '_df_f_all_pseudoICA.mat'],'df_f_ROI_pseudoICA');

            LickBoutOn_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            LickRate_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            MvmOn_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            MvmLeverSpeed_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));            
            CueMiddle_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            RewardOn_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            PunishOn_Matrix = cell(size(df_f_ROI_pseudoICA{1},2),size(df_f_ROI_pseudoICA{1},1));
            
            for curr_block = 1:length(df_f_ROI_pseudoICA{1})            
                if isempty(FrameIndex.MvmOn{curr_session}{curr_block})
                    for roiNum = 1:size(df_f_ROI_pseudoICA,2) % for each roi
                        df_f_ROI_pseudoICA{roiNum}{curr_block} = [];
                    end
                    continue
                end
                % Lick
                LickBoutOn_Matrix{curr_block,1} = CR_DesignPredictorMatrix(FrameIndex.LickBoutOn{curr_session}{curr_block},9000,Shift.LickBoutOn_start,Shift.LickBoutOn_end);
                LickRate_Matrix{curr_block,1} = Traces.LickRate{curr_session}{curr_block};

                % Mvm
                MvmOn_Matrix{curr_block,1} = CR_DesignPredictorMatrix(FrameIndex.MvmOn{curr_session}{curr_block},9000,Shift.MvmOn_start,Shift.MvmOn_end);
                MvmLeverSpeed_Matrix{curr_block,1} = Traces.LeverSpeed{curr_session}{curr_block};

                % Cue, rewrad, and punishment signal signal
                CueMiddle_Matrix{curr_block,1} = CR_DesignPredictorMatrix(FrameIndex.Cue{curr_session}{curr_block},9000,Shift.CueMiddle_start,Shift.CueMiddle_end);            
                RewardOn_Matrix{curr_block,1} = CR_DesignPredictorMatrix(FrameIndex.Reward{curr_session}{curr_block},9000,Shift.RewardOn_start,Shift.RewardOn_end);
                PunishOn_Matrix{curr_block,1} = CR_DesignPredictorMatrix(FrameIndex.PunishOn{curr_session}{curr_block},9000,Shift.PunishOn_start,Shift.PunishOn_end); 

            end
            LickBoutOn_Matrix = cell2mat(LickBoutOn_Matrix);
            LickRate_Matrix = cell2mat(LickRate_Matrix);
            MvmOn_Matrix = cell2mat(MvmOn_Matrix);
            MvmLeverSpeed_Matrix = cell2mat(MvmLeverSpeed_Matrix);
            CueMiddle_Matrix = cell2mat(CueMiddle_Matrix);
            RewardOn_Matrix = cell2mat(RewardOn_Matrix);
            PunishOn_Matrix = cell2mat(PunishOn_Matrix);

            % GLM predictor, no PunishOn as it will cause fitting failure when all zero
            predictors_all{1,curr_session} = [CueMiddle_Matrix,RewardOn_Matrix,...
                MvmOn_Matrix,MvmLeverSpeed_Matrix, ...
                LickBoutOn_Matrix,LickRate_Matrix]';

            % df/f
            for roiNum = 1:size(df_f_ROI_pseudoICA,2) % for each roi
                df_f_all{roiNum}{1,curr_session} = cell2mat(df_f_ROI_pseudoICA{roiNum});
            end
        end

        clear LickBoutOn_Matrix LickBoutOff_Matrix LickRate_Matrix MvmOn_Matrix MvmOff_Matrix MvmLeverPos_Matrix MvmLeverSpeed_Matrix CueOn_Matrix CueMiddle_Matrix RewardOn_Matrix PunishOn_Matrix
        clear df_f_ROI_pseudoICA

        Stage_sessions = {[1],[2,3,4],[5,6,7,8],[9,10,11]};
        Stages = {'Naive','Early','Middle','Late'};
        ROIs = {'rVisual','raS1BC','laS1BC','lM1','PPC','lVisual','rS1HL','rM1','lS1HL','lpS1BC','aRSC','rpS1BC','pRSC','lM1_S1FL','rM1_S1FL','M2'};

        for stage_ii = 1:length(Stage_sessions)
            
            predictors = cell2mat(predictors_all(Stage_sessions{stage_ii}))';        
            df_f = [];
            for roiNum = 1:size(df_f_all,2) % for each roi
                df_f(:,roiNum) = cell2mat(df_f_all{roiNum}(Stage_sessions{stage_ii}));
            end  
%             check whether the matrix is full rank
%             [~, fullQRR] = qr(bsxfun(@rdivide,predictors,sqrt(sum(predictors.^2))),0); %orthogonalize normalized design matrix
%             figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
%             axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
%             temp = ~(abs(diag(fullQRR)) > max(size(predictors)) * eps(fullQRR(1)));
            
            % Full model
            disp([Animal ' ' Stages{stage_ii} ' Full']);
            tic
            for cv_ii = 1:cv_num
 
                clear training_index
                training_index = true(1,size(predictors,1));
                test_size = length(training_index)/cv_num;
                training_index([(cv_ii-1)*test_size+1:cv_ii*test_size]) = false;
                test_index = ~training_index;

                [lambdas, betas] = ridgeMML(df_f(training_index,:), predictors(training_index,:), false);
                Bias(:,:,cv_ii) = betas(1,:);
                Betas(:,:,cv_ii) = betas(2:end,:); 
                Lambda(:,:,cv_ii) = lambdas;   
                temp_fit_trace = predictors*betas(2:end,:)+betas(1,:);
                [MSE_train(cv_ii,:),MSE_test(cv_ii,:),VarExp_train(cv_ii,:),VarExp_test(cv_ii,:),Corr_train(cv_ii,:),Corr_test(cv_ii,:)] = ...
                    CR_EvaluatePostridgeMML(df_f, temp_fit_trace, training_index);
                FitTrace_training(:,:,cv_ii) = temp_fit_trace(training_index,:);
                FitTrace_test(:,:,cv_ii) = temp_fit_trace(~training_index,:);
                TrueTrace_training(:,:,cv_ii) = df_f(training_index,:);
                TrueTrace_test(:,:,cv_ii) = df_f(~training_index,:);
            end            
            toc
            CV10.Full.MSE_train{stage_ii} = MSE_train;
            CV10.Full.MSE_test{stage_ii} = MSE_test;
            CV10.Full.VarExp_train{stage_ii} = VarExp_train;
            CV10.Full.VarExp_test{stage_ii} = VarExp_test;
            CV10.Full.Corr_train{stage_ii} = Corr_train;
            CV10.Full.Corr_test{stage_ii} = Corr_test;
            CV10.Full.Betas{stage_ii} = Betas;
            CV10.Full.Bias{stage_ii} = Bias;
            CV10.Full.Lambda{stage_ii} = Lambda;
            CV10.Full.FitTrace_training{stage_ii} = FitTrace_training;
            CV10.Full.FitTrace_test{stage_ii} = FitTrace_test;
            CV10.Full.TrueTrace_training{stage_ii} = TrueTrace_training;
            CV10.Full.TrueTrace_test{stage_ii} = TrueTrace_test;

            clear MSE_train MSE_test VarExp_train VarExp_test Corr_train Corr_test Betas Bias Lambda
            clear FitTrace_training FitTrace_test TrueTrace_training TrueTrace_test
            
            if ~full_only % keep one predictor and shuffle the rest one-by-one
                predictor_names = unique(predictors_tag);
                for predictor_ii = 1:length(predictor_names)
                    curr_predictor = predictor_names{predictor_ii};
                    predictor_index = cellfun(@(x) contains(x, curr_predictor), predictors_tag);
                    shuffle_index = randperm(size(predictors,1));
                    predictors_shuffle = predictors;
                    predictors_shuffle(:,~predictor_index) = predictors(shuffle_index,~predictor_index);
                    
                    disp([Animal ' ' Stages{stage_ii} ' ' curr_predictor ' Only']);
                    tic
                    for cv_ii = 1:cv_num 
                        clear training_index
                        training_index = true(1,size(predictors,1));
                        test_size = length(training_index)/cv_num;
                        training_index([(cv_ii-1)*test_size+1:cv_ii*test_size]) = false;
                        test_index = ~training_index;

                        [lambdas, betas] = ridgeMML(df_f(training_index,:), predictors_shuffle(training_index,:), false);
                        Bias(:,:,cv_ii) = betas(1,:);
                        Betas(:,:,cv_ii) = betas(2:end,:); 
                        Lambda(:,:,cv_ii) = lambdas;   
                        temp_fit_trace = predictors_shuffle*betas(2:end,:)+betas(1,:);
                        [MSE_train(cv_ii,:),MSE_test(cv_ii,:),VarExp_train(cv_ii,:),VarExp_test(cv_ii,:),Corr_train(cv_ii,:),Corr_test(cv_ii,:)] = ...
                            CR_EvaluatePostridgeMML(df_f, temp_fit_trace, training_index);
                        FitTrace_training(:,:,cv_ii) = temp_fit_trace(training_index,:);
                        FitTrace_test(:,:,cv_ii) = temp_fit_trace(~training_index,:);
                        TrueTrace_training(:,:,cv_ii) = df_f(training_index,:);
                        TrueTrace_test(:,:,cv_ii) = df_f(~training_index,:);
                    end            
                    toc
                    CV10.(curr_predictor).MSE_train{stage_ii} = MSE_train;
                    CV10.(curr_predictor).MSE_test{stage_ii} = MSE_test;
                    CV10.(curr_predictor).VarExp_train{stage_ii} = VarExp_train;
                    CV10.(curr_predictor).VarExp_test{stage_ii} = VarExp_test;
                    CV10.(curr_predictor).Corr_train{stage_ii} = Corr_train;
                    CV10.(curr_predictor).Corr_test{stage_ii} = Corr_test;
                    CV10.(curr_predictor).Betas{stage_ii} = Betas;
                    CV10.(curr_predictor).Bias{stage_ii} = Bias;
                    CV10.(curr_predictor).Lambda{stage_ii} = Lambda;
                    CV10.(curr_predictor).FitTrace_training{stage_ii} = FitTrace_training;
                    CV10.(curr_predictor).FitTrace_test{stage_ii} = FitTrace_test;
                    CV10.(curr_predictor).TrueTrace_training{stage_ii} = TrueTrace_training;
                    CV10.(curr_predictor).TrueTrace_test{stage_ii} = TrueTrace_test;

                    clear MSE_train MSE_test VarExp_train VarExp_test Corr_train Corr_test Betas Bias Lambda
                    clear FitTrace_training FitTrace_test TrueTrace_training TrueTrace_test predictors_shuffle
                end
            end
            
            % ********** SHUFFLE **********
            % Full model
            disp([Animal ' ' Stages{stage_ii} ' Full SHUFFLE']);
            shuffle_index = randperm(size(predictors,1));
            predictors_shuffle = predictors(shuffle_index,:);
            tic
            for cv_ii = 1:cv_num
 
                clear training_index
                training_index = true(1,size(predictors,1));
                test_size = length(training_index)/cv_num;
                training_index([(cv_ii-1)*test_size+1:cv_ii*test_size]) = false;
                test_index = ~training_index;

                [lambdas, betas] = ridgeMML(df_f(training_index,:), predictors_shuffle(training_index,:), false);
                Bias(:,:,cv_ii) = betas(1,:);
                Betas(:,:,cv_ii) = betas(2:end,:); 
                Lambda(:,:,cv_ii) = lambdas;   
                temp_fit_trace = predictors_shuffle*betas(2:end,:)+betas(1,:);
                [MSE_train(cv_ii,:),MSE_test(cv_ii,:),VarExp_train(cv_ii,:),VarExp_test(cv_ii,:),Corr_train(cv_ii,:),Corr_test(cv_ii,:)] = ...
                    CR_EvaluatePostridgeMML(df_f, temp_fit_trace, training_index);
                FitTrace_training(:,:,cv_ii) = temp_fit_trace(training_index,:);
                FitTrace_test(:,:,cv_ii) = temp_fit_trace(~training_index,:);
                TrueTrace_training(:,:,cv_ii) = df_f(training_index,:);
                TrueTrace_test(:,:,cv_ii) = df_f(~training_index,:);
            end            
            toc
            shuffle_CV10.Full.MSE_train{stage_ii} = MSE_train;
            shuffle_CV10.Full.MSE_test{stage_ii} = MSE_test;
            shuffle_CV10.Full.VarExp_train{stage_ii} = VarExp_train;
            shuffle_CV10.Full.VarExp_test{stage_ii} = VarExp_test;
            shuffle_CV10.Full.Corr_train{stage_ii} = Corr_train;
            shuffle_CV10.Full.Corr_test{stage_ii} = Corr_test;
            shuffle_CV10.Full.Betas{stage_ii} = Betas;
            shuffle_CV10.Full.Bias{stage_ii} = Bias;
            shuffle_CV10.Full.Lambda{stage_ii} = Lambda;
            shuffle_CV10.Full.FitTrace_training{stage_ii} = FitTrace_training;
            shuffle_CV10.Full.FitTrace_test{stage_ii} = FitTrace_test;
            shuffle_CV10.Full.TrueTrace_training{stage_ii} = TrueTrace_training;
            shuffle_CV10.Full.TrueTrace_test{stage_ii} = TrueTrace_test;

            clear MSE_train MSE_test VarExp_train VarExp_test Corr_train Corr_test Betas Bias Lambda
            clear FitTrace_training FitTrace_test TrueTrace_training TrueTrace_test predictors_shuffle
            
            if ~full_only % shuffle predictors one-by-one
                predictor_names = unique(predictors_tag);
                for predictor_ii = 1:length(predictor_names)
                    curr_predictor = predictor_names{predictor_ii};
                    predictor_index = cellfun(@(x) contains(x, curr_predictor), predictors_tag);
                    shuffle_index = randperm(size(predictors,1));
                    predictors_shuffle = predictors;
                    predictors_shuffle(:,predictor_index) = predictors(shuffle_index,predictor_index);
                    
                    disp([Animal ' ' Stages{stage_ii} ' ' curr_predictor ' SHUFFLE']);
                    tic
                    for cv_ii = 1:cv_num
 
                        clear training_index
                        training_index = true(1,size(predictors,1));
                        test_size = length(training_index)/cv_num;
                        training_index([(cv_ii-1)*test_size+1:cv_ii*test_size]) = false;
                        test_index = ~training_index;

                        [lambdas, betas] = ridgeMML(df_f(training_index,:), predictors_shuffle(training_index,:), false);
                        Bias(:,:,cv_ii) = betas(1,:);
                        Betas(:,:,cv_ii) = betas(2:end,:); 
                        Lambda(:,:,cv_ii) = lambdas;   
                        temp_fit_trace = predictors_shuffle*betas(2:end,:)+betas(1,:);
                        [MSE_train(cv_ii,:),MSE_test(cv_ii,:),VarExp_train(cv_ii,:),VarExp_test(cv_ii,:),Corr_train(cv_ii,:),Corr_test(cv_ii,:)] = ...
                            CR_EvaluatePostridgeMML(df_f, temp_fit_trace, training_index);
                        FitTrace_training(:,:,cv_ii) = temp_fit_trace(training_index,:);
                        FitTrace_test(:,:,cv_ii) = temp_fit_trace(~training_index,:);
                        TrueTrace_training(:,:,cv_ii) = df_f(training_index,:);
                        TrueTrace_test(:,:,cv_ii) = df_f(~training_index,:);
                    end            
                    toc
                    shuffle_CV10.(curr_predictor).MSE_train{stage_ii} = MSE_train;
                    shuffle_CV10.(curr_predictor).MSE_test{stage_ii} = MSE_test;
                    shuffle_CV10.(curr_predictor).VarExp_train{stage_ii} = VarExp_train;
                    shuffle_CV10.(curr_predictor).VarExp_test{stage_ii} = VarExp_test;
                    shuffle_CV10.(curr_predictor).Corr_train{stage_ii} = Corr_train;
                    shuffle_CV10.(curr_predictor).Corr_test{stage_ii} = Corr_test;
                    shuffle_CV10.(curr_predictor).Betas{stage_ii} = Betas;
                    shuffle_CV10.(curr_predictor).Bias{stage_ii} = Bias;
                    shuffle_CV10.(curr_predictor).Lambda{stage_ii} = Lambda;
                    shuffle_CV10.(curr_predictor).FitTrace_training{stage_ii} = FitTrace_training;
                    shuffle_CV10.(curr_predictor).FitTrace_test{stage_ii} = FitTrace_test;
                    shuffle_CV10.(curr_predictor).TrueTrace_training{stage_ii} = TrueTrace_training;
                    shuffle_CV10.(curr_predictor).TrueTrace_test{stage_ii} = TrueTrace_test;

                    clear MSE_train MSE_test VarExp_train VarExp_test Corr_train Corr_test Betas Bias Lambda
                    clear FitTrace_training FitTrace_test TrueTrace_training TrueTrace_test predictors_shuffle
                end
            end
            
        end
        if ~exist(save_folder)
            mkdir(save_folder)
        end
        file_name = [save_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_ridgeMML_CV10.mat'];
        save(file_name,'CV10','shuffle_CV10','predictors_all',...
            'predictors_tag','df_f_all','Shift','ROIs','Stage_sessions','Stages','-v7.3');

        clear CV10 shuffle_CV10 shuffle_index predictors_all df_f_all predictors predictors_tag

    end       
end

%% Check corss corr, beta and bias
close all
clc

clearvars -except VIP SOM PV Stats ColorMap

Model = 'cv10_Model_7';
cd(['Z:\People\Chi\WFLP_IN\GLM\' Model]);

INs = {'VIP','SOM','PV'};

tic
for curr_IN = 1:length(INs)
    
    clearvars -except VIP SOM PV INs curr_IN Model Stats ColorMap beta_exe 
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
        Animal = Animals{curr_animal};
        disp([IN ' ' Animal]);
        data_folder = ['Z:\People\Chi\WFLP_IN\GLM\' Model];
        file_name = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_ridgeMML_CV10.mat'];
        load(file_name,'CV10','shuffle_CV10','predictors_tag','-mat');
        
        fields = fieldnames(CV10);
        for field_ii = 1:length(fields)
            field = fields{field_ii};
            for ii = 1:4
                for roi = 1:16
                    CV10_all.(field).Betas{roi}{ii}(curr_animal,:) = nanmean(CV10.(field).Betas{ii}(:,reordered_module(roi),:),3);
                    CV10_all.(field).Bias{roi}(curr_animal,ii) = nanmean(CV10.(field).Bias{ii}(:,reordered_module(roi),:),3);
                    CV10_all.(field).Corr_train{roi}(curr_animal,ii) = nanmean(CV10.(field).Corr_train{ii}(:,reordered_module(roi)),1);
                    CV10_all.(field).Corr_test{roi}(curr_animal,ii) = nanmean(CV10.(field).Corr_test{ii}(:,reordered_module(roi)),1);
                end
            end
        end
        
        fields = fieldnames(shuffle_CV10);
        for field_ii = 1:length(fields)
            field = fields{field_ii};
            for ii = 1:4
                for roi = 1:16
                    shuffle_CV10_all.(field).Betas{roi}{ii}(curr_animal,:) = nanmean(shuffle_CV10.(field).Betas{ii}(:,reordered_module(roi),:),3);
                    shuffle_CV10_all.(field).Bias{roi}(curr_animal,ii) = nanmean(shuffle_CV10.(field).Bias{ii}(:,reordered_module(roi),:),3);
                    shuffle_CV10_all.(field).Corr_train{roi}(curr_animal,ii) = nanmean(shuffle_CV10.(field).Corr_train{ii}(:,reordered_module(roi)),1);
                    shuffle_CV10_all.(field).Corr_test{roi}(curr_animal,ii) = nanmean(shuffle_CV10.(field).Corr_test{ii}(:,reordered_module(roi)),1);
                end
            end
        end       
        
        clear CV10 shuffle_CV10
    end
    
    fields = fieldnames(shuffle_CV10_all);
    for field_ii = 1:length(fields)
        field = fields{field_ii};
        for roi = 1:length(ROI)
            Unique.(field).Corr_test{roi} = CV10_all.Full.Corr_test{roi}-shuffle_CV10_all.(field).Corr_test{roi};
        end
    end
    
    switch IN
        case 'VIP'
            VIP.CV10_all = CV10_all;
            VIP.shuffle_CV10_all = shuffle_CV10_all;
            VIP.Unique = Unique;
            VIP.Animals = Animals;
        case 'SOM'
            SOM.CV10_all = CV10_all;
            SOM.shuffle_CV10_all = shuffle_CV10_all;
            SOM.Unique = Unique;
            SOM.Animals = Animals;
        case 'PV'
            PV.CV10_all = CV10_all;
            PV.shuffle_CV10_all = shuffle_CV10_all;
            PV.Unique = Unique;
            PV.Animals = Animals;
    end
    clear CV10_all shuffle_CV10_all Unique
    
    % Stats
    switch IN
        case 'VIP'
            temp_var_all = VIP.CV10_all;
            temp_var_all_2 = VIP.shuffle_CV10_all;
            temp_var_all_3 = VIP.Unique;
        case 'SOM'
            temp_var_all = SOM.CV10_all;
            temp_var_all_2 = SOM.shuffle_CV10_all;
            temp_var_all_3 = SOM.Unique;
        case 'PV'
            temp_var_all = PV.CV10_all;
            temp_var_all_2 = PV.shuffle_CV10_all;
            temp_var_all_3 = PV.Unique;
    end
    
    % Corr
    for roi = 1:16
        fields = fieldnames(temp_var_all);
        for ii = 1:length(fields)
            field = fields{ii};
            % stats
            temp_var = temp_var_all.(field).Corr_train{roi};
            Sessions = repmat([1:4],size(temp_var,1),1);
            Sessions = Sessions(:);
            Animals_test = [];
            for curr_animal = 1:length(Animals)
                Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
            end
            Animals_test = Animals_test(:);
            y = temp_var(:);
            tbl = table(Animals_test,Sessions,y);
            tbl.Animals_test = nominal(tbl.Animals_test);
            lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
            [beta,betanames,stats] = fixedEffects(lme);
            pValue_corr_test.(field)(:,roi) = stats.pValue(2);
        end
    end
    for ii = 1:length(fields)
        field = fields{ii};
        [CV10_FDR_pValue_corr_train.(field)] = mafdr(pValue_corr_test.(field),'BHFDR', true);
    end    
    clear pValue_corr_test
    
    for roi = 1:16
        fields = fieldnames(temp_var_all);
        for ii = 1:length(fields)
            field = fields{ii};
            % stats
            temp_var = temp_var_all.(field).Corr_test{roi};
            Sessions = repmat([1:4],size(temp_var,1),1);
            Sessions = Sessions(:);
            Animals_test = [];
            for curr_animal = 1:length(Animals)
                Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
            end
            Animals_test = Animals_test(:);
            y = temp_var(:);
            tbl = table(Animals_test,Sessions,y);
            tbl.Animals_test = nominal(tbl.Animals_test);
            lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
            [beta,betanames,stats] = fixedEffects(lme);
            pValue_corr_test.(field)(:,roi) = stats.pValue(2);
        end
    end
    for ii = 1:length(fields)
        field = fields{ii};
        [CV10_FDR_pValue_corr_test.(field)] = mafdr(pValue_corr_test.(field),'BHFDR', true);
    end    
    clear pValue_corr_test
    
    % Unique contribution
    for roi = 1:16
        fields = fieldnames(temp_var_all_3);
        for ii = 1:length(fields)
            field = fields{ii};
            % stats
            temp_var = temp_var_all_3.(field).Corr_test{roi};
            
            Sessions = repmat([1:4],size(temp_var,1),1);
            Sessions = Sessions(:);
            Animals_test = [];
            for curr_animal = 1:length(Animals)
                Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
            end
            Animals_test = Animals_test(:);
            y = temp_var(:);
            tbl = table(Animals_test,Sessions,y);
            tbl.Animals_test = nominal(tbl.Animals_test);
            lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
            [beta,betanames,stats] = fixedEffects(lme);
            pValue_corr_test.(field)(:,roi) = stats.pValue(2);
        end
    end
    for ii = 1:length(fields)
        field = fields{ii};
        [unique_CV10_FDR_pValue_corr_test.(field)] = mafdr(pValue_corr_test.(field),'BHFDR', true);
    end
    clear pValue_corr_test
        
    switch IN
        case 'VIP'
            Stats.VIP.CV10_FDR_pValue_corr_train = CV10_FDR_pValue_corr_train;
            Stats.VIP.CV10_FDR_pValue_corr_test = CV10_FDR_pValue_corr_test;
            Stats.VIP.unique_CV10_FDR_pValue_corr_test = unique_CV10_FDR_pValue_corr_test;
        case 'SOM'
            Stats.SOM.CV10_FDR_pValue_corr_train = CV10_FDR_pValue_corr_train;
            Stats.SOM.CV10_FDR_pValue_corr_test = CV10_FDR_pValue_corr_test;
            Stats.SOM.unique_CV10_FDR_pValue_corr_test = unique_CV10_FDR_pValue_corr_test;
        case 'PV'
            Stats.PV.CV10_FDR_pValue_corr_train = CV10_FDR_pValue_corr_train;
            Stats.PV.CV10_FDR_pValue_corr_test = CV10_FDR_pValue_corr_test;
            Stats.PV.unique_CV10_FDR_pValue_corr_test = unique_CV10_FDR_pValue_corr_test;
    end
end
toc

warning off
[CX] = cbrewer('div','PRGn',10);
ColorMap.VIP.color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
ColorMap.VIP.color_line = nanmean(ColorMap.VIP.color_value([2,3],:),1);
[CX] = cbrewer('div','PRGn',10);
ColorMap.SOM.color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];  
ColorMap.SOM.color_line = nanmean(ColorMap.SOM.color_value([2,3],:),1);
[CX] = cbrewer('div','PuOr',12);
ColorMap.PV.color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)]; 
ColorMap.PV.color_line = nanmean(ColorMap.PV.color_value([2,3],:),1);
warning on

save(['GLM_' Model '.mat'],'ColorMap','Stats','VIP','SOM','PV','ROI','Ordered_ROI','reordered_module','-v7.3');

% plot corr
figure; set(gcf,'color','w','position',[100 100 1000 800]); hold on;
fields = fieldnames(VIP.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder
colors_field = [0,0,0;...
    166,73,90;...
    217,156,161;...
    60,94,115;...
    124,158,166;...
    217,181,150;...
    166,120,93]/255;
for roi = 1:16
    subplot(4,4,roi); hold on;        
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = VIP.CV10_all.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.VIP.CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = -VIP.Unique.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.VIP.unique_CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0 x_control(end)+1]);xticks([]);
    title(['VIP ' Ordered_ROI{roi}],'fontsize',8);
    ylim([-0.45 0.45]);yticks([-0.4:0.2:0.4]);
    yticklabels({'0.4','0.2','0','0.2','0.4'});
    ylabel('R^2');
end
saveas(gcf,['VIP_CVRsquare_ridgeMML.fig'],'fig'); pause(1);
saveas(gcf,['VIP_CVRsquare_ridgeMML.png'],'png');

figure; set(gcf,'color','w','position',[100 100 1000 800]); hold on;
fields = fieldnames(SOM.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder
colors_field = [0,0,0;...
    166,73,90;...
    217,156,161;...
    60,94,115;...
    124,158,166;...
    217,181,150;...
    166,120,93]/255;
for roi = 1:16
    subplot(4,4,roi); hold on;        
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = SOM.CV10_all.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.SOM.CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = -SOM.Unique.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.SOM.unique_CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0 x_control(end)+1]);xticks([]);
    title(['SOM ' Ordered_ROI{roi}],'fontsize',8);
    ylim([-0.45 0.45]);yticks([-0.4:0.2:0.4]);
    yticklabels({'0.4','0.2','0','0.2','0.4'});
    ylabel('R^2');
end
saveas(gcf,['SOM_CVRsquare_ridgeMML.fig'],'fig'); pause(1);
saveas(gcf,['SOM_CVRsquare_ridgeMML.png'],'png');

figure; set(gcf,'color','w','position',[100 100 1000 800]); hold on;
fields = fieldnames(PV.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder
colors_field = [0,0,0;...
    166,73,90;...
    217,156,161;...
    60,94,115;...
    124,158,166;...
    217,181,150;...
    166,120,93]/255;
for roi = 1:16
    subplot(4,4,roi); hold on;        
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = PV.CV10_all.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.PV.CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    for ii_field = 1:length(fields)
        field = fields{ii_field};
        temp_var = -PV.Unique.(field).Corr_test{roi};   
        temp_1 = nanmean(temp_var,1);
        temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
        x_control = [1:4]+(ii_field-1)*5;
        if Stats.PV.unique_CV10_FDR_pValue_corr_test.(field)(roi) < 0.05
            bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
        else
            bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
        end    
        for jj = 1:4
            line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
        end
    end
    xlim([0 x_control(end)+1]);xticks([]);
    title(['PV ' Ordered_ROI{roi}],'fontsize',8);
    ylim([-0.45 0.45]);yticks([-0.4:0.2:0.4]);
    yticklabels({'0.4','0.2','0','0.2','0.4'});
    ylabel('R^2');
end
saveas(gcf,['PV_CVRsquare_ridgeMML.fig'],'fig'); pause(1);
saveas(gcf,['PV_CVRsquare_ridgeMML.png'],'png');
close all

% average across cortical regions
fields = fieldnames(VIP.CV10_all);
for ii_field = 1:length(fields)
    field = fields{ii_field};
    for ii = 1:4
        for roi = 1:16
            VIP.CV10_StageOrg.(field).Individual(ii,roi,:) = VIP.CV10_all.(field).Corr_test{roi}(:,ii);
            VIP.Unique_StageOrg.(field).Individual(ii,roi,:) = VIP.Unique.(field).Corr_test{roi}(:,ii);
        end
    end
    VIP.CV10_StageOrg.(field).CorticalMean = squeeze(nanmean(VIP.CV10_StageOrg.(field).Individual,2))';
    VIP.Unique_StageOrg.(field).CorticalMean = squeeze(nanmean(VIP.Unique_StageOrg.(field).Individual,2))';
    
    % stats
    temp_var = VIP.CV10_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.VIP.CorticalMean.CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);
    
    temp_var = VIP.Unique_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.VIP.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);    
end

fields = fieldnames(SOM.CV10_all);
for ii_field = 1:length(fields)
    field = fields{ii_field};
    for ii = 1:4
        for roi = 1:16
            SOM.CV10_StageOrg.(field).Individual(ii,roi,:) = SOM.CV10_all.(field).Corr_test{roi}(:,ii);
            SOM.Unique_StageOrg.(field).Individual(ii,roi,:) = SOM.Unique.(field).Corr_test{roi}(:,ii);
        end
    end
    SOM.CV10_StageOrg.(field).CorticalMean = squeeze(nanmean(SOM.CV10_StageOrg.(field).Individual,2))';
    SOM.Unique_StageOrg.(field).CorticalMean = squeeze(nanmean(SOM.Unique_StageOrg.(field).Individual,2))';
    
    % stats
    temp_var = SOM.CV10_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.SOM.CorticalMean.CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);
    
    temp_var = SOM.Unique_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.SOM.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);   
end

fields = fieldnames(PV.CV10_all);
for ii_field = 1:length(fields)
    field = fields{ii_field};
    for ii = 1:4
        for roi = 1:16
            PV.CV10_StageOrg.(field).Individual(ii,roi,:) = PV.CV10_all.(field).Corr_test{roi}(:,ii);
            PV.Unique_StageOrg.(field).Individual(ii,roi,:) = PV.Unique.(field).Corr_test{roi}(:,ii);
        end
    end
    PV.CV10_StageOrg.(field).CorticalMean = squeeze(nanmean(PV.CV10_StageOrg.(field).Individual,2))';
    PV.Unique_StageOrg.(field).CorticalMean = squeeze(nanmean(PV.Unique_StageOrg.(field).Individual,2))';
    
    % stats
    temp_var = PV.CV10_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.PV.CorticalMean.CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);
    
    temp_var = PV.Unique_StageOrg.(field).CorticalMean;
    Sessions = repmat([1:4],size(temp_var,1),1);
    Sessions = Sessions(:);
    Animals_test = [];
    for curr_animal = 1:size(temp_var,1)
        Animals_test = [Animals_test;ones(1,size(temp_var,2))*curr_animal];
    end
    Animals_test = Animals_test(:);
    y = temp_var(:);
    tbl = table(Animals_test,Sessions,y);
    tbl.Animals_test = nominal(tbl.Animals_test);
    lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
    [beta,betanames,stats] = fixedEffects(lme);
    Stats.PV.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) = stats.pValue(2);   
end

colors_field = [0,0,0;...
    166,73,90;...
    217,156,161;...
    60,94,115;...
    124,158,166;...
    217,181,150;...
    166,120,93]/255;
figure; set(gcf,'color','w','position',[100 100 800 300]); hold on;
subplot(1,3,1); hold on;
fields = fieldnames(VIP.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder   
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = VIP.CV10_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.VIP.CorticalMean.CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = -VIP.Unique_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.VIP.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
xlim([0 x_control(end)+1]);xticks([]);
title('VIP cortical mean','fontsize',8);
ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);
yticklabels({'0.4','0.2','0','0.2','0.4'});
ylabel('R^2');

subplot(1,3,2); hold on;
fields = fieldnames(SOM.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder   
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = SOM.CV10_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.SOM.CorticalMean.CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = -SOM.Unique_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.SOM.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
xlim([0 x_control(end)+1]);xticks([]);
title('SOM cortical mean','fontsize',8);
ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);
yticklabels({'0.4','0.2','0','0.2','0.4'});
ylabel('R^2');

subplot(1,3,3); hold on;
fields = fieldnames(PV.CV10_all);
fields = fields([1,5,6,3,4,2,7]); % reorder   
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = PV.CV10_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.PV.CorticalMean.CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
for ii_field = 1:length(fields)
    field = fields{ii_field};
    temp_var = -PV.Unique_StageOrg.(field).CorticalMean;
    temp_1 = nanmean(temp_var,1);
    temp_2 = nanstd(temp_var,[],1)./sqrt(sum(~isnan(temp_var),1));
    x_control = [1:4]+(ii_field-1)*5;
    if Stats.PV.CorticalMean.unique_CV10_FDR_pValue_corr_test.(field) < 0.05
        bar(x_control,temp_1,'barwidth',1,'edgecolor','w','facecolor',colors_field(ii_field,:));
    else
        bar(x_control,temp_1,'edgecolor',colors_field(ii_field,:),'facecolor','none','barwidth',1);
    end    
    for jj = 1:4
        line([x_control(jj) x_control(jj)],[temp_1(jj)-temp_2(jj),temp_1(jj)+temp_2(jj)],'color',colors_field(ii_field,:));
    end
end
xlim([0 x_control(end)+1]);xticks([]);
title('PV cortical mean','fontsize',8);
ylim([-0.4 0.4]);yticks([-0.4:0.2:0.4]);
yticklabels({'0.4','0.2','0','0.2','0.4'});
ylabel('R^2');

saveas(gcf,['CorticalMean_CVRsquare_ridgeMML.fig'],'fig'); pause(1);
saveas(gcf,['CorticalMean_CVRsquare_ridgeMML.png'],'png');
close all