% Extract full model from all regression to save loading time
close all
clc

clearvars -except VIP SOM PV Stats ColorMap

Model = 'cv10_Model_7';
cd(['Z:\People\Chi\WFLP_IN\GLM\' Model]);

INs = {'VIP','SOM','PV'};

for curr_IN = 1:length(INs)
    
    clearvars -except VIP SOM PV INs curr_IN Model Stats ColorMap cv_num 
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
        tic
        disp([IN ' ' Animal]);
        data_folder = ['Z:\People\Chi\WFLP_IN\GLM\' Model];
        file_name = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_ridgeMML_CV10.mat'];
        load(file_name,'CV10','shuffle_CV10','predictors_tag','-mat');
        CV10 = CV10.Full;
        shuffle_CV10 = shuffle_CV10.Full;
        file_name_2 = [data_folder filesep IN '_' Initial '_' Animal '_df_f_all_pseudoICA_GLM_stage_ridgeMML_CV10_FullOnly.mat'];
        save(file_name_2,'CV10','shuffle_CV10','predictors_tag','-v7.3');
        clear CV10 shuffle_CV10 predictors_tag
        toc
    end
end


        
        
        
        
        