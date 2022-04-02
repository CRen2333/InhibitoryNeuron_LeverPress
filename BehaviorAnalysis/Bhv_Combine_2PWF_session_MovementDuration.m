%% 
clearvars -except IN PV SOM VIP SOM_2P VIP_2P

% Initial = 'CR';
% Animals = {'3161016-O','3161018-R','3233232-L','3233232-O','3233232-R','3233233-O','3233233-L','3233233-R','3491479-L','3491479-LR','3547207-LR'};
% IN = 'PV';
% 
% Initial = 'CR';
% Animals = {'3183958-L','3183958-R','3183959-LL','3218181-L','3218181-O','3218181-R','3438521-O','3438521-L','3438521-R','3438521-LR','3453262-O','3453262-L'};
% IN = 'SOM';

Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};
IN = 'VIP';

for ii = 1:length(Animals)
    tic
    Animal = Animals{ii};
    disp(Animal)
    cd(['Z:\People\Chi\WFLP_IN' filesep IN filesep Initial '_' Animal filesep 'ForBhv' filesep 'MovAnalysis']);
    load([Initial '_' Animal '_CuedMov.mat'],'CuedMov_SingleAnimal');
    
    for jj = 1:length(CuedMov_SingleAnimal) %usefuldays
        if ~isempty(CuedMov_SingleAnimal{jj}) && ~isempty(CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All)
            temp = CuedMov_SingleAnimal{jj}.Cued_MovOnset_Info_All;
            index_reward = temp(:,5)==1;
            CueDuration{ii}{jj,1} = temp(index_reward,7)-temp(index_reward,6);
            MovDuration{ii}{jj,1} = temp(index_reward,3);
        else
            CueDuration{ii}{jj,1} = [];
            MovDuration{ii}{jj,1} = [];
        end
    end
    clear CuedMov_SingleAnimal     
end
switch IN
    case 'PV'        
        PV.Animals = Animals;
        PV.CueDuration = CueDuration;
        PV.MovDuration = MovDuration;
    case 'SOM'        
        SOM.Animals = Animals;
        SOM.CueDuration = CueDuration;
        SOM.MovDuration = MovDuration;
    case 'VIP'        
        VIP.Animals = Animals;
        VIP.CueDuration = CueDuration;
        VIP.MovDuration = MovDuration;
end
clear CueDuration MovDuration

% 2P
clearvars -except IN PV SOM VIP SOM_2P VIP_2P

% IN = 'VIP';
% Animals = {'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','WL_3526641-R','WL_3526642-L','WL_3526642-R'};
    
IN = 'SOM';
Animals = {'KP_3459921_1','KP_3461990_1','WL_3526578-O','WL_3526580-O','WL_3547273-LR','WL_3547273-R','CR_3786142-L','CR_3887041-L','CR_3887041-R','CR_3936483-O'};
    
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal '\MovAnalysis\' Animal '_CuedMov_All.mat'],'CuedMov_SingleAnimal');
    
    for curr_day = 1:length(CuedMov_SingleAnimal)
        if isempty(CuedMov_SingleAnimal{curr_day})
            CueDuration{curr_animal}{curr_day,1} = [];
            MovDuration{curr_animal}{curr_day,1} = [];
            continue
        end
        
        temp_info = CuedMov_SingleAnimal{curr_day}.Cued_MovOnset_Info_All;
        rewared_mvm_index = temp_info(:,5)==1;
        CueDuration{curr_animal}{curr_day,1} = temp_info(rewared_mvm_index,7)-temp_info(rewared_mvm_index,6);
        MovDuration{curr_animal}{curr_day,1} = temp_info(rewared_mvm_index,3);
    end
    clear CuedMov_SingleAnimal     
end
switch IN
    case 'SOM'        
        SOM_2P.Animals = Animals;
        SOM_2P.CueDuration = CueDuration;
        SOM_2P.MovDuration = MovDuration;
    case 'VIP'        
        VIP_2P.Animals = Animals;
        VIP_2P.CueDuration = CueDuration;
        VIP_2P.MovDuration = MovDuration;
end
clear CueDuration MovDuration
cd(['Z:\People\Chi\TwoP_IN\BhaviorAnalysis\Combine_2PWF']);
save('Combine_2PWF_ForBhv_CueMvmDur.mat','VIP','SOM','VIP_2P','SOM_2P','-v7.3');

Cue_all = cell(1,4);
RwdMvm_all = cell(1,4);
for ii = 1:length(SOM.CueDuration)
    Cue_all{1} = [Cue_all{1};cell2mat(SOM.CueDuration{ii}([1:2]))];
    Cue_all{2} = [Cue_all{2};cell2mat(SOM.CueDuration{ii}([3:8]))];
    Cue_all{3} = [Cue_all{3};cell2mat(SOM.CueDuration{ii}([9:15]))];
    Cue_all{4} = [Cue_all{4};cell2mat(SOM.CueDuration{ii}([16:end]))];
    RwdMvm_all{1} = [RwdMvm_all{1};cell2mat(SOM.MovDuration{ii}([1:2]))];
    RwdMvm_all{2} = [RwdMvm_all{2};cell2mat(SOM.MovDuration{ii}([3:8]))];
    RwdMvm_all{3} = [RwdMvm_all{3};cell2mat(SOM.MovDuration{ii}([9:15]))];
    RwdMvm_all{4} = [RwdMvm_all{4};cell2mat(SOM.MovDuration{ii}([16:end]))];
end
for ii = 1:length(VIP.CueDuration)
    Cue_all{1} = [Cue_all{1};cell2mat(VIP.CueDuration{ii}([1:2]))];
    Cue_all{2} = [Cue_all{2};cell2mat(VIP.CueDuration{ii}([3:8]))];
    Cue_all{3} = [Cue_all{3};cell2mat(VIP.CueDuration{ii}([9:15]))];
    Cue_all{4} = [Cue_all{4};cell2mat(VIP.CueDuration{ii}([16:end]))];
    RwdMvm_all{1} = [RwdMvm_all{1};cell2mat(VIP.MovDuration{ii}([1:2]))];
    RwdMvm_all{2} = [RwdMvm_all{2};cell2mat(VIP.MovDuration{ii}([3:8]))];
    RwdMvm_all{3} = [RwdMvm_all{3};cell2mat(VIP.MovDuration{ii}([9:15]))];
    RwdMvm_all{4} = [RwdMvm_all{4};cell2mat(VIP.MovDuration{ii}([16:end]))];
end   
for ii = 1:length(SOM_2P.CueDuration)
    Cue_all{1} = [Cue_all{1};cell2mat(SOM_2P.CueDuration{ii}([1:2]))];
    Cue_all{2} = [Cue_all{2};cell2mat(SOM_2P.CueDuration{ii}([3:8]))];
    Cue_all{3} = [Cue_all{3};cell2mat(SOM_2P.CueDuration{ii}([9:15]))];
    Cue_all{4} = [Cue_all{4};cell2mat(SOM_2P.CueDuration{ii}([16:end]))];
    RwdMvm_all{1} = [RwdMvm_all{1};cell2mat(SOM_2P.MovDuration{ii}([1:2]))];
    RwdMvm_all{2} = [RwdMvm_all{2};cell2mat(SOM_2P.MovDuration{ii}([3:8]))];
    RwdMvm_all{3} = [RwdMvm_all{3};cell2mat(SOM_2P.MovDuration{ii}([9:15]))];
    RwdMvm_all{4} = [RwdMvm_all{4};cell2mat(SOM_2P.MovDuration{ii}([16:end]))];
end
for ii = 1:length(VIP_2P.CueDuration)
    Cue_all{1} = [Cue_all{1};cell2mat(VIP_2P.CueDuration{ii}([1:2]))];
    Cue_all{2} = [Cue_all{2};cell2mat(VIP_2P.CueDuration{ii}([3:8]))];
    Cue_all{3} = [Cue_all{3};cell2mat(VIP_2P.CueDuration{ii}([9:15]))];
    Cue_all{4} = [Cue_all{4};cell2mat(VIP_2P.CueDuration{ii}([16:end]))];
    RwdMvm_all{1} = [RwdMvm_all{1};cell2mat(VIP_2P.MovDuration{ii}([1:2]))];
    RwdMvm_all{2} = [RwdMvm_all{2};cell2mat(VIP_2P.MovDuration{ii}([3:8]))];
    RwdMvm_all{3} = [RwdMvm_all{3};cell2mat(VIP_2P.MovDuration{ii}([9:15]))];
    RwdMvm_all{4} = [RwdMvm_all{4};cell2mat(VIP_2P.MovDuration{ii}([16:end]))];
end       

save('Combine_2PWF_ForBhv_CueMvmDur.mat','Cue_all','RwdMvm_all','-append');

color_value = cbrewer('seq','PuBuGn',9);
color_value = color_value([4,5,7,9],:);
edges = [0:0.2:10];
figure; hold on; set(gcf,'color','w','pos',[200 200 500 200])
for ii = 1:4
    temp_n = histcounts(Cue_all{ii},edges,'Normalization','probability');
    temp_n = [temp_n,temp_n(end)];
    temp_median = nanmedian(Cue_all{ii});
    stairs(temp_n,'color',color_value(ii,:));
    xlim([1 51]);
    xticks([1:5:51]); xticklabels([0:1:10]);
    ylim([0 0.35]); yticks([0:0.05:0.35]); yticklabels([0:5:35]);
%     line([temp_median*10 temp_median*10],ylim,'color',color_value(ii,:),'linestyle',':');
    xlabel('Duration (s)');
    ylabel('Fraction (%)');
end
saveas(gcf,'Combine_CueDuration_distribution.fig','fig');
saveas(gcf,'Combine_CueDuration_distribution.png','png');
saveas(gcf,'Combine_CueDuration_distribution.pdf','pdf');

color_value = cbrewer('seq','PuBuGn',9);
color_value = color_value([4,5,7,9],:);
edges = [0:0.2:10];
figure; hold on; set(gcf,'color','w','pos',[200 200 500 200])
for ii = 1:4
    temp_n = histcounts(RwdMvm_all{ii},edges,'Normalization','probability');
    temp_n = [temp_n,temp_n(end)];
    temp_median = nanmedian(RwdMvm_all{ii});
    stairs(temp_n,'color',color_value(ii,:));
    xlim([1 51]);
    xticks([1:5:51]); xticklabels([0:1:10]);
    ylim([0 0.15]); yticks([0:0.05:0.15]); yticklabels([0:5:25]);
%     line([temp_median*10 temp_median*10],ylim,'color',color_value(ii,:),'linestyle',':');
    xlabel('Duration (s)');
    ylabel('Fraction (%)');
end
saveas(gcf,'Combine_RwdMvmDuration_distribution.fig','fig');
saveas(gcf,'Combine_RwdMvmDuration_distribution.png','png');
saveas(gcf,'Combine_RwdMvmDuration_distribution.pdf','pdf');