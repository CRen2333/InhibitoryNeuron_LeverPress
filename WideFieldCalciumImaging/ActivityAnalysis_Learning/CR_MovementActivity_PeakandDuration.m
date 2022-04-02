% Just use activity during 2 sec after movement onset
close all
clc

clearvars -except VIP SOM PV Stats

INs = {'VIP','SOM','PV'};

for curr_IN = 1:length(INs)
    
    clearvars -except VIP SOM PV INs curr_IN Stats ColorMap
    IN = INs{curr_IN};
    load(['Z:\People\Chi\WFLP_IN\' IN '\Craniotomy\GAP500LEEWAY150_THY1MASK\' IN '_Activity_THY1MASK_OPExclude.mat'], 'Cued_ROI_df_f_arranged');
    for roi = 1:16
        for curr_animal = 1:size(Cued_ROI_df_f_arranged{roi},2)
            for curr_session = 1:size(Cued_ROI_df_f_arranged{roi}{curr_animal},2)
                temp_matrix = Cued_ROI_df_f_arranged{roi}{curr_animal}{curr_session};
                temp_matrix = temp_matrix-repmat(nanmean(temp_matrix(:,5:9),2),1,76);
                Cued_ROI_df_f_arranged_sub{roi}{curr_animal}{curr_session} = temp_matrix;
                PeakValue{roi}{curr_animal}{curr_session} = nanmax(temp_matrix,[],2)';
                [~,I] = nanmax(temp_matrix,[],2);
                PeakValue_PosEnd{roi}{curr_animal}{curr_session} = (I>73)';
                halfPeak = PeakValue{roi}{curr_animal}{curr_session}/2;
                % interp to improve temporal resolution
                x = 1:76;
                xq = [1:0.1:76];
                temp_matrix_2 = [];
                bi_temp_matrix_2 = [];
                for kk = 1:size(temp_matrix,1)
                    temp_matrix_2(kk,:) = interp1(x,temp_matrix(kk,:),xq);
                    bi_temp_matrix_2(kk,:) = temp_matrix_2(kk,:)>halfPeak(kk);
                    % close gap < 100 ms and remove events < 100 ms
%                     bi_on = find(diff([0,bi_temp_matrix_2(kk,:),0])==1);
%                     bi_off = find(diff([0,bi_temp_matrix_2(kk,:),0])==-1);
%                     for jj = 1:length(bi_on)-1
%                         if bi_on(jj+1)-bi_off(jj)<30
%                             bi_temp_matrix_2(kk,bi_off(jj):bi_on(jj+1)) = 1;
%                         end
%                     end
%                     bi_on = find(diff([0,bi_temp_matrix_2(kk,:),0])==1);
%                     bi_off = find(diff([0,bi_temp_matrix_2(kk,:),0])==-1);
%                     for jj = 1:length(bi_on)
%                         if bi_off(jj)-bi_on(jj)<30
%                             bi_temp_matrix_2(kk,bi_on(jj):bi_off(jj)) = 0;
%                         end
%                     end
                    bi_temp_matrix_2(kk,1:81) = 0;
                    bi_temp_matrix_2 = bi_temp_matrix_2(:,1:751);
                    temp_on = find(diff(bi_temp_matrix_2(kk,:))==1,1,'first')+1;
                    if ~any(bi_temp_matrix_2(kk,temp_on:end)==0)
                        Duration_TillEnd{roi}{curr_animal}{curr_session}(1,kk) = 1;
                    else
                        Duration_TillEnd{roi}{curr_animal}{curr_session}(1,kk) = 0;
                    end
                end
                temp_on = find(diff(bi_temp_matrix_2,[],1)==1,1,'first')+1;
                
                Duration{roi}{curr_animal}{curr_session} = sum(bi_temp_matrix_2,2)'/299.8;
            end
        end
    end
    
    Stage_sessions = {[1],[2,3,4],[5,6,7,8],[9,10,11]};
    Stages = {'Naive','Early','Middle','Late'};
    for roi = 1:16
        for curr_animal = 1:size(Cued_ROI_df_f_arranged{roi},2)
            for stage_ii = 1:4
                PeakValue_Stage{roi}{curr_animal}{stage_ii} = cell2mat(PeakValue{roi}{curr_animal}(Stage_sessions{stage_ii}));
                PeakValue_PosEnd_Stage{roi}{curr_animal}{stage_ii} = cell2mat(PeakValue_PosEnd{roi}{curr_animal}(Stage_sessions{stage_ii}));
                Duration_Stage{roi}{curr_animal}{stage_ii} = cell2mat(Duration{roi}{curr_animal}(Stage_sessions{stage_ii}));
                Duration_TillEnd_Stage{roi}{curr_animal}{stage_ii} = cell2mat(Duration_TillEnd{roi}{curr_animal}(Stage_sessions{stage_ii}));
                PeakValue_Stage_mean{roi}(curr_animal,stage_ii) = nanmean(PeakValue_Stage{roi}{curr_animal}{stage_ii});
                Duration_Stage_median{roi}(curr_animal,stage_ii) = nanmedian(Duration_Stage{roi}{curr_animal}{stage_ii});
            end
        end
    end
    
    switch IN
        case 'VIP'
            VIP.PeakValue_Stage = PeakValue_Stage;
            VIP.Duration_Stage = Duration_Stage;
            VIP.PeakValue_PosEnd_Stage = PeakValue_PosEnd_Stage;
            VIP.Duration_TillEnd_Stage = Duration_TillEnd_Stage;
            VIP.PeakValue_Stage_mean = PeakValue_Stage_mean;
            VIP.Duration_Stage_median = Duration_Stage_median;
        case 'SOM'
            SOM.PeakValue_Stage = PeakValue_Stage;
            SOM.Duration_Stage = Duration_Stage;
            SOM.PeakValue_PosEnd_Stage = PeakValue_PosEnd_Stage;
            SOM.Duration_TillEnd_Stage = Duration_TillEnd_Stage;
            SOM.PeakValue_Stage_mean = PeakValue_Stage_mean;
            SOM.Duration_Stage_median = Duration_Stage_median;
        case 'PV'
            PV.PeakValue_Stage = PeakValue_Stage;
            PV.Duration_Stage = Duration_Stage;
            PV.PeakValue_PosEnd_Stage = PeakValue_PosEnd_Stage;
            PV.Duration_TillEnd_Stage = Duration_TillEnd_Stage;
            PV.PeakValue_Stage_mean = PeakValue_Stage_mean;
            PV.Duration_Stage_median = Duration_Stage_median;
    end
    
    % plot
    ROI = {'r-Visual','r-aS1BC','l-aS1BC','l-M1','PPC','l-Visual','r-S1HL','r-M1','l-S1HL','l-pS1BC','aRSC','r-pS1BC','pRSC','l-M1/S1FL','r-M1/S1FL','M2'};
    Ordered_ROI = {'M2','l-M1/S1FL','r-M1/S1FL','l-M1','r-M1','l-aS1BC','r-aS1BC','l-S1HL','r-S1HL','PPC','l-pS1BC','r-pS1BC','aRSC','pRSC','l-Visual','r-Visual'};
    for ii = 1:length(ROI)
        temp = cellfun(@(x) strcmp(x, Ordered_ROI{ii}), ROI);
        reordered_module(ii) = find(temp==1);
        clear temp
    end
    
    % Duration > 2 fraction
    for roi = 1:16
        for curr_animal = 1:size(Duration_TillEnd_Stage{roi},2)
            for stage_ii = 1:4
                Duration_TillEnd_fraction{stage_ii}(curr_animal,roi) = sum(Duration_TillEnd_Stage{roi}{curr_animal}{stage_ii})/length(Duration_TillEnd_Stage{roi}{curr_animal}{stage_ii});
            end
        end
    end
    Duration_TillEnd_fraction{stage_ii} = Duration_TillEnd_fraction{stage_ii}(:,reordered_module);
    
    position.Naive = [1:5:5*length(ROI)];
    position.Early = [2:5:5*length(ROI)];
    position.Middle = [3:5:5*length(ROI)];
    position.Late = [4:5:5*length(ROI)];
    switch IN
        case 'VIP'
            [CX] = cbrewer('div','PRGn',10);
            color_value = [CX(end-3,:); CX(end-2,:); CX(end-1,:); CX(end,:)];
            VIP.color_value = color_value;
        case 'SOM'
            [CX] = cbrewer('div','PRGn',10);
            color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)];  
            SOM.color_value = color_value;
        case 'PV'
            [CX] = cbrewer('div','PuOr',12);
            color_value = [CX(4,:); CX(3,:); CX(2,:); CX(1,:)]; 
            PV.color_value = color_value;
    end
    
    % reorganize
    clear temp_matrix temp_mean temp_std
    for ii = 1:4
        for roi = 1:16
            temp_matrix{ii}(:,roi) = PeakValue_Stage_mean{roi}(:,ii);
        end
        temp_mean(ii,:) = nanmean(temp_matrix{ii});
        temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
    end
        
    figure('position',[100,100,1200,200],'Color','w')
    hold on;
    bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
    bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
    bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
    bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
    for curr_ROI = 1:length(ROI)
        line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
        line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
        line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
        line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
    end
    xlim([-1,position.Late(end)+2]); ylim([0,0.045])
    set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
    xlabel('Cortical module','FontSize',10);ylabel('Peak df/f','FontSize',10)
    box off
    title([IN ' post-movement-onset activity Peak'],'FontSize',10)
    
    % Statistics, REM
    for curr_ROI = 1:length(ROI)
        temp_var = PeakValue_Stage_mean{curr_ROI};
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
        pValue = stats.pValue;
        p_rem_ROI(:,curr_ROI) = pValue(2);
    end
    p_rem_ROI_reorder = p_rem_ROI(reordered_module);
    [FDR] = mafdr(p_rem_ROI_reorder,'BHFDR', true);

    position.P_ROI = (position.Early + position.Middle)/2;
    hold on
    for curr_ROI = 1:length(ROI)
        if FDR(:,curr_ROI) < 0.0001 
            text(position.P_ROI(curr_ROI),0.042,'****','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.001
            text(position.P_ROI(curr_ROI),0.042,'***','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.01
            text(position.P_ROI(curr_ROI),0.042,'**','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.05
            text(position.P_ROI(curr_ROI),0.042,'*','HorizontalAlignment','center');
        end
    end

    savefig(gcf,[IN '_PostMovActivityPeak_subtract_baseline_OPExclude.fig']);
    saveas(gcf,[IN '_PostMovActivityPeak_subtract_baseline_OPExclude.png']);
    print([IN '_PostMovActivityPeak_subtract_baseline_OPExclude.pdf'],'-dpdf','-bestfit'); pause(1);
    
    clear temp_matrix temp_mean temp_std
    for ii = 1:4
        for roi = 1:16
            temp_matrix{ii}(:,roi) = Duration_Stage_median{roi}(:,ii);
        end
        temp_mean(ii,:) = nanmean(temp_matrix{ii});
        temp_std(ii,:) = nanstd(temp_matrix{ii})./sqrt(sum(~isnan(temp_matrix{ii})));
    end
        
    figure('position',[100,100,1200,200],'Color','w')
    hold on;
    bar(position.Naive,temp_mean(1,reordered_module),0.20,'FaceColor',color_value(1,:),'LineStyle','none')
    bar(position.Early,temp_mean(2,reordered_module),0.20,'FaceColor',color_value(2,:),'LineStyle','none')
    bar(position.Middle,temp_mean(3,reordered_module),0.20,'FaceColor',color_value(3,:),'LineStyle','none')
    bar(position.Late,temp_mean(4,reordered_module),0.20,'FaceColor',color_value(4,:),'LineStyle','none')
    for curr_ROI = 1:length(ROI)
        line([position.Naive(curr_ROI),position.Naive(curr_ROI)],[temp_mean(1,reordered_module(curr_ROI))-temp_std(1,reordered_module(curr_ROI)),temp_mean(1,reordered_module(curr_ROI))+temp_std(1,reordered_module(curr_ROI))],'Color',color_value(1,:),'LineWidth',1);
        line([position.Early(curr_ROI),position.Early(curr_ROI)],[temp_mean(2,reordered_module(curr_ROI))-temp_std(2,reordered_module(curr_ROI)),temp_mean(2,reordered_module(curr_ROI))+temp_std(2,reordered_module(curr_ROI))],'Color',color_value(2,:),'LineWidth',1);
        line([position.Middle(curr_ROI),position.Middle(curr_ROI)],[temp_mean(3,reordered_module(curr_ROI))-temp_std(3,reordered_module(curr_ROI)),temp_mean(3,reordered_module(curr_ROI))+temp_std(3,reordered_module(curr_ROI))],'Color',color_value(3,:),'LineWidth',1);
        line([position.Late(curr_ROI),position.Late(curr_ROI)],[temp_mean(4,reordered_module(curr_ROI))-temp_std(4,reordered_module(curr_ROI)),temp_mean(4,reordered_module(curr_ROI))+temp_std(4,reordered_module(curr_ROI))],'Color',color_value(4,:),'LineWidth',1);    
    end
    xlim([-1,position.Late(end)+2]); ylim([0,1.2])
    set(gca,'FontSize',8,'TickLength',[0.01 0.01],'LineWidth',1,'XTick',(position.Early+position.Middle)/2,'XTickLabel',Ordered_ROI);
    xlabel('Cortical module','FontSize',10);ylabel('df/f FWHM','FontSize',10)
    box off
    title([IN ' post-movement-onset activity Duration'],'FontSize',10)
    
    % Statistics, REM
    for curr_ROI = 1:length(ROI)
        temp_var = Duration_Stage_median{curr_ROI};
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
        pValue = stats.pValue;
        p_rem_ROI(:,curr_ROI) = pValue(2);
    end
    p_rem_ROI_reorder = p_rem_ROI(reordered_module);
    [FDR] = mafdr(p_rem_ROI_reorder,'BHFDR', true);

    position.P_ROI = (position.Early + position.Middle)/2;
    hold on
    for curr_ROI = 1:length(ROI)
        if FDR(:,curr_ROI) < 0.0001 
            text(position.P_ROI(curr_ROI),1.18,'****','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.001
            text(position.P_ROI(curr_ROI),1.18,'***','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.01
            text(position.P_ROI(curr_ROI),1.18,'**','HorizontalAlignment','center');
        elseif FDR(:,curr_ROI) < 0.05
            text(position.P_ROI(curr_ROI),1.18,'*','HorizontalAlignment','center');
        end
    end

    savefig(gcf,[IN '_PostMovActivityDuration_subtract_baseline_OPExclude.fig']);
    saveas(gcf,[IN '_PostMovActivityDuration_subtract_baseline_OPExclude.png']);
    print([IN '_PostMovActivityDuration_subtract_baseline_OPExclude.pdf'],'-dpdf','-bestfit'); pause(1);
    
end

save('MovementActivty_Peak_Duration.mat','VIP','SOM','PV','-v7.3');

