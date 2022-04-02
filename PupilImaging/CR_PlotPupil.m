%% Plot trace first
load('TaskAnimal_Pupil_2sec.mat')
for ii = 1:length(PV.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre)
    temp = PV.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}{1};
    PV.PupilTrace_4Stages_Naive(ii,:) = nanmean(temp,2)';
    temp = cell2mat(PV.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(2:4));
    PV.PupilTrace_4Stages_Early(ii,:) = nanmean(temp,2)';
    temp = cell2mat(PV.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(6:8));
    PV.PupilTrace_4Stages_Middle(ii,:) = nanmean(temp,2)';
    temp = cell2mat(PV.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(9:11));
    PV.PupilTrace_4Stages_Late(ii,:) = nanmean(temp,2)';
end
for ii = 1:length(VIP.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre)
    temp = VIP.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}{1};
    VIP.PupilTrace_4Stages_Naive(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(2:4));
    VIP.PupilTrace_4Stages_Early(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(6:8));
    VIP.PupilTrace_4Stages_Middle(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(9:11));
    VIP.PupilTrace_4Stages_Late(ii,:) = nanmean(temp,2)';
end
for ii = 1:length(SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre)
    temp = SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}{1};
    SOM.PupilTrace_4Stages_Naive(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(2:4));
    SOM.PupilTrace_4Stages_Early(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(6:8));
    SOM.PupilTrace_4Stages_Middle(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre{ii}(9:11));
    SOM.PupilTrace_4Stages_Late(ii,:) = nanmean(temp,2)';
end


for ii = 1:length(VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre)
    temp = VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}{1};
    VIP_CTRL.PupilTrace_4Stages_Naive(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(2:4));
    VIP_CTRL.PupilTrace_4Stages_Early(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(6:8));
    VIP_CTRL.PupilTrace_4Stages_Middle(ii,:) = nanmean(temp,2)';
    temp = cell2mat(VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(9:11));
    VIP_CTRL.PupilTrace_4Stages_Late(ii,:) = nanmean(temp,2)';
end
for ii = 1:length(SOM_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre)
    temp = SOM_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}{1};
    SOM_CTRL.PupilTrace_4Stages_Naive(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(2:4));
    SOM_CTRL.PupilTrace_4Stages_Early(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(6:8));
    SOM_CTRL.PupilTrace_4Stages_Middle(ii,:) = nanmean(temp,2)';
    temp = cell2mat(SOM_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre{ii}(9:11));
    SOM_CTRL.PupilTrace_4Stages_Late(ii,:) = nanmean(temp,2)';
end


CX_pupil = cbrewer('div','RdBu',10);
colorvalue = [CX_pupil(end-3,:); CX_pupil(end-2,:); CX_pupil(end-1,:); CX_pupil(end,:)];
figure
set(gcf,'position',[50,50,200,200]); hold on;
temp = [VIP.PupilTrace_4Stages_Naive;SOM.PupilTrace_4Stages_Naive;PV.PupilTrace_4Stages_Naive];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(1,:),'linewidth',1);
temp = [VIP.PupilTrace_4Stages_Early;SOM.PupilTrace_4Stages_Early;PV.PupilTrace_4Stages_Early];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(2,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(2,:),'linewidth',1);
temp = [VIP.PupilTrace_4Stages_Middle;SOM.PupilTrace_4Stages_Middle;PV.PupilTrace_4Stages_Middle];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(3,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(3,:),'linewidth',1);
temp = [VIP.PupilTrace_4Stages_Late;SOM.PupilTrace_4Stages_Late;PV.PupilTrace_4Stages_Late];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(4,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(4,:),'linewidth',1);
ylim([-0.14,0.08]);xlim([1 38]);xticks([8,23,38]);xticklabels({'0','1','2'});
xlabel('Time (sec)'); ylabel({'norm. Pupil Dia.'});
line([8,8],ylim,'color','k','linestyle',':','linewidth',2);
axis square;
savefig(gcf,'TaskAnimal_PupilTrace_2sec.fig');pause(1);
saveas(gcf,'TaskAnimal_PupilTrace_2sec.png');pause(1);
saveas(gcf,'TaskAnimal_PupilTrace_2sec.pdf');pause(1);

figure
set(gcf,'position',[50,50,250,230]); hold on;
temp = [VIP_CTRL.PupilTrace_4Stages_Naive;SOM_CTRL.PupilTrace_4Stages_Naive];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(1,:));
temp = [VIP_CTRL.PupilTrace_4Stages_Early;SOM_CTRL.PupilTrace_4Stages_Early];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(2,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(2,:));
temp = [VIP_CTRL.PupilTrace_4Stages_Middle;SOM_CTRL.PupilTrace_4Stages_Middle];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(3,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(3,:));
temp = [VIP_CTRL.PupilTrace_4Stages_Late;SOM_CTRL.PupilTrace_4Stages_Late];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(4,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(4,:));
ylim([-0.08,0.08]);xlim([1 53]);xticks([8,23,38]);xticklabels({'0','1','2'});
yticks([-0.08:0.04:0.08]);
xlabel('Time from mvm. onset (sec)'); ylabel({'Pupil Dia.'; 'norm. to Naive pre-mvm. onset'});
line([8,8],ylim,'color','k','linestyle',':');
savefig(gcf,'NoTaskAnimal_PupilTrace.fig');
saveas(gcf,'NoTaskAnimal_PupilTrace.png');
saveas(gcf,'NoTaskAnimal_PupilTrace.pdf');

%% Pre, post, delta
for ii = 1:length(SOM.PupilDia_RwdMovOnset_norm_Cued_subtract_NaivePre)
    temp_1 = SOM.PupilDia_PreRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial{ii,1};
    SOM.PupilAve_4Stages_Pre(ii,1) = nanmean(temp_1);
    temp_2 = SOM.PupilDia_AfterRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial{ii,1};
    SOM.PupilAve_4Stages_Post(ii,1) = nanmean(temp_2);
    SOM.PupilAve_4Stages_Delta(ii,1) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(SOM.PupilDia_PreRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,2:4));
    SOM.PupilAve_4Stages_Pre(ii,2) = nanmean(temp_1);
    temp_2 = cell2mat(SOM.PupilDia_AfterRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,2:4));
    SOM.PupilAve_4Stages_Post(ii,2) = nanmean(temp_2);
    SOM.PupilAve_4Stages_Delta(ii,2) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(SOM.PupilDia_PreRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,6:8));
    SOM.PupilAve_4Stages_Pre(ii,3) = nanmean(temp_1);
    temp_2 = cell2mat(SOM.PupilDia_AfterRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,6:8));
    SOM.PupilAve_4Stages_Post(ii,3) = nanmean(temp_2);
    SOM.PupilAve_4Stages_Delta(ii,3) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(SOM.PupilDia_PreRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,9:11));
    SOM.PupilAve_4Stages_Pre(ii,4) = nanmean(temp_1);
    temp_2 = cell2mat(SOM.PupilDia_AfterRwdMovOnset_norm_Cued_subtract_NaivePre_EachTrial(ii,9:11));
    SOM.PupilAve_4Stages_Post(ii,4) = nanmean(temp_2);
    SOM.PupilAve_4Stages_Delta(ii,4) = nanmean(temp_2-temp_1);
end

for ii = 1:length(VIP_CTRL.PupilDia_QuaMovOnset_norm_subtract_NaivePre)
    temp_1 = VIP_CTRL.PupilDia_PreQuaMovOnset_norm_subtract_NaivePre_EachTrial{ii,1};
    VIP_CTRL.PupilAve_4Stages_Pre(ii,1) = nanmean(temp_1);
    temp_2 = VIP_CTRL.PupilDia_AfterQuaMovOnset_norm_subtract_NaivePre_EachTrial{ii,1};
    VIP_CTRL.PupilAve_4Stages_Post(ii,1) = nanmean(temp_2);
    VIP_CTRL.PupilAve_4Stages_Delta(ii,1) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(VIP_CTRL.PupilDia_PreQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,2:4));
    VIP_CTRL.PupilAve_4Stages_Pre(ii,2) = nanmean(temp_1);
    temp_2 = cell2mat(VIP_CTRL.PupilDia_AfterQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,2:4));
    VIP_CTRL.PupilAve_4Stages_Post(ii,2) = nanmean(temp_2);
    VIP_CTRL.PupilAve_4Stages_Delta(ii,2) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(VIP_CTRL.PupilDia_PreQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,6:8));
    VIP_CTRL.PupilAve_4Stages_Pre(ii,3) = nanmean(temp_1);
    temp_2 = cell2mat(VIP_CTRL.PupilDia_AfterQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,6:8));
    VIP_CTRL.PupilAve_4Stages_Post(ii,3) = nanmean(temp_2);
    VIP_CTRL.PupilAve_4Stages_Delta(ii,3) = nanmean(temp_2-temp_1);
    
    temp_1 = cell2mat(VIP_CTRL.PupilDia_PreQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,9:11));
    VIP_CTRL.PupilAve_4Stages_Pre(ii,4) = nanmean(temp_1);
    temp_2 = cell2mat(VIP_CTRL.PupilDia_AfterQuaMovOnset_norm_subtract_NaivePre_EachTrial(ii,9:11));
    VIP_CTRL.PupilAve_4Stages_Post(ii,4) = nanmean(temp_2);
    VIP_CTRL.PupilAve_4Stages_Delta(ii,4) = nanmean(temp_2-temp_1);
end

% 2 sec
temp_1 = nanmean(PV.PupilTrace_4Stages_Naive(:,1:7),2);
temp_2 = nanmean(PV.PupilTrace_4Stages_Naive(:,8:38),2);
PV.PupilAve_4Stages_Delta_2sec(:,1) = temp_2-temp_1;
temp_1 = nanmean(PV.PupilTrace_4Stages_Early(:,1:7),2);
temp_2 = nanmean(PV.PupilTrace_4Stages_Early(:,8:38),2);
PV.PupilAve_4Stages_Delta_2sec(:,2) = temp_2-temp_1;
temp_1 = nanmean(PV.PupilTrace_4Stages_Middle(:,1:7),2);
temp_2 = nanmean(PV.PupilTrace_4Stages_Middle(:,8:38),2);
PV.PupilAve_4Stages_Delta_2sec(:,3) = temp_2-temp_1;
temp_1 = nanmean(PV.PupilTrace_4Stages_Late(:,1:7),2);
temp_2 = nanmean(PV.PupilTrace_4Stages_Late(:,8:38),2);
PV.PupilAve_4Stages_Delta_2sec(:,4) = temp_2-temp_1;

temp_1 = nanmean(SOM.PupilTrace_4Stages_Naive(:,1:7),2);
temp_2 = nanmean(SOM.PupilTrace_4Stages_Naive(:,8:38),2);
SOM.PupilAve_4Stages_Delta_2sec(:,1) = temp_2-temp_1;
temp_1 = nanmean(SOM.PupilTrace_4Stages_Early(:,1:7),2);
temp_2 = nanmean(SOM.PupilTrace_4Stages_Early(:,8:38),2);
SOM.PupilAve_4Stages_Delta_2sec(:,2) = temp_2-temp_1;
temp_1 = nanmean(SOM.PupilTrace_4Stages_Middle(:,1:7),2);
temp_2 = nanmean(SOM.PupilTrace_4Stages_Middle(:,8:38),2);
SOM.PupilAve_4Stages_Delta_2sec(:,3) = temp_2-temp_1;
temp_1 = nanmean(SOM.PupilTrace_4Stages_Late(:,1:7),2);
temp_2 = nanmean(SOM.PupilTrace_4Stages_Late(:,8:38),2);
SOM.PupilAve_4Stages_Delta_2sec(:,4) = temp_2-temp_1;

temp_1 = nanmean(VIP.PupilTrace_4Stages_Naive(:,1:7),2);
temp_2 = nanmean(VIP.PupilTrace_4Stages_Naive(:,8:38),2);
VIP.PupilAve_4Stages_Delta_2sec(:,1) = temp_2-temp_1;
temp_1 = nanmean(VIP.PupilTrace_4Stages_Early(:,1:7),2);
temp_2 = nanmean(VIP.PupilTrace_4Stages_Early(:,8:38),2);
VIP.PupilAve_4Stages_Delta_2sec(:,2) = temp_2-temp_1;
temp_1 = nanmean(VIP.PupilTrace_4Stages_Middle(:,1:7),2);
temp_2 = nanmean(VIP.PupilTrace_4Stages_Middle(:,8:38),2);
VIP.PupilAve_4Stages_Delta_2sec(:,3) = temp_2-temp_1;
temp_1 = nanmean(VIP.PupilTrace_4Stages_Late(:,1:7),2);
temp_2 = nanmean(VIP.PupilTrace_4Stages_Late(:,8:38),2);
VIP.PupilAve_4Stages_Delta_2sec(:,4) = temp_2-temp_1;

figure
set(gcf,'position',[50,50,200,200]); hold on;
% temp = [VIP.PupilAve_4Stages_Pre;SOM.PupilAve_4Stages_Pre;PV.PupilAve_4Stages_Pre];
% temp = [VIP.PupilAve_4Stages_Post;SOM.PupilAve_4Stages_Post;PV.PupilAve_4Stages_Post];
temp = [VIP.PupilAve_4Stages_Delta;SOM.PupilAve_4Stages_Delta;PV.PupilAve_4Stages_Delta];
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
for ii = 1:size(temp,1)
    plot(temp(ii,:),'color',[0.5 0.5 0.5]);
end
plot(temp_1,'color',nanmean(colorvalue(2:3,:)),'linewidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',nanmean(colorvalue(2:3,:)),'linewidth',2);
%     plot(ii,temp_1(ii),'marker','.','markersize',20,'color',colorvalue(4,:));
end
xlim([0.5 4.5]);xticks([1:4]);xticklabels({'Naive','Early','Middle','Late'});
% ylim([-0.12 0.07]);yticks([-0.1:0.05:0.05]);ylabel({'Pupil Dia.'; 'norm. to Naive pre-mvm. onset'});
ylabel({'delta nom. Pupil Dia.'}); ylim([0 0.06]);
yticks([0:0.02:0.06]);
axis square;

Sessions = repmat([1:4],size(temp,1),1);
Sessions = Sessions(:);
Animals_test = [];
for curr_animal = 1:size(temp,1)
    Animals_test = [Animals_test;ones(1,size(temp,2))*curr_animal];
end
Animals_test = Animals_test(:);
y = temp(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;



x = [ones(1,size(temp,1)),ones(1,size(temp,1))*2,ones(1,size(temp,1))*3,ones(1,size(temp,1))*4];
y = temp(:);
temp_p = regstats(y,x);
temp_p = temp_p.tstat.pval(2);
p = friedman(temp,1)
text(3,0.055,'***','fontsize',16,'color',colorvalue(4,:));
% savefig(gcf,'TaskAnimal_PupilPre.fig');
% saveas(gcf,'TaskAnimal_PupilPre.png');
% saveas(gcf,'TaskAnimal_PupilPre.pdf');
% savefig(gcf,'TaskAnimal_PupilPost.fig');
% saveas(gcf,'TaskAnimal_PupilPost.png');
% saveas(gcf,'TaskAnimal_PupilPost.pdf');
% savefig(gcf,'TaskAnimal_PupilDelta.fig');
% saveas(gcf,'TaskAnimal_PupilDelta.png');
% saveas(gcf,'TaskAnimal_PupilDelta.pdf');
% 
% save('TaskAnimal_Pupil_2sec.mat','PV','SOM','VIP','-append');

