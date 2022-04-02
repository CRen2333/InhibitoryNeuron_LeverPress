%% Plot Figure S1A, time from cue to rewarded movement onset
load('FigureS1AS1BS1C.mat')
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.Cue_to_CuedRewardedMov;
temp_var(Combine.TrialNum<3) = nan;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:22
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 22.5]);xticks([5:5:20]);xlabel('Session');
ylim([0 3]); yticks([0:1:3]); ylabel('Time (s)');
set(gca,'XColor','k','YColor','k');
line([2.5 2.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([8.5 8.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([16.5 16.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,2.5,['Session: ' num2str(pValue(2))],'fontsize',6);

%% Plot Figure S1B, time from rewarded movement onset to reward
load('FigureS1AS1BS1C.mat')
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.RwdMVM_Rwd;
temp_var(Combine.TrialNum<3) = nan;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:22
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 22.5]);xticks([5:5:20]);xlabel('Session');
ylim([0.3 0.7]); yticks([0.3:0.1:0.7]); ylabel('Time (s)');
set(gca,'XColor','k','YColor','k');
line([2.5 2.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([8.5 8.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([16.5 16.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,0.65,['Session: ' num2str(pValue(2))],'fontsize',6);

%% Plot Figure S1C, time from cue to reward
load('FigureS1AS1BS1C.mat')
figure; hold on; set(gcf,'color','w','position',[200,200,210,200]);
temp_var = Combine.Cue_to_Reward;
temp_var(Combine.TrialNum<3) = nan;
temp_11 = nanmean(temp_var);
temp_12 = nanstd(temp_var)./sqrt(sum(~isnan(temp_var)));
plot(temp_11,'color','k','linewidth',1);
for ii = 1:22
    line([ii ii],[temp_11(ii)-temp_12(ii),temp_11(ii)+temp_12(ii)],'color','k','linewidth',1);
end
xlim([0.5 22.5]);xticks([5:5:20]);xlabel('Session');
ylim([0 4]); yticks([0:1:4]); ylabel('Time (s)');
set(gca,'XColor','k','YColor','k');
line([2.5 2.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([8.5 8.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
line([16.5 16.5],ylim,'linestyle',':','color',[0.5,0.5,0.5]);
axis square
% Stats
Sessions = repmat([1:22],size(temp_var,1),1);
Sessions = Sessions(:);
Animals_test = repmat([1:size(temp_var,1)]',1,22);
Animals_test = Animals_test(:);
y = temp_var(:);
tbl = table(Animals_test,Sessions,y);
tbl.Animals_test = nominal(tbl.Animals_test);
lme = fitlme(tbl,'y ~ 1 + Sessions + (1|Animals_test) + (Sessions-1|Animals_test)');
[beta,betanames,stats] = fixedEffects(lme);
pValue = stats.pValue;
p_val.LeverCorr_Reward_across = pValue(2);
text(3,0.5,['Session: ' num2str(pValue(2))],'fontsize',6);

%% Plot Figure S1D, distribution of rewarded movement duration
load('FigureS1D.mat')
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
    xlabel('Duration (s)');
    ylabel('Fraction (%)');
end

%% Plot figure S1E, example traces
% Example: CR_3218181-O, 100 Hz
% left to right, top to bottom
% #1
ii = 3; jj = 44;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
% #2
ii = 2; jj = 11;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
% #3
ii = 2; jj = 25;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
% #4
ii = 3; jj = 6;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
% #5
ii = 2; jj = 37;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
% #6
ii = 2; jj = 35;
temp_trace = CuedMov_SingleAnimal{1,ii}.LeverTrace_Reward(:,jj);
temp_trace = resample(temp_trace,1,100); % 10k --> 100 Hz
rwdmvm_index = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(:,5)==1;
rwd_time = CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,7)-CuedMov_SingleAnimal{ii}.Cued_MovOnset_Info_All(rwdmvm_index,2);
rwd_time = rwd_time(jj)*100;
figure; set(gcf,'color','w','position',[100,100,200,200]); hold on;
plot(temp_trace,'color',[0.7 0.7 0.7],'linewidth',0.5);
plot(rwd_time+50,-0.5,'b.');
ylim([-1.5 -0.2]);xlim([10 240]);
line([51 51],ylim,'color','k','linestyle',':');
line(xlim,[-0.77 -0.77],'color','b','linestyle',':');
line(xlim,[-0.77-0.25 -0.77-0.25],'color','b','linestyle',':');
xticks([]);
