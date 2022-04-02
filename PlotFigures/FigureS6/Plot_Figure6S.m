%% Plot Figure S6B
load('FigureS6B.mat');
% Trace
figure
set(gcf,'position',[50,50,200,200]); hold on;
temp = PupilTrace_4Stages.Naive;
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(1,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(1,:),'linewidth',1);
temp = PupilTrace_4Stages.Early;
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(2,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(2,:),'linewidth',1);
temp = PupilTrace_4Stages.Middle;
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(3,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(3,:),'linewidth',1);
temp = PupilTrace_4Stages.Late;
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
h = area([(temp_1-temp_2)',(2*temp_2)']);
set(h(1),'EdgeColor','none','FaceColor','none');
set(h(2),'EdgeColor','none','FaceColor',colorvalue(4,:),'FaceAlpha',0.3);
plot(temp_1,'color',colorvalue(4,:),'linewidth',1);
ylim([-0.14,0.08]);xlim([1 38]);xticks([8,23,38]);xticklabels({'0','1','2'});
xlabel('Time (s)'); ylabel({'Norm. pupil diameter'});
line([8,8],ylim,'color','k','linestyle',':','linewidth',2);
axis square;

% Pupil dilation during movements
figure
set(gcf,'position',[50,50,200,200]); hold on;
temp = PupilDilation_4Stages;
temp_1 = nanmean(temp);
temp_2 = nanstd(temp)./sqrt(sum(~isnan(temp(:,1))));
for ii = 1:size(temp,1)
    plot(temp(ii,:),'color',[0.5 0.5 0.5]);
end
plot(temp_1,'color',nanmean(colorvalue(2:3,:)),'linewidth',2);
for ii = 1:4
    line([ii,ii],[temp_1(ii)-temp_2(ii),temp_1(ii)+temp_2(ii)],'color',nanmean(colorvalue(2:3,:)),'linewidth',2);
end
xlim([0.5 4.5]);xticks([1:4]);xticklabels({'Naive','Early','Middle','Late'});
ylabel({'Norm. pupil dilation'}); ylim([0 0.06]);
yticks([0:0.02:0.06]);
axis square;
% Stats
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