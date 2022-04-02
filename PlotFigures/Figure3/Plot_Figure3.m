load('Figure3.mat');

%% Plot Figure 3B-3D, example true and GLM-fitted traces for M1
colors_field = [0,0,0;...
    63.75*1.5,63.75*1.5,63.75*1.5;...
    166,73,90;...
    60,94,115;...
    217,181,150;...
    166,120,93]/255;
Stages = {'Naive','Early','Middle','Late'};
fields = {'True','Full','Mvm','Lick','Cue','Reward'}; % reorder

% VIP
figure; set(gcf,'color','w','position',[100 100 800 200]); hold on;
roi = 8; % right M1
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

% SOM
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

% PV
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

%% Plot Figure 3E, Activity during movements, true and GLM-fitted
% VIP
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
% SOM
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
% PV
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

%% Plot Figure 3F, Changes of activity during movements relative to naive stage, true and GLM-fitted
% VIP
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
% SOM
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
% PV
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
