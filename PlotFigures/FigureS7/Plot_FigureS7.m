%% Plot Figure S7C
% Example: CR_3672020-O_190429_ROI_Traces.mat
colors = distinguishable_colors(size(roi_trace_df_2{1},1),[1,1,1;0,0,0]);
figure; hold on; set(gcf,'color','white','position',[50 50 600 500]); colormap(colors);
area(-active_frames*200,'facecolor',[0.8 0.8 0.8]);
for ii_cell = 1:size(roi_trace_df_2{1},1)
    Trace = roi_trace_df_2{1}(ii_cell,:);
    Trace = medfilt1(Trace,5);
    Trace = zscore(Trace);
    plot(Trace-7*(ii_cell-1),'color',colors(ii_cell,:));
end
xlim([1 length(active_frames)]);
ylim([-190 15]);


