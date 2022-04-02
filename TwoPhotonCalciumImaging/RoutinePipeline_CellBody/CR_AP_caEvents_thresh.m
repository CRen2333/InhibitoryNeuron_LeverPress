function [caEvents]= CR_AP_caEvents_thresh(df_matrix,thresh,method,framerate,noise_est_window_min,figureOn)
% caEvent_matrix = AP_caEvents_thresh(df_matrix,thresh,method)
% Create a matrix of calcium events from raw df/f traces
% 'thresh' - noise multiplier threshold
%
% method: 
% 0 (default) = zero inactive portions of the trace
% 1 = zero the trace when df/f is falling
% 2 = retain max-min df of increasing active df/f portions
% 
% modify:
% HL 16-09-06
% also return the threshold as a measure of SNR
%%

if nargin < 3
    method = 0;
end
num_cells = size(df_matrix,1);
noiseEst = nan(num_cells,1);
% Discard cells that contain all NaNs and zero occasional NaNs (is this ok?)
discard_cells = false(num_cells,1);
discard_cells(all(isnan(df_matrix),2)) = true;
df_matrix(isnan(df_matrix)) = 0;

use_cells = find(~discard_cells);

% Initialize calcium events matrix
caEvents = nan(size(df_matrix));

below_zero_trace_all = nan(size(df_matrix));
below_zero_trace_all(df_matrix < 0) = df_matrix (df_matrix < 0);
clear noise_est_all
% Try dynamic noise_est
noise_est_window = round(noise_est_window_min*60*framerate/2);
for curr_frame = 1:size(df_matrix,2)
    if curr_frame <= noise_est_window
        noise_est_all(:,curr_frame) = nanstd([below_zero_trace_all(:,1:min(curr_frame+noise_est_window,size(df_matrix,2))) -below_zero_trace_all(:,1:min(curr_frame+noise_est_window,size(df_matrix,2)))],[],2);
    elseif curr_frame > size(df_matrix,2)- noise_est_window
        noise_est_all(:,curr_frame) = nanstd([below_zero_trace_all(:,max(1,curr_frame-noise_est_window):end) -below_zero_trace_all(:,max(1,curr_frame-noise_est_window):end)],[],2);
    else
        noise_est_all(:,curr_frame) = nanstd([below_zero_trace_all(:,curr_frame-noise_est_window:curr_frame+noise_est_window) -below_zero_trace_all(:,curr_frame-noise_est_window:curr_frame+noise_est_window)],[],2);
    end
end
% refine edge
noise_est_all(:,1:round(framerate/2)) = repmat(noise_est_all(:,round(framerate/2)+1),1,round(framerate/2));
noise_est_all(:,end-round(framerate/2)+1:end) = repmat(noise_est_all(:,end-round(framerate/2)),1,round(framerate/2));

if figureOn
    figure;
    hold on;
end
for curr_cell = use_cells'
    
    curr_trace = df_matrix(curr_cell,:);
    noise_est = noise_est_all(curr_cell,:);
    % plot(curr_trace-5*(curr_cell-1),'color',[0.5,0.5,0.5]);
    temp_smooth_30 = smooth(curr_trace,30,'loess')';

    %     noise_est = mean(abs(df_matrix(curr_cell,:) - temp_smooth_30));
    % NEW NOISE EST: the std of (mirrored) below-zero df/f values
    
    % return noise_est 16-9-6 HL
    % noiseEst(curr_cell,1) = noise_est;
    
    thresh_lo = temp_smooth_30 > noise_est*1;
    thresh_hi = temp_smooth_30 > noise_est*thresh;
    
    % fill in hi-threshold portions where lo threshold is not crossed (dips
    % down but not to baseline, so continuously active portion)
    
    % find edges of long-smooth above-thresh periods
    thresh_lo_start = diff([0 thresh_lo 0]) == 1;
    thresh_lo_stop = diff([0 thresh_lo 0]) == -1;
    thresh_hi_start = diff([0 thresh_hi 0]) == 1;
    thresh_hi_stop = diff([0 thresh_hi 0]) == -1;
    
    thresh_hi_start_idx = find(thresh_hi_start);
    thresh_hi_stop_idx = find(thresh_hi_stop);
    
    thresh_lo_hi_smooth_idx = arrayfun(@(x,y) ...
        x-find(thresh_lo_start(x:-1:1),1)+1:y, ...
        thresh_hi_start_idx,thresh_hi_stop_idx,'uni',false);
    
    % don't bother trying to find activity that starts before imaging or
    % continues after last frame
    exclude_activity = cellfun(@(x) any(x(x <= 1 | ...
        x >= length(curr_trace))),thresh_lo_hi_smooth_idx);
    
    % refine start times of activity to when df/f goes above hi thresh
    thresh_lo_hi_raw_idx = cellfun(@(x) x(1) +...
        find(curr_trace(x:end) > noise_est(x:end)*thresh,1) - 1: ...
        x(end),thresh_lo_hi_smooth_idx(~exclude_activity),'uni',false);
    
    % again, filter out events that overlap with beginning or end
    exclude_activity_2 = cellfun(@(x) any(x(x <= 1 | ...
        x >= length(curr_trace))),thresh_lo_hi_raw_idx);
    thresh_lo_hi_raw_idx(exclude_activity_2) = [];
    
    % activity duration
    activity_duration = cellfun(@length, thresh_lo_hi_raw_idx);
    thresh_lo_hi_raw_idx(activity_duration<framerate/2) = [];
    
    % find continuous active portions after this process      
    active_trace = zeros(size(curr_trace));
    active_trace(horzcat(thresh_lo_hi_raw_idx{:})) = true;
    
    active_idx = arrayfun(@(x,y) x:y,find(diff([0 active_trace 0]) == 1), ...
        find(diff([0 active_trace 0]) == -1),'uni',false);
    
%     aaa = nan(size(curr_trace));
%     aaa(logical(active_trace)) = curr_trace(logical(active_trace));
%     plot(aaa-5*(curr_cell-1),'r');
    
    % Create trace from chosen method
    switch method
        % 0 (default) = retain over-threshold, zero others
        case 0
            over_thresh_frames = horzcat(active_idx{:});
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = curr_trace(over_thresh_frames);
            
        % 1 = retain over-threshold & increasing, zero others
        case 1
            thresh_lo_hi_rising = cellfun(@(x) x(diff([0 temp_smooth_30(x)]) > 0), ...
                active_idx,'uni',false);
            
            over_thresh_frames = horzcat(thresh_lo_hi_rising{:});
            
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = curr_trace(:,over_thresh_frames);
            
        % 2 = retain ddf of over-threshold & increasing portions
        case 2
            % find rising portions
            thresh_lo_hi_rising = cellfun(@(x) x(diff([0 temp_smooth_30(x)]) > 0), ...
                active_idx,'uni',false);
            % split rising portions
            active_rising_trace = zeros(size(curr_trace));
            active_rising_trace(horzcat(thresh_lo_hi_rising{:})) = true;            
            active_rising_idx = arrayfun(@(x,y) x:y,find(diff([0 active_rising_trace 0]) == 1), ...
                find(diff([0 active_rising_trace 0]) == -1),'uni',false);
    
            thresh_lo_hi_rising_ddf = cellfun(@(x) repmat((max(curr_trace(x)) - ...
                min(curr_trace(x))),1,length(x)),active_rising_idx,'uni',false);
            
            over_thresh_frames = horzcat(active_rising_idx{:});
            over_thresh_values = horzcat(thresh_lo_hi_rising_ddf{:});
            
            caEvents(curr_cell,:) = 0;
            caEvents(curr_cell,over_thresh_frames) = over_thresh_values;
    end
    if figureOn
        plot(curr_trace-5*(curr_cell-1),'color',[0.5,0.5,0.5]);
        aaa = nan(size(curr_trace));
        aaa(logical(active_trace)) = curr_trace(logical(active_trace));
        plot(aaa-5*(curr_cell-1),'r');
        plot(caEvents(curr_cell,:)-5*(curr_cell-1),'color','b');
    end
        
end










