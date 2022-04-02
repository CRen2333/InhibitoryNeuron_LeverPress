function [roi_trace,roi_trace_bg] = CR_AP_Rig3_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp)
%[roi_trace,roi_trace_bg] = AP_getConcatTrace_continuous_batch(roi_filename,tiff_path,interp_sub,local_comp)
%
% Get fluorescence traces of ROIs (and background ROIs, if available) of
% movies which were taken continuously. Assumes new TIFF header format
% where frame number is logged.
%
% interp_sub - get fluorescence interpolated across ROIs as background
%
% local_comp - if being done locally, copy file to local drive temporarily
% and load/process

% Get TIFF filenames
tiff_dir = dir([tiff_path filesep '*.tif']);
tiff_files = sort({tiff_dir.name});

% Get summary file
fns_summary = dir(fullfile(tiff_path,'*summary.mat'));
fns_summary = sort({fns_summary.name});
if length(fns_summary) ~= length(tiff_files)
   error('summary file number does not match tiff file number');
end

% append directory to tiff filename
tiff_filenames = cellfun(@(x) [tiff_path filesep x], ...
    tiff_files,'UniformOutput',false);
fns_summary = cellfun(@(x) [tiff_path filesep x], ...
    fns_summary,'UniformOutput',false);

% if on local computer, create directory for temporary files
if local_comp
    local_dir = 'C:\temp_files';
    local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
    if ~isdir(local_dir)
        mkdir(local_dir)
    end
    % delete local temp file if it exists
    if exist(local_filename)
        delete(local_filename)
    end
end

% load ROIs
load(roi_filename,'-MAT')

% get size of frame from first frame of first file and remove zero columns
% and rows due to motion correction
imageinfo=imfinfo(tiff_filenames{1},'tiff');
temp_summary_file = load(fns_summary{1},'-mat');
if isfield(temp_summary_file, 'warp')
    M = sum(temp_summary_file.nonzero_column);
    N = sum(temp_summary_file.nonzero_row);
    transform = 'affine';
else
    M=imageinfo(1).Width;
    N=imageinfo(1).Height;
end

% get total number of frames from the last frame of the last file
imageinfo = imfinfo(tiff_filenames{end},'tiff');
[img_parameter img_value] = strread( ...
    imageinfo(end).ImageDescription,'%s %s', 'delimiter','=\n');
frame_tag_idx = cellfun(@(x) ~isempty(strfind(x,'frameNumbers')),img_parameter);
total_frames = str2num(img_value{frame_tag_idx});

% create masks from ROIs
cellMask = false(N*M,length(polygon.ROI));
for n_polygon = 1:length(polygon.ROI);
    %[HL 2016-5-17] check if ROI_mask exist => use the fixed masks by make_bg_ring.m
    if isfield(polygon,'ROI_mask')
        cellMask(:,n_polygon) = polygon.ROI_mask{n_polygon}(:); 
    elseif ~isfield(polygon,'autosort');
        temp_mask = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
        cellMask(:,n_polygon) = temp_mask(:);
    elseif isfield(polygon,'autosort');
        temp_mask = polygon.autosort(:,:,n_polygon);
        cellMask(:,n_polygon) = temp_mask(:);
    end
end
% Normalize mask values to prepare for matrix multiplation
cellMask_norm = sparse(bsxfun(@times,cellMask,1./sum(cellMask,1)));%sparse matrix can speed up the process

% Grab mask for background ROIs, if available
if isfield(polygon,'bgROI')
    for n_polygon = 1:length(polygon.bgROI);
        cellMask_bg(:,n_polygon) = polygon.bgROI_mask{n_polygon}(:);
    end
    % Normalize mask values to prepare for matrix multiplation
    cellMask_bg_norm = sparse(bsxfun(@times,cellMask_bg,1./sum(cellMask_bg,1)));
end

% get the traces for the ROIs
disp('Getting and concatenating traces');
roi_trace = nan(length(polygon.ROI),total_frames);
roi_trace_bg = nan(length(polygon.ROI),total_frames);

% frame_tag_previous = 0;
for i = 1:length(tiff_filenames);
    
    % if on local computer, copy image to temporary file, otherwise load
    if ~local_comp
        img_filename = [tiff_filenames{i}];
    else
        disp('Copying movie to local drive...')
        tic
        copyfile(tiff_filenames{i},local_filename);
        toc
        disp('Done.')
        img_filename = local_filename;
    end
    % summary file
    temp_summary_file = load(fns_summary{i},'-mat');

    % get traces from current tiff file
    imageinfo=imfinfo(img_filename,'tiff');
    numframes=length(imageinfo);
%     M=imageinfo(1).Width;
%     N=imageinfo(1).Height;
    
    % get number of each frame from header
    [img_parameter img_value] = arrayfun(@(x) strread( ...
        imageinfo(x).ImageDescription,'%s %s', 'delimiter','=\n'),1:length(imageinfo),'uni',false);
    frame_tag_idx = cellfun(@(x) cellfun(@(x) ...
        ~isempty(strfind(x,'frameNumbers')),x),img_parameter,'uni',false);
    frame_tag = cellfun(@(value,idx) str2num(value{idx}),img_value,frame_tag_idx);
%     frame_tag = frame_tag + frame_tag_previous;
%     frame_tag_previous = frame_tag(end);
    
    % Load and get activity frame-by-frame to allow for any computer to do
    % many of these simultanously (as opposed to loading in whole movie,
    % which is probably faster)
    tic
    disp(['Loading file, getting activity (' img_filename ')']);
    if interp_sub
        allROI_mask = reshape(any(cellMask,2),M,N);
    end
    for loadframe = 1:numframes
        % Load current frame
        curr_frame = [];
        curr_frame = imread(img_filename,'tiff',loadframe,'Info',imageinfo);
        
        % warp
        if isfield(temp_summary_file, 'warp')
            curr_frame = curr_frame(temp_summary_file.nonzero_row,temp_summary_file.nonzero_column);
            if temp_summary_file.piece == 1
                curr_frame = spatial_interp(double(curr_frame), temp_summary_file.warp, 'cubic', transform, 1:M, 1:N);
            end
        end
        % Get activity from current frame
        roi_trace(:,frame_tag(loadframe)) = transpose(double(curr_frame(:)')*cellMask_norm);
        
        %         %%%%%%%%%%%%%%%% TEMPORARY TROUBLESHOOT
        %         if any(isnan(roi_trace(:,frame_tag(loadframe))))
        %             keyboard
        %         end

        % Get background activity from current frame if selected
        if interp_sub
            curr_frame_interp = double(curr_frame);
            curr_frame_interp(allROI_mask) = NaN;
            curr_frame_interp = inpaint_nans(curr_frame_interp,4);
            roi_trace_bg(:,frame_tag(loadframe)) = transpose(curr_frame_interp(:)'*cellMask_norm);
        elseif isfield(polygon,'bgROI')
            roi_trace_bg(:,frame_tag(loadframe)) = transpose(double(curr_frame(:)')*cellMask_bg_norm);
        end

    end
    disp('Done')
    toc
    
    % If on local computer, get rid of temporary file
    if local_comp
        delete(local_filename);
    end
end

disp('Finished getting ROI traces');

