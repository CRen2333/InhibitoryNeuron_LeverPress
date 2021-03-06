
% Get image information
tiff_fullfilename = get(handles.tiffFile,'string');
imageinfo=imfinfo(tiff_fullfilename,'tiff');
numframes=length(imageinfo);
M=imageinfo(1).Width;
N=imageinfo(1).Height;

% Create mask for selected ROI, get trace
for n_polygon = 1:length(polygon.ROI);
    if ~isfield(polygon,'autosort');
        cellMask(:,:,n_polygon) = poly2mask(polygon.ROI{n_polygon}(:,1)',polygon.ROI{n_polygon}(:,2)',N,M);
    elseif isfield(polygon,'autosort');
        cellMask(:,:,n_polygon) = polygon.autosort(:,:,n_polygon);
    end
end

w = waitbar(0,'Saving ROI Traces');
for frame=1:numframes
    im_s=double(imread(tiff_fullfilename,frame));
    for n_polygon = 1:length(polygon.ROI)
        normactivityMap = (im_s.*cellMask(:,:,n_polygon))./sum(sum(cellMask(:,:,n_polygon)));
        traceVal(n_polygon,1) = sum(normactivityMap(:));
    end
    roi_trace(:,frame) = traceVal;
    waitbar(frame/numframes,w,'Saving ROI Traces'); % not with parfor
end

% Don't normalize on save at the moment

% waitbar(50,w,'Normalizing Traces')
% norm_fit_obj = gmdistribution.fit(cellTrace',2); % there are two hidden curves
% dx = (max(cellTrace)-min(cellTrace))/length(cellTrace);
% x = min(cellTrace):dx:max(cellTrace);
% x = x(1:end-1); % it adds an extra, not sure why but can fix later
% norm_fit = pdf(norm_fit_obj,x');
% baseline_firing = x(norm_fit == max(norm_fit));

% % change to df/f make peakin relative to df
% cellTrace = cellTrace./baseline_firing;

close(w);

    

