function CR_warp_Npatches_RH4tscc(islocal, foldername, local2savefoldername, rl_RemovingPix, block_size, overlap_npix, n_split4warpinit)
% warp_routine (foldername, local2savefoldername)
% routine warp function run after motion correction
% It looks for motioncorrected_tiff folder under foldername, warp frames in _align file to the
% first frame, and save .warp file in the folder of input. 
% It generate a sum_50_warp.tif file of warped images, using summary files
% It also apply warp to each tif file under motioncorrected_tiff folder and
% save info into summary file [if warped, and warp used]
% INPUT:
%   foldername: -s, the data folder contains a 'motioncorrected_tiff' folder
%   local2savefoldername: -s, optional, local folder to save warped sum file if
%   provided
% OUTPUT:
%   none
%%

if ispc && islocal
   %make temp folder and files
   local_dir = 'C:\temp_files';
   local_filename = [local_dir filesep 'temp_img_' datestr(now,'yymmddHHMMSSFFF') '.tif'];
   if ~isdir(local_dir)
       mkdir(local_dir)
   end
   % delete local temp file if it exists
   if exist(local_filename,'file')
       delete(local_filename)
   end
end

mkdir(local2savefoldername, 'summed');

%find align file % get names to save
fn_align = dir(fullfile(foldername,'summed','*align.tif'));
fn_sum = dir(fullfile(foldername,'summed','*summed_50.tif'));
if length(fn_align)~=1 || length(fn_sum)~=1
    disp(fn_align); disp(fn_sum);
    error('align file is not the only one')
end
align_name = fn_align.name;
fn_sum_name = fn_sum.name;
fn_align = fullfile(foldername,'summed',align_name);
% fn_align_save = [fn_align(1:end-4),'_warped.tif'];
fn_align_save = fullfile(local2savefoldername,'summed',[align_name(1:end-4),'_warped.tif']);
% fn_warp2save = fullfile(foldername,[fn_sum.name(1:end-4),'_warped.warp']);
fn_warp2save = fullfile(local2savefoldername,[fn_sum_name(1:end-4),'_warped.warp']);
fn_sum = fullfile(foldername,'summed',fn_sum_name);
% fn_sum2save = [fn_sum(1:end-4),'_warped.tif'];
fn_sum2save = fullfile(local2savefoldername,'summed',[fn_sum_name(1:end-4),'_warped.tif']);

%default parameters for warping
n_ch = 1;
ch = 1;
% n_shrink = 2;
% margin_ratio = 0.125;
n_depth = 1;
n_step = 50;
transform = 'affine';
normalize_r1 = 1;
noramlize_r2 = 10; % default 32. Decrease this number if image quality is bad (dark). For example, 10.
% noramlize_r2 = 32; % Use when you use non-shrunk images for estimating
% affine transformation?
noramlize_offset = 50;

%load tif files
disp('load align file')
tic
if ispc && islocal %copy to local
    disp('Copying align to local drive...')
    img_filename = fn_align;
    disp(img_filename);
    copyfile(img_filename,local_filename);
    disp('Done.')
    img_filename = local_filename;
    disp(img_filename);
    [stack_align,info_align,frame_tag] = read_tiff(img_filename, ch, n_ch);
else
    [stack_align,info_align,frame_tag] = read_tiff(fn_align, ch, n_ch);
end

stack_align = stack_align(:,1+rl_RemovingPix:end-rl_RemovingPix,:); % remove n pixels from left and right edges. 

zero_zone = min(stack_align(:,:,2:end-1),[],3); % zero-zone created by motion correction, but excluding 1st and last
nonzero_row = sum(zero_zone,2)~=0;
nonzero_column = sum(zero_zone,1)~=0;
stack_align = stack_align(nonzero_row,nonzero_column,:); % remove zero-zone


toc

% minpixvalue = min(stack_align(:,:,2:end-1),[],3);
% prctile(minpixvalue(:),50);

ny = size(stack_align,1);
nx = size(stack_align,2);
nz = size(stack_align,3);

% overlap_npix = round((overlap_pix-1)/2);
qN_x = cell(block_size,block_size);
for i1 = 1:block_size
    qN_x{i1,1} = 1:ceil(nx/block_size)+overlap_npix;
    for i2 = 2:block_size-1
        qN_x{i1,i2} = (i2-1)*ceil(nx/block_size)-overlap_npix:i2*ceil(nx/block_size)+overlap_npix;
    end
    qN_x{i1,block_size} = (block_size-1)*ceil(nx/block_size)-overlap_npix:nx;
end
qN_y = cell(block_size,block_size);
for i2 = 1:block_size
    qN_y{1,i2} = 1:ceil(ny/block_size)+overlap_npix;
    for i1 = 2:block_size-1
        qN_y{i1,i2} = (i1-1)*ceil(ny/block_size)-overlap_npix:i1*ceil(ny/block_size)+overlap_npix;
    end
    qN_y{block_size,i2} = (block_size-1)*ceil(ny/block_size)-overlap_npix:ny;
end

% first get the warp transform, then median filter it to use
disp('warping align image')
tic
%get warp transform
warp_cell = cell(block_size,block_size,nz);
% for i1 = 1:block_size
%     for i2 = 1:block_size
%         for i3 = 1:nz
%             warp_cell{i1,i2,i3} = NaN(2,3);
%         end
%     end
% end

if mod(n_split4warpinit,2) ~= 1
    error('Please switch to split code');
end

template_ave_length = 10;
template_im = mean(stack_align(:,:,2:2+template_ave_length),3); % default (:,:,1)

template_shrink = NaN(ny, nx, n_split4warpinit);
template_shrink = imnormalize(template_im,normalize_r1,noramlize_r2,noramlize_offset);
% template_shrink = normalize_image(template_im);

% mask = ones(size(template_shrink,1),size(template_shrink,2));
moving_shrink = NaN(ny,nx,nz);
for i = 1:nz
%     if i < 17
%         moving_shrink(:,:,i) = int16(imnormalize(mean(stack_align(:,:,2:16),3),normalize_r1,noramlize_r2,noramlize_offset));
% %        moving_shrink(:,:,i) = int16(normalize_image(mean(stack_align(:,:,2:16),3)));
%     elseif i > nz - 15
%         moving_shrink(:,:,i) = int16(imnormalize(mean(stack_align(:,:,nz-15:nz-1),3),normalize_r1,noramlize_r2,noramlize_offset));
% %        moving_shrink(:,:,i) = int16(normalize_image(mean(stack_align(:,:,nz-15:nz-1),3)));
%     else
%         moving_shrink(:,:,i) = int16(imnormalize(mean(stack_align(:,:,i-10:i+10),3),normalize_r1,noramlize_r2,noramlize_offset));
% %        moving_shrink(:,:,i) = int16(normalize_image(mean(stack_align(:,:,i-10:i+10),3)));
%     end

    if i < 6
        moving_shrink(:,:,i) = int16(imnormalize(median(stack_align(:,:,2:12),3),normalize_r1,noramlize_r2,noramlize_offset));
    elseif i > nz - 5
        moving_shrink(:,:,i) = int16(imnormalize(median(stack_align(:,:,nz-10:nz),3),normalize_r1,noramlize_r2,noramlize_offset));
    else
        moving_shrink(:,:,i) = int16(imnormalize(median(stack_align(:,:,i-5:i+5),3),normalize_r1,noramlize_r2,noramlize_offset));
    end
    
%     if i < 6
%         moving_shrink(:,:,i) = int16(normalize_image(double(median(stack_align(:,:,2:12),3))));
%     elseif i > nz - 5
%         moving_shrink(:,:,i) = int16(normalize_image(double(median(stack_align(:,:,nz-10:nz),3))));
%     else
%         moving_shrink(:,:,i) = int16(normalize_image(double(median(stack_align(:,:,i-5:i+5),3))));
%     end
end

% [columnsInImage rowsInImage] = meshgrid(1:nx, 1:ny);
%     centerX = round(nx/2);
%     centerY = round(ny/2);
%     radius = round(nx/2);
%     mask = (rowsInImage - centerY).^2 ...
%         + (columnsInImage - centerX).^2 <= radius.^2;
% mask = ones(ny,nx);

affine_init = [1 0 0; 0 1 0];
warp4template = cell(block_size, block_size, n_split4warpinit);
for i1 = 1:block_size
    for i2 = 1:block_size
        warp4template{i1,i2,n_split4warpinit} = affine_init;
    end
end
rho = NaN(block_size,block_size,nz);

warp_range = nan(n_split4warpinit,2);
for i = 1:n_split4warpinit
    if i ==1
        warp_range(i,:) = [1,i*ceil(nz/n_split4warpinit)];
%     elseif i == n_split4warpinit
%         warp_range(i,:) = [(i-1)*ceil(nz/n_split4warpinit)+1, nz];
%     else
%         warp_range(i,:) = [(i-1)*ceil(nz/n_split4warpinit)+1, i*ceil(nz/n_split4warpinit)];
    end
end
    
n_warp2use4template = 7;
BadWarp_threshold = 50;
% out_parts = cell(block_size, block_size, n_split4warpinit);
for nth_warp = 1%[n_split4warpinit/2:-1:1,n_split4warpinit/2+1:n_split4warpinit]
    
    for i3 = warp_range(nth_warp,1):warp_range(nth_warp,2)
        parfor i1 = 1:block_size
            for i2 = 1:block_size
%                 [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
%                     ecc(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},i3),template_shrink(qN_y{i1,i2},qN_x{i1,i2},nth_warp),n_depth,n_step, transform, warp4template{i1,i2,nth_warp},mask(qN_y{i1,i2},qN_x{i1,i2}));
%                 [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
%                     ecc(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},i3),template_shrink(qN_y{i1,i2},qN_x{i1,i2},n_split4warpinit/2),n_depth,n_step, transform, warp4template{i1,i2,nth_warp}, mask(qN_y{i1,i2},qN_x{i1,i2}));
                  
                  [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
                    ecc_old(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},i3),template_shrink(qN_y{i1,i2},qN_x{i1,i2},n_split4warpinit),n_depth,n_step, transform, warp4template{i1,i2,nth_warp}(1:2,1:3));
%                 [~, warp_cell{i1,i2,i3}, ~, rho(i1,i2,i3)]= ...
%                     ecc_RH(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},i3),template_shrink(qN_y{i1,i2},qN_x{i1,i2},n_split4warpinit),n_depth,n_step, transform, warp4template{i1,i2,nth_warp}(1:2,1:3));
            end
        end
    end
    
    for i1 = 1:block_size
        for i2 = 1:block_size
            for i3 = warp_range(nth_warp,1):warp_range(nth_warp,2)
                if sum(sum(abs(warp_cell{i1,i2,i3})>BadWarp_threshold))
                    rho(i1,i2,i3) = NaN;
                end
                if isempty(warp_cell{i1,i2,i3})
                    warp_cell{i1,i2,i3} = NaN(2,3);
                end
                if warp_cell{i1,i2,i3}(1,1)<0.6 || warp_cell{i1,i2,i3}(2,2)<0.6 || warp_cell{i1,i2,i3}(1,1)>1.4 || warp_cell{i1,i2,i3}(2,2)>1.4 || rho(i1,i2,i3) <= 0.5
                    warp_cell{i1,i2,i3} = NaN(2,3);
                end
                if size(warp_cell{i1,i2,i3},1) == 3
                    warp_cell{i1,i2,i3}(3,:) = [];
                end
            end
        end
    end
    
    if (nth_warp ~= 1 && nth_warp ~= n_split4warpinit)
        for i1 = 1:block_size
            for i2 = 1:block_size
                if nth_warp <= n_split4warpinit/2
                    temp_id = find(~isnan(sum(sum(cell2mat(warp_cell(i1,i2,warp_range(nth_warp,1):warp_range(nth_warp,2))),1),2))==1,n_warp2use4template,'first');
                    if ~isempty(temp_id)
                        warp4template{i1,i2,nth_warp-1} = nanmedian(cell2mat(warp_cell(i1,i2,temp_id + warp_range(nth_warp,1) - 1)),3);
%                         out_parts{i1,i2,nth_warp-1} = spatial_interp(nanmean(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},(1:template_ave_length) + warp_range(nth_warp,1) - 1),3), warp4template{i1,i2,nth_warp-1}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                    else
                        warp4template{i1,i2,nth_warp-1} = warp4template{i1,i2,nth_warp};
%                         out_parts{i1,i2,nth_warp-1} = spatial_interp(nanmean(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},(1:template_ave_length) + warp_range(nth_warp,1) - 1),3), warp4template{i1,i2,nth_warp-1}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                    end
%                     out_parts{i1,i2,nth_warp-1}(out_parts{i1,i2,nth_warp-1}==0) = NaN;
                    if isempty(warp4template{i1,i2,nth_warp-1})
                        warp4template{i1,i2,nth_warp-1} = warp4template{i1,i2,nth_warp};
                    end
                else
                    temp_id = find(~isnan(sum(sum(cell2mat(warp_cell(i1,i2,warp_range(nth_warp,1):warp_range(nth_warp,2))),1),2))==1,n_warp2use4template,'last');
                    if ~isempty(temp_id)
                        warp4template{i1,i2,nth_warp+1} = nanmedian(cell2mat(warp_cell(i1,i2,temp_id + warp_range(nth_warp,1) - 1)),3);
%                         out_parts{i1,i2,nth_warp+1} = spatial_interp(nanmean(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},(1:template_ave_length) + warp_range(nth_warp,2) - template_ave_length),3), warp4template{i1,i2,nth_warp+1}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                    else
                        warp4template{i1,i2,nth_warp+1} = warp4template{i1,i2,nth_warp};
%                         out_parts{i1,i2,nth_warp+1} = spatial_interp(nanmean(moving_shrink(qN_y{i1,i2},qN_x{i1,i2},(1:template_ave_length) + warp_range(nth_warp,2) - template_ave_length),3), warp4template{i1,i2,nth_warp+1}, transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                    end
%                     out_parts{i1,i2,nth_warp+1}(out_parts{i1,i2,nth_warp+1}==0) = NaN;
                    if isempty(warp4template{i1,i2,nth_warp+1})
                        warp4template{i1,i2,nth_warp+1} = warp4template{i1,i2,nth_warp};
                    end
                end
            end
        end

%         if nth_warp <= n_split4warpinit/2
%             temp = cell2mat(out_parts(:,:,nth_warp-1));
%         else
%             temp = cell2mat(out_parts(:,:,nth_warp+1));
%         end
%         for i1 = 2:block_size
%             temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:) = ...
%                 nanmean(cat(3,temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:)...
%                 , temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:)),3);
%             temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:) = [];
%         end
%         for i2 = 2:block_size
%             temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end)) = ...
%                 nanmean(cat(3,temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end))...
%                 , temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1))),3);
%             temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1)) = [];
%         end
%         if nth_warp <= n_split4warpinit/2
%             template_shrink(:,:,nth_warp-1) = temp;
%         else
%             template_shrink(:,:,nth_warp+1) = temp;
%         end
    end
end

% for i =1:6
% figure,imshow(template_shrink(:,:,i),[])
% end


HalfMedianLength = 10;
for i3 = 1:nz
    for i1 = 1:block_size
        for i2 = 1:block_size
            if i3 <= HalfMedianLength
                correlated_id = rho(i1,i2,1:2*HalfMedianLength)>0.5;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,correlated_id)),3);
            elseif i3 >= nz-(HalfMedianLength-1)
                correlated_id = rho(i1,i2,end-(2*HalfMedianLength-1):end)>0.5;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,(nz-(2*HalfMedianLength-1)-1)+find(correlated_id))),3);
            else
                correlated_id = rho(i1,i2,i3-HalfMedianLength:i3+HalfMedianLength)>0.5;
                warp_cell{i1,i2,i3} = nanmedian(cell2mat(warp_cell(i1,i2,(i3-HalfMedianLength-1)+find(correlated_id))),3);  
            end
        end
    end
end

SuddenJumpThreshold = 10;
for i1 = 1:block_size
    for i2 = 1:block_size
        temp = find(sum(sum(abs(diff(cell2mat(warp_cell(i1,i2,:)),1,3)),1),2)>SuddenJumpThreshold);
        if ~isempty(temp)
            for i3 = 1:length(temp)
                warp_cell{i1,i2,temp(i3)} = NaN(2,3);
                warp_cell{i1,i2,temp(i3)+1} = NaN(2,3);
            end
        end
    end
end

for i3 = 1:nz
    for i1 = 1:block_size
        for i2 = 1:block_size
            if isempty(warp_cell{i1,i2,i3})
                warp_cell{i1,i2,i3} = NaN(2,3);
            end
        end
    end
end


for i1 = 1:block_size
    for i2 = 1:block_size
        temp = cell2mat(warp_cell(i1,i2,:));
        try
            for i3_1 = 1:2
                for i3_2 = 1:3
                    temp(i3_1,i3_2,:) = interp1(find(~isnan(temp(i3_1,i3_2,:))),permute(temp(i3_1,i3_2,~isnan(temp(i3_1,i3_2,:))),[3,1,2]),1:nz,'linear','extrap');
                end
            end
        catch
            temp = repmat(affine_init,[1,1,nz]);
        end
        for i3 = 1:nz
            warp_cell{i1,i2,i3} = temp(:,:,i3);
        end
    end
end


% apply warp
out_align = cell(block_size,block_size,nz);
for i3 = 1:nz
    for i1 = 1:block_size
        for i2 = 1:block_size
            out_align{i1,i2,i3} = spatial_interp(double(stack_align(qN_y{i1,i2},qN_x{i1,i2},i3)), warp_cell{i1,i2,i3}, 'cubic', transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
            out_align{i1,i2,i3}(out_align{i1,i2,i3}==0) = NaN;
        end
    end
end

stack_align_warp = NaN(ny,nx,nz);
for i3 = 1:nz
    temp = cell2mat(out_align(:,:,i3));
    for i1 = 2:block_size
        temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:) = ...
            nanmean(cat(3,temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:)...
            , temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:)),3);
        temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:) = [];
    end
    for i2 = 2:block_size
        temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end)) = ...
            nanmean(cat(3,temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end))...
            , temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1))),3);
        temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1)) = [];
    end
    
    temp(isnan(temp))=0;
    stack_align_warp(:,:,i3) = int16(temp); 
end



toc

% save align file and copy to server (if local)
disp('save wapred align file and warp file')
tic
write_tiff(fn_align_save,int16(stack_align_warp),info_align);
save(fn_warp2save,'warp_cell','rho');
toc 


%% generate sum movive use summary file; summary file is fast
%to load no need for local copy. 
disp('Warp and generate summed movie')
fns_summary = dir(fullfile(foldername,'*summary.mat'));
if length(fns_summary) ~= nz
   error('warp number does not match file number');
end
% if any(exist(fn_sum2save,'file')) % check if done
%     disp('Warped summed 50 has been generate, skip')
% else
% load all summary files
tic
summed_c = cell(1,1,nz);
for i=1:nz
    summary_temp=load(fullfile(foldername,fns_summary(i).name));
    summed_c{i} = summary_temp.summed(:,1+rl_RemovingPix:end-rl_RemovingPix,:); % removed n pix from left and right. 
    summed_c{i} = summed_c{i}(nonzero_row,nonzero_column,:); % remove zero-zone
end
toc
%warp from the 2nd to the last
disp('Apply warp') % NOTE: summary files is not modified by warp
tic
parfor i3 = 1:nz
     %apply warping to each frame
    for i_fr = 1:size(summed_c{i3},3)
        out_temp = cell(block_size,block_size);
        for i1 = 1:block_size
            for i2 = 1:block_size
                out_temp{i1,i2} = spatial_interp(double(summed_c{i3}(qN_y{i1,i2},qN_x{i1,i2},i_fr)), warp_cell{i1,i2,i3}, 'cubic', transform, 1:length(qN_x{i1,i2}), 1:length(qN_y{i1,i2}));
                out_temp{i1,i2}(out_temp{i1,i2}==0) = NaN;
            end
        end
        
        out_temp = cell2mat(out_temp);
        for i1 = 2:block_size
            out_temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:) = ...
                nanmean(cat(3,out_temp(qN_y{i1-1,1}(end)-(2*overlap_npix+1)+1:qN_y{i1-1,1}(end),:)...
                , out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:)),3);
            out_temp(qN_y{i1-1,1}(end)+1:qN_y{i1-1,1}(end)+(2*overlap_npix+1),:) = [];
        end
        for i2 = 2:block_size
            out_temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end)) = ...
                nanmean(cat(3,out_temp(:,qN_x{1,i2-1}(end)-(2*overlap_npix+1)+1:qN_x{1,i2-1}(end))...
                , out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1))),3);
            out_temp(:,qN_x{1,i2-1}(end)+1:qN_x{1,i2-1}(end)+(2*overlap_npix+1)) = [];
        end

        out_temp(isnan(out_temp))=0;
        
        summed_c{i3}(:,:,i_fr) = int16(out_temp); 
    end
end
summed = cell2mat(summed_c);
toc
%save
disp('save warped movie')
tic
%TiffWriter(int16(summed),fn_sum2save,16); %  code from Suite2P
write_tiff(fn_sum2save,summed);
toc
% end

%% apply warp to tiff files
fns_tiff = dir(fullfile(foldername,'*.tif'));
if length(fns_summary) ~= length(fns_tiff)
   error('tiff number does not match summary file number');
else

end
if length(fns_tiff) ~= nz
   error('warp number does not match file number');
end

%apply the warp to each file and save
parfor i = 1:length(fns_tiff)
    applywarp_Npatches(fns_tiff(i).name,fns_summary(i).name,foldername,local2savefoldername,[],ch,n_ch,warp_cell(:,:,i),transform,rl_RemovingPix,nonzero_row,nonzero_column,qN_x,qN_y,overlap_npix,block_size);
end

%% save warped sum movie to local folder if needed
% if nargin > 1
%     [~,local_sum_fn2save] = fileparts(fn_sum2save);
%     disp('save warped sum movie to local')
%     tic
%     copyfile(fn_sum2save,fullfile(local2savefoldername, [local_sum_fn2save,'.tif']));
%     toc
% end

disp('done WARP ROUTINE')
%% old

%     tkns=regexp(fns(1).name,'([^/]*)_summary.mat','once','tokens');
%     fn_root = tkns{1};
%     ffn = fullfile(save_dir,'summed',sprintf('%s_summed_%d.tif',fn_root,summary_temp.n_sum));
%     ffn_align = fullfile(save_dir,'summed',sprintf('%s_summed_align.tif',fn_root));
% need to verify warping

