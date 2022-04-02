function [fastz_flag] = CR_SplitFastz(fd_dir)
% function [fastz_flag] = SplitFastz(fd_dir)
% 
% updated 2018-5-4
%%
% fd_dir = 'Z:\Data\ImagingRig1\180417\HL262\L23_F1_Ladder';
% field_n = [1 2];
% 2018-4-18 here just work with 2 fields
% for i_f = field_n
%     mkdir(fullfile(fd_dir, ['F',num2str(i_f)]));
% end
%%
tiff_filenames = dir(fullfile(fd_dir,'*.tif'));
tiff_filenames = cat(1,{tiff_filenames(:).name});
tiff_filenames = sort(tiff_filenames);

curr_stack_filename = tiff_filenames{1};

im_info = imfinfo(fullfile(fd_dir, curr_stack_filename));
numframes = length(im_info);

[img_parameter,img_value] = strread(im_info(1).ImageDescription,'%s %s', 'delimiter','=\n');

% Confirm that stack is fast-z: if not, return
fastz_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.fastZEnable')),img_parameter,'UniformOutput',0));
fastz = str2num(img_value{fastz_indx});
if ~fastz
    warning('File is not fastz, returned');
    fastz_flag = 0;
    return;
else
    fastz_flag = 1;
end

% Get channel number
n_channel_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.channelsSave')),img_parameter,'UniformOutput',0));
n_channel = length(str2num(img_value{n_channel_indx}));

% Find number of slices
numslices_indx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'scanimage.SI4.stackNumSlices')),img_parameter,'UniformOutput',0));
numslices = str2num(img_value{numslices_indx});

for i_f = 1:numslices
    mkdir(fullfile(fd_dir, ['Z',num2str(i_f)]));
end

% Go through tiff files, parellelize 
if ispc
    parfor i_file = 1:length(tiff_filenames)
        CR_SplitFastz_single(fullfile(fd_dir, tiff_filenames{i_file}), numslices); %, fr_idx_slice, numframes);
    end
else    
    for i_file = 1:length(tiff_filenames)
        CR_SplitFastz_single(fullfile(fd_dir, tiff_filenames{i_file}), numslices); %, fr_idx_slice, numframes);
    end
end
%current recording mode is F1_green F1_red F2_green F2_red .....
% into each slice  fr1 channel 1:n, fr2 channel 1:n ...


