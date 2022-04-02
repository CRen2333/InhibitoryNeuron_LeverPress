function CR_Routine_MC_Daily(rig, date_run,animal,session,signal_channels,isPiezo)

% INPUT: date_run -str, if no input, use current date
% NOTE: this is for 2 ch recording
%%
path_n = path;
if ~isempty(strfind(path_n,'Z:\People\Chi\MatLab_Coding\Code_from_Aki'))
    rmpath(genpath('Z:\People\Chi\MatLab_Coding\Code_from_Aki'));
end
if ~isempty(strfind(path_n,'C:\Lab\Code\Code_from_Aki'))
    rmpath(genpath('C:\Lab\Code\Code_from_Aki'));
end
if ~isempty(strfind(path_n,'Z:\People\Haixin\MATLAB\aki_rank_MC'))
    rmpath(genpath('Z:\People\Haixin\MATLAB\aki_rank_MC'));
end
if ~isempty(strfind(path_n,'Z:\People\Haixin\MATLAB\matlab_motion_correct-master'))
    rmpath(genpath('Z:\People\Haixin\MATLAB\matlab_motion_correct-master'));
end
if ~isempty(strfind(path_n,'Z:\People\Chi\MatLab_Coding\Code_from_Aki\matlab_motion_correct-master_axon'))
    rmpath(genpath('Z:\People\Chi\MatLab_Coding\Code_from_Aki\matlab_motion_correct-master_axon'));
end
path_n = path;
if isempty(strfind(path_n,'Z:\People\Chi\MatLab_Coding\Code_from_Aki\matlab_motion_correct-master'))
    addpath(genpath('Z:\People\Chi\MatLab_Coding\Code_from_Aki\matlab_motion_correct-master'));
end

if nargin < 1
    date_run= datestr(date,'yymmdd');
end
disp(['Go through rec files for:', date_run]);

switch rig
    case 'MOM'
        Data_fd = 'Z:\Data\ImagingRig1'; % Raw image
%         Data_fd = 'F:\Data\BackUp162\ImagingRig1';
        % Data_server = 'Z:\People\Chi\TwoP_IN\MotionCorrection'; % Target image
    case 'BScope1'
        Data_fd = 'Z:\Data\ImagingRig3'; % Raw image
%         Data_fd = 'G:\My Drive\DBH'; % Raw image
%         Data_fd = 'F:\Data\BackUp160\ImagingRig3';
    case 'BScope3'
        Data_fd = 'Z:\Data\ImagingRig5\Data'; % Raw image
%         Data_fd = 'F:\Data\BackUp160\ImagingRig3';
end
Data_server = 'F:\Data\MotionCorrection';
% Data_server = 'Z:\People\Chi\PlanarSti_2P';

if nargin < 2
    animal = dir(fullfile(Data_fd,date_run,'CR*'));
    idx2keep = find(arrayfun(@(x) x.isdir, animal));
    animal = cat(2, {animal(idx2keep).name});
end

if isempty(animal)
    disp('no animals')
    return;
end

% Split fastZ
if isPiezo
    if exist(fullfile(Data_fd, date_run, animal, session, 'Image', 'FastzSplit.txt'),'file') ~=2
        disp('Splitting fastz...')
        [fastz_flag] = CR_SplitFastz(fullfile(Data_fd, date_run, animal, session, 'Image')); %
        % generate a split file to label it is done
        fID = fopen(fullfile(Data_fd, date_run, animal, session, 'Image', 'FastzSplit.txt'),'w');
        fclose(fID);
        disp('Done!')
    else
        disp('Already split')
    end
else % To make things simple
    disp('No piezo.')
    mkdir(fullfile(Data_fd, date_run, animal, session, 'Image', 'Z1'));
    files = dir(fullfile(Data_fd, date_run, animal, session, 'Image'));
    files = {files(cellfun(@(x) ~isempty(strfind(x, 'tif')), {files.name})).name};
    for n_file = 1:length(files)
        movefile(fullfile(Data_fd, date_run, animal, session, 'Image', files{n_file}), fullfile(Data_fd, date_run, animal, session, 'Image', 'Z1'));
    end
end

data_dir_path = [];
people_dir_path = [];
ch_align = []; ch_signal = [];

temp_ses_fd = dir(fullfile(Data_fd, date_run, animal, session, 'Image', 'Z*')); % Z
idx2keep = find(arrayfun(@(x) x.isdir, temp_ses_fd));
temp_ses_fd = temp_ses_fd(idx2keep);

for i_f = 1:length(temp_ses_fd) %
    if length(signal_channels) == 1
        if signal_channels == 1 % green
            data_dir_path = cat(2,data_dir_path, ...
                {fullfile(Data_fd,date_run,animal,session,'Image',temp_ses_fd(i_f).name)});
            people_dir_path = cat(2,people_dir_path, ...
                {fullfile(Data_server,animal,date_run,session,temp_ses_fd(i_f).name, 'motioncorrected_tiff')});
            ch_align = cat(2, ch_align, 1);
            ch_signal = cat(2, ch_signal, 1);
        elseif signal_channels == 2 % red
            data_dir_path = cat(2,data_dir_path, ...
                {fullfile(Data_fd,date_run,animal,session,'Image',temp_ses_fd(i_f).name)});
            people_dir_path = cat(2,people_dir_path, ...
                {fullfile(Data_server,animal,date_run,session,['red_',temp_ses_fd(i_f).name], 'motioncorrected_tiff')});
            ch_align = cat(2, ch_align, 1);
            ch_signal = cat(2, ch_signal, 1);
        end
    elseif length(signal_channels) == 2 % red as reference
        data_dir_path = cat(2,data_dir_path, ...
                {fullfile(Data_fd,date_run,animal,session,'Image',temp_ses_fd(i_f).name)});
        people_dir_path = cat(2,people_dir_path, ...
            {fullfile(Data_server,animal,date_run,session,temp_ses_fd(i_f).name, 'motioncorrected_tiff')});
        ch_align = cat(2, ch_align, 2);
        ch_signal = cat(2, ch_signal, 1);
        data_dir_path = cat(2,data_dir_path, ...
            {fullfile(Data_fd,date_run,animal,session,'Image',temp_ses_fd(i_f).name)});
        people_dir_path = cat(2,people_dir_path, ...
            {fullfile(Data_server,animal,date_run,session,['red_',temp_ses_fd(i_f).name], 'motioncorrected_tiff')});
        ch_align = cat(2, ch_align, 2);
        ch_signal = cat(2, ch_signal, 2);
    end
end


for ii = 1:length(data_dir_path)
    fprintf('%s\n %s\n align %d signal %d\n\n\n\n', data_dir_path{ii},people_dir_path{ii},ch_align(ii),ch_signal(ii));
end


disp('RUNNING MC')
for ii = 1:length(data_dir_path)
    batch_motion_correct_dir(data_dir_path{ii},people_dir_path{ii},ch_align(ii),ch_signal(ii));%     
end
% batch_mc_single(data_dir_path{1},people_dir_path{1},ch_align(1),ch_signal(1));
