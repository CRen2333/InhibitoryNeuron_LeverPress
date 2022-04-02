%% Get image frames aligned to behavirol events, registered to reference image (usually session 1)
Initial = 'CR';
IN = 'VIP';
Animals = {'3438544-R'};

Image_Range_Cue = [-3:6]; % 100ms before and 200ms after
Image_Range_MovOnset = [-15:60]; % 500ms before and 2000ms after
Image_Range_Reward = [-6:30]; % 200ms before and 1000ms after

% execute or not
Cue_exe = true;
Mov_exe = true;
Reward_exe = true;

% Clear
clearvars -except IN Initial Animals Image_Range_MovOnset Image_Range_Cue Image_Range_Reward Cue_exe Mov_exe Reward_exe Bregma_label
close all
clc

% Set animal name and date
for curr_animal = 1:length(Animals)    
    
    clear Bregma_label
    
    Animal = Animals{curr_animal};   
    switch Animal
        % VIP
        case '3183935-L'
            Bregma_label = 'Fake'; % 'REAL': real coordinates based on bone landmarks; 'Fake': based on estimation, only for first a few animals who didn't have bone pictures;
        case '3183970-L'
            Bregma_label = 'Fake'; % 'REAL': real bregma; 'Fake';
        case '3183970-O'
            Bregma_label = 'Fake';
        case '3183972-L'
            Bregma_label = 'Fake';
        case '3184010-O'
            Bregma_label = 'REAL';
        case '3184011-LL'
            Bregma_label = 'Fake';
        case '3184012-R'
            Bregma_label = 'Fake';
        case '3218183-O'
            Bregma_label = 'REAL'; % 'REAL': real bregma; 'Fake';
        case '3218183-R'
            Bregma_label = 'REAL';
        case '3161016-O'
            Bregma_label = 'Fake';
        case '3161018-R'
            Bregma_label = 'Fake'; 
        case '3183959-LL'
            Bregma_label = 'Fake';
        case '3183959-LR'
            Bregma_label = 'Fake';
        case '3183958-L'
            Bregma_label = 'REAL';
        case '3183958-R'
            Bregma_label = 'REAL';
        case '3218181-R'
            Bregma_label = 'REAL';
        case '3218181-L'
            Bregma_label = 'REAL';
        case '3218181-O'
            Bregma_label = 'REAL';
        case '3233232-O'
            Bregma_label = 'REAL';
        case '3233232-L'
            Bregma_label = 'REAL';
        case '3233232-R'
            Bregma_label = 'REAL';
        case '3233233-O'
            Bregma_label = 'REAL';
        case '3233233-L'
            Bregma_label = 'REAL';
        case '3233233-R'
            Bregma_label = 'REAL';
        case '3271320-L'
            Bregma_label = 'REAL';
        case '3333408-O'
            Bregma_label = 'REAL';
        case '3333408-L'
            Bregma_label = 'REAL';
        case '3333408-R'
            Bregma_label = 'REAL';
        case '3373693-O'
            Bregma_label = 'REAL';
        case '3373693-L'
            Bregma_label = 'REAL';
        case '3373693-R'
            Bregma_label = 'REAL';
        case '3373693-LR'
            Bregma_label = 'REAL';
        case '3438483-L'
            Bregma_label = 'REAL';
        case '3438500-L'
            Bregma_label = 'REAL';
        case '3438500-O'
            Bregma_label = 'REAL';
        case '3438544-O'
            Bregma_label = 'REAL';
        case '3438544-L'
            Bregma_label = 'REAL';
        case '3438544-R'
            Bregma_label = 'REAL';
        case '3495100-R'
            Bregma_label = 'REAL';
        % SOM
        case '3358883-O'
            Bregma_label = 'REAL';
        case '3358883-R'
            Bregma_label = 'REAL';
        case '3358884-L'
            Bregma_label = 'REAL';
        case '3358884-R'
            Bregma_label = 'REAL';
        case '3438521-O'
            Bregma_label = 'REAL';
        case '3438521-L'
            Bregma_label = 'REAL';
        case '3438521-R'
            Bregma_label = 'REAL';
        case '3438521-LR'
            Bregma_label = 'REAL';
        case '3453262-O'
            Bregma_label = 'REAL';
        case '3453262-L'
            Bregma_label = 'REAL';
        % PV
        case '3258531-O'
            Bregma_label = 'REAL';
        case '3258531-L'
            Bregma_label = 'REAL';
        case '3491479-L'
            Bregma_label = 'REAL';
        case '3491479-LR'
            Bregma_label = 'REAL';
        case '3491479-R'
            Bregma_label = 'REAL';
        case '3547207-LR'
            Bregma_label = 'REAL';
        % CHAT
        case '3420509-L'
            Bregma_label = 'REAL';
        case '3420509-R'
            Bregma_label = 'REAL';
        case '3358845-O'
            Bregma_label = 'REAL';
        case '3358845-LR'
            Bregma_label = 'REAL';
        % VIP_GFP
        case '4383182-O'
            Bregma_label = 'REAL';
        case '4383182-L'
            Bregma_label = 'REAL';
        case '4383183-O'
            Bregma_label = 'REAL';
    end
    
    %% Get behavior & Image frame information
    if ispc
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'MovAnalysis']);
        load([Initial '_' Animal '_FrameIndex_Info'],'IndexInfo');
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'WarpedTiff']);
        load([Initial '_' Animal '_WarpedTiff'],'tformSimilarity');
    elseif isunix
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'MovAnalysis']);
        load([Initial '_' Animal '_FrameIndex_Info'],'IndexInfo');
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'WarpedTiff']);
        load([Initial '_' Animal '_WarpedTiff'],'tformSimilarity');
    end

    %% Create folder
    if ispc
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal]);
    elseif isunix
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal]);
    end
    mkdir(cd,'EventAligned_Gap500');
    mkdir('EventAligned_Gap500','CueAligned');
    mkdir('EventAligned_Gap500','MovOnsetAligned');
    mkdir('EventAligned_Gap500','RewardAligned');

    %% Data folder
    if ispc
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'df_f']) 
    elseif isunix
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'df_f']) 
    end
    
    %% Get all imaging date folder
    All_file_list = dir(cd);
    Image_folder_list = {All_file_list(cellfun(@(x) ~isempty(strfind(x,'17'))||~isempty(strfind(x,'18'))||~isempty(strfind(x,'21')), {All_file_list.name})).name};
    Image_folder_list = sort(Image_folder_list);
    if ismember(Animal,{'3258531-O','3258531-L'})
        task_start_date = '171129';
    elseif ismember(Animal,{'3271320-L','3333408-O','3333408-L','3333408-R'})
        task_start_date = '171105';
    elseif ismember(Animal,{'3358884-L','3358884-R'})
        task_start_date = '171213';
    elseif ismember(Animal,{'3358883-O','3358883-R'})
        task_start_date = '171214';
    elseif ismember(Animal,{'3373693-O','3373693-L'})
        task_start_date = '180111';
    elseif ismember(Animal,{'3373693-LR','3373693-R'})
        task_start_date = '180117';
    elseif ismember(Animal,{'3420509-L','3420509-R'})
        task_start_date = '180206';
    else
        task_start_date = Image_folder_list{1};
    end
    
    index = find(ismember(Image_folder_list,task_start_date));
    Image_folder_list = Image_folder_list(index:end);
    Im_Session = length(Image_folder_list);
   
    % Bregma
    switch Bregma_label
        case 'REAL'
            if ispc
                % load information of valid pixels, bregma position
                load('C:\Lab\Projects\WideFieldLeverPress\Data\RefBregma.mat');
            elseif isunix
                load('/usr/local/lab/Projects/WideFieldLeverPress/Data/RefBregma.mat');
            end
        case 'Fake'
            Bregma_Ref = [64,44];
    end
    % Load Mask, align to first day
    if ismember(Animal,{'3373693-R','3373693-LR'})
        load(['180110' filesep Initial '_180110_' Animal '_01(2).coordinatePixel'], '-mat');
        load(['180110' filesep Initial '_180110_' Animal '_01(2).pixel'], '-mat');
    elseif ismember(Animal,{'3420509-L','3420509-R'})
        load(['180129' filesep Initial '_180129_' Animal '_01(2).coordinatePixel'], '-mat');
        load(['180129' filesep Initial '_180129_' Animal '_01(2).pixel'], '-mat');
    else
        load([Image_folder_list{1} filesep Initial '_' Image_folder_list{1} '_' Animal '_01(2).coordinatePixel'], '-mat');
        load([Image_folder_list{1} filesep Initial '_' Image_folder_list{1} '_' Animal '_01(2).pixel'], '-mat');
    end
    % Mask
    PixelIndex = true(16384,1);
    PixelIndex(roiPixelNum,1) = false;
    
    for curr_session = 1:min(length(IndexInfo),length(Image_folder_list))
        tic
        if ispc
            cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'df_f'])
        elseif isunix
            cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'df_f'])
        end
        disp([Initial '_' Animal ' on ' Image_folder_list{curr_session} '...'])        
        % load df_f files
        load([Image_folder_list{curr_session} filesep Initial '_' Image_folder_list{curr_session} '_' Animal '_df_f_all.mat']);
        % Check df_f_all's dim
        Check = length(size(df_f_all));
        if Check == 3
            disp(['Reformating df_f_all on ' Image_folder_list{curr_session} '...'])
            temp_all = df_f_all;
            clear df_f_all;
            for block = 1:size(temp_all,3)
                df_f_all{1,block} = temp_all{block};
            end
        end

        % Initialize
        Im_Block = length(df_f_all);
        temp_Conc_Cue = cell(size(df_f_all));
        temp_Conc_Mov = cell(size(df_f_all));
        temp_Conc_Reward = cell(size(df_f_all));
        
        for curr_block = 1:Im_Block
            if ~isempty(IndexInfo{curr_session})
                if ~isempty(IndexInfo{curr_session}{curr_block})
                    if Cue_exe
                        Temp = IndexInfo{curr_session}{curr_block};
                        Temp = Temp(Temp(:,3)==1,:); % Only look at cued
                        for curr_trial = 1:size(Temp,1)
                            Temp_Im = [];
                            curr_frame = round(Temp(curr_trial,2).*29.98)+Image_Range_Cue;
                            if curr_frame(1)<0 || curr_frame(end)>9000
                                Temp(curr_trial,4) = -inf;
                                Temp(curr_trial,6) = -inf;
                                continue
                            end
                            Temp_Im = df_f_all{1,curr_block}(:,curr_frame);
                            temp_Conc_Cue{1,curr_block} = horzcat(temp_Conc_Cue{1,curr_block},Temp_Im);
                        end
                        Temp = Temp(Temp(:,4)~=-inf,:);
                        index_responded{curr_block,1} = double(~isnan(Temp(:,4)));
                        index_rewarded{curr_block,1} = double(~isnan(Temp(:,6)));
                    end
                    if Mov_exe
                        Temp = IndexInfo{curr_session}{curr_block};
                        Temp = Temp(~isnan(Temp(:,6)),:); % Only look at rewarded, doesn't matter whether cued or not
                        for curr_trial = 1:size(Temp,1)
                            Temp_Im = [];
                            curr_frame = round(Temp(curr_trial,7).*29.98)+Image_Range_MovOnset;
                            if curr_frame(1)<0 || curr_frame(end)>9000
                                Temp(curr_trial,3) = -inf;
                                Temp(curr_trial,13) = -inf;
                                continue
                            end
                            Temp_Im = df_f_all{1,curr_block}(:,curr_frame);
                            temp_Conc_Mov{1,curr_block} = horzcat(temp_Conc_Mov{1,curr_block},Temp_Im);
                        end        
                        Temp = Temp(Temp(:,3)~=-inf,:);
                        Temp = Temp(Temp(:,13)~=-inf,:);
                        index_cued{curr_block,1} = double(~isnan(Temp(:,3)));
                        index_Catch_Mov{curr_block,1} = double(Temp(:,13)~=0);
                    end
                    if Reward_exe
                        Temp = IndexInfo{curr_session}{curr_block};
                        Temp = Temp(~isnan(Temp(:,6)),:); % Only look at rewarded, doesn't matter whether cued or not
                        for curr_trial = 1:size(Temp,1)
                            Temp_Im = [];
                            curr_frame = round(Temp(curr_trial,6).*29.98)+Image_Range_Reward;
                            if curr_frame(1)<0 || curr_frame(end)>9000
                                Temp(curr_trial,13) = -inf;
                                continue
                            end
                            Temp_Im = df_f_all{1,curr_block}(:,curr_frame);
                            temp_Conc_Reward{1,curr_block} = horzcat(temp_Conc_Reward{1,curr_block},Temp_Im);
                        end        
                        Temp = Temp(Temp(:,13)~=-inf,:);
                        index_Catch_Reward{curr_block,1} = double(Temp(:,13)~=0);
                    end
                else
                    if Cue_exe
                        temp_Conc_Cue{1,curr_block} = [];
                        index_responded{curr_block,1} = [];
                        index_rewarded{curr_block,1} = [];
                    end
                    if Mov_exe
                        temp_Conc_Mov{1,curr_block} = [];
                        index_cued{curr_block,1} = [];
                        index_Catch_Mov{curr_block,1} = [];
                    end
                    if Reward_exe
                        temp_Conc_Reward{1,curr_block} = [];
                        index_Catch_Reward{curr_block,1} = [];
                    end
                end
            end
        end
        % Save
        if Cue_exe
            AlignedIm_Cue.Cue_Conc = cell2mat(temp_Conc_Cue);                
            if ~isempty(AlignedIm_Cue.Cue_Conc)
                if curr_session ~= 1
                    AlignedIm_Cue.Cue_Conc = WarpImage(AlignedIm_Cue.Cue_Conc, 128, tformSimilarity{curr_session});
                end
                AlignedIm_Cue.Cue_Conc(PixelIndex,:) = 0;
                AlignedIm_Cue.Cue_Conc = AlignWithBregma(AlignedIm_Cue.Cue_Conc, coordinate, Bregma_Ref);
                AlignedIm_Cue.Index_Responded = logical(cell2mat(index_responded));
                AlignedIm_Cue.Index_Rewarded = logical(cell2mat(index_rewarded));
            else
                AlignedIm_Cue.Cue_Conc = [];
                AlignedIm_Cue.Index_Responded = [];
                AlignedIm_Cue.Index_Rewarded = [];
            end                
            curr_filename = [Initial '_' Image_folder_list{curr_session} '_' Animal '_CueAligned.mat'];
            if ispc
                cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500\CueAligned'])
            elseif isunix
                cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500/CueAligned'])
            end
            save(curr_filename, 'AlignedIm_Cue','Image_Range_Cue','-v7.3');
        end
        if Mov_exe
            AlignedIm_MovOnset.MovOnset_Conc = cell2mat(temp_Conc_Mov); 
            if ~isempty(AlignedIm_MovOnset.MovOnset_Conc)
                if curr_session ~= 1
                    AlignedIm_MovOnset.MovOnset_Conc = WarpImage(AlignedIm_MovOnset.MovOnset_Conc, 128, tformSimilarity{curr_session});
                end
                AlignedIm_MovOnset.MovOnset_Conc(PixelIndex,:) = 0;
                AlignedIm_MovOnset.MovOnset_Conc = AlignWithBregma(AlignedIm_MovOnset.MovOnset_Conc, coordinate, Bregma_Ref);
                AlignedIm_MovOnset.Index_Cued = logical(cell2mat(index_cued));
                AlignedIm_MovOnset.Index_Catch = logical(cell2mat(index_Catch_Mov));
            else
                AlignedIm_MovOnset.MovOnset_Conc = [];
                AlignedIm_MovOnset.Index_Cued = [];
                AlignedIm_MovOnset.Index_Catch = [];
            end   
            curr_filename = [Initial '_' Image_folder_list{curr_session} '_' Animal '_MovOnsetAligned.mat'];
            if ispc
                cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500\MovOnsetAligned'])
            elseif isunix
                cd(['/usr/local/lab/People/Chi/WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500/MovOnsetAligned'])
            end  
            save(curr_filename, 'AlignedIm_MovOnset','Image_Range_MovOnset','-v7.3');
        end        
        if Reward_exe
            AlignedIm_Reward.Reward_Conc = cell2mat(temp_Conc_Reward);
            if ~isempty(AlignedIm_Reward.Reward_Conc)
                if curr_session ~= 1
                    AlignedIm_Reward.Reward_Conc = WarpImage(AlignedIm_Reward.Reward_Conc, 128, tformSimilarity{curr_session});
                end
                AlignedIm_Reward.Reward_Conc(PixelIndex,:) = 0;
                AlignedIm_Reward.Reward_Conc = AlignWithBregma(AlignedIm_Reward.Reward_Conc, coordinate, Bregma_Ref);
                AlignedIm_Reward.Index_Catch = logical(cell2mat(index_Catch_Reward));
            else
                AlignedIm_Reward.Reward_Conc = [];
                AlignedIm_Reward.Index_Catch = [];
            end 
            curr_filename = [Initial '_' Image_folder_list{curr_session} '_' Animal '_RewardAligned.mat'];
            if ispc
                cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500/RewardAligned'])
            elseif isunix
                cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500\RewardAligned'])
            end
            save(curr_filename, 'AlignedIm_Reward','Image_Range_Reward','-v7.3');
        end     
        
        clear df_f_all temp_Conc_Cue temp_Conc_Mov temp_Conc_Mov index_responded index_rewarded index_cued index_Catch_Mov index_Catch_Reward
        clear AlignedIm_Cue AlignedIm_MovOnset AlignedIm_Reward
        toc
    end

    disp('Finish all imaging sessions! \^o^/')
end

disp('Finish all animals')

