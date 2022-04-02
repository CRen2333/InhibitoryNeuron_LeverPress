function CR_PupilDiameter_VideoSegmentation_2P(Animals,Dates)%% Select file

if ischar(Animals)
    Animals = {Animals};
end
if ischar(Dates)
    Dates = {Dates};
end

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    Initial = Animal(1:2);
    for curr_date = 1:length(Dates)
        clearvars -except Animals Dates curr_animal Animal curr_date Initial
        Date = Dates{curr_date};
        if ispc
            Data_folder = fullfile('Z:\Data\ImagingRig1',Date,Animal);
        else
            Data_folder = fullfile('/usr/local/lab/Data/ImagingRig1',Date,Animal);
        end
        Sessions = dir(Data_folder);
        Sessions = {Sessions(cellfun(@(x) ~isempty(strfind(x, 'Rec')), {Sessions.name})).name};        
        for curr_session = 1:length(Sessions)
            Session = Sessions{curr_session};
            if ispc
                DataPath = fullfile('Z:\Data\ImagingRig1',Date,Animal,Session,'Pupil');
                TargetPath = fullfile('Z:\People\Chi\TwoP_IN\PupilFitting',Animal,Date,Session);
            else
                DataPath = fullfile('/usr/local/lab/Data/ImagingRig1',Date,Animal,Session,'Pupil');
                TargetPath = fullfile('/usr/local/lab/People/Chi/TwoP_IN/PupilFitting',Animal,Date,Session);
            end
            if ~exist(TargetPath)
                mkdir(TargetPath);
            end
            Existing_files = dir(TargetPath);
            Existing_files = {Existing_files(cellfun(@(x) ~isempty(strfind(x, Animal(3:end))), {Existing_files.name})).name};
            if any(cellfun(@(x) ~isempty(strfind(x, 'Baseline')), Existing_files))
                disp('Baseline has been imported');
                Baseline_done = true;
            else
                Baseline_done = false;
            end
            if any(cellfun(@(x) ~isempty(strfind(x, 'Seg')), Existing_files))
                Existing_recording_files = Existing_files(cellfun(@(x) ~isempty(strfind(x, 'Seg')), Existing_files));
                Existing_recording_files = sort(Existing_recording_files);
                for n_files = 1:length(Existing_recording_files)
                    disp([Existing_recording_files{n_files} ' Baseline has been imported']);
                end
                n_exist = length(Existing_recording_files);
            else
                n_exist = 0;
            end
            
            All_Avi = dir(DataPath);
            Dark_Avi = {All_Avi(cellfun(@(x) ~isempty(strfind(x, 'Dark')), {All_Avi.name})).name};
            Recording_Avi = {All_Avi(cellfun(@(x) ~isempty(strfind(x, '000')), {All_Avi.name})).name};
            if length(Recording_Avi) ~= 1
                warning('Check pupil recording');
                return
            end

            %% Read baseline file
            if ~Baseline_done
                disp([Animal ' ' Date ': Load ' Dark_Avi{1}]);
                tic;
                readerobj = VideoReader(fullfile(DataPath,Dark_Avi{1}));
                clear mov Mouse_Mov
                TotalFrameNum = readerobj.NumberOfFrames;
                videoWidth = readerobj.Width;
                videoHeight = readerobj.Height;
                for curr_frame = 1:TotalFrameNum
                    vidFrames = read(readerobj,curr_frame); 
                   if (mod(curr_frame,100) == 0)
                       disp(['Frame: ', num2str(curr_frame) '/' num2str(TotalFrameNum)]);
                   end
                   mov = [];
                   mov = double(vidFrames(:,:,1));
                   Mouse_mov(:,:,curr_frame) = mov;
                   clear mov
                end
                save([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_' Session '_Baseline_PupilDia_CR'],...
                    'Animal','Initial','Date','DataPath','TargetPath','Mouse_mov','-v7.3');
                toc;
            end
            
            %% Read recording movie
            disp([Animal ' ' Date ': Load ' Recording_Avi{1}]);
            tic;
            readerobj = VideoReader(fullfile(DataPath,Recording_Avi{1}));
            clear mov Mouse_Mov
            TotalFrameNum = readerobj.NumberOfFrames;
            videoWidth = readerobj.Width;
            videoHeight = readerobj.Height;
            TotalFrameNum = readerobj.NumberOfFrames;
            videoWidth = readerobj.Width;
            videoHeight = readerobj.Height;
            % 1000 frames/file
            filenum = floor(TotalFrameNum/1000);
            SelectFramePer100 = [];
            frame_count = n_exist*10;
            if exist([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_SelectFramePer100'])
                load([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_SelectFramePer100'],'SelectFramePer100','-mat');
            end
            for n = n_exist+1:filenum
                tic
                disp(['Segmentation ' num2str(n)]);
                clear mov Mouse_Mov
                for ii = 1:1000
                    curr_frame = (n-1)*1000+ii;
                    vidFrames = read(readerobj,curr_frame); 
                    if (mod(curr_frame,1000) == 0)
                        disp(['Frame: ', num2str(curr_frame) '/' num2str(TotalFrameNum)]);
                    end
                    mov = [];
                    mov = double(vidFrames(:,:,1));
                    Mouse_mov(:,:,ii) = mov;
                    if (mod(curr_frame,100) == 0)
                        frame_count = frame_count + 1;
                        SelectFramePer100(:,:,frame_count) = mov;
                    end
                end
                if n < 10
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_00' num2str(n)];
                elseif n < 100
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_0' num2str(n)];
                else
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_' num2str(n)];
                end
                save([TargetPath filesep curr_filename],'Mouse_mov','-v7.3');
                toc
                if n == 1
                    save([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_SelectFramePer100'],'SelectFramePer100','-v7.3');
                else
                    save([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_SelectFramePer100'],'SelectFramePer100','-append');
                end
            end
            if mod(TotalFrameNum,1000) == 0
                continue
            else
                disp(['Segmentation ' num2str(n+1)]);
                clear mov Mouse_mov
                tic            
                for ii = 1:TotalFrameNum-filenum*1000    
                    curr_frame = ii+filenum*1000;
                    vidFrames = read(readerobj,curr_frame); 
                    if (mod(curr_frame,1000) == 0)
                        disp(['Frame: ', num2str(curr_frame) '/' num2str(TotalFrameNum)]);
                    end
                    mov = [];
                    mov = double(vidFrames(:,:,1));
                    Mouse_mov(:,:,ii) = mov;
                    if (mod(curr_frame,100) == 0)
                           frame_count = frame_count + 1;
                           SelectFramePer100(:,:,frame_count) = mov;
                    end
                end
                if n+1 < 10
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_00' num2str(n+1)];
                elseif n+1 < 100
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_0' num2str(n+1)];
                else
                    curr_filename = [Recording_Avi{1}(1:end-4) '_Seg_' num2str(n+1)];
                end
                save([TargetPath filesep curr_filename],'Mouse_mov','-v7.3');
                toc
                clear mov
                save([TargetPath filesep Initial '_' Date '_' Animal(4:end) '_SelectFramePer100'],'SelectFramePer100','-append');
            end
            disp('Finish reading');
        end
    end
end
