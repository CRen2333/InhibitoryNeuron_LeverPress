%% Select file
clear all
close all
clc

cd(['Z:\Data\OtherData\WideFieldRig']);
Initial = 'CR';

IN = 'VIP';
Animals = {'3373693-O','3373693-L'};
Dates = {'180124','180125','180126','180127','180128','180129','180130','180131','180201'};

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    for curr_date = 1:length(Dates)
        SelectFramePer100 = [];
        frame_count = 0;
        Date = Dates{curr_date};

        PathName = ['Z:\Data\OtherData\WideFieldRig\' filesep Date filesep Initial '_' Animal filesep 'Pupil'];
        TargetPath = ['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'Pupil' filesep Date];
        if ~exist(TargetPath)
            mkdir(TargetPath);
        end

        cd(PathName);
        All_Avi = dir('*.avi');
        All_Avi = {All_Avi.name};
        All_Avi = All_Avi(1:end-1); % Ignore dark
        % Segementation
        for curr_avi = 1:length(All_Avi)
            FileName = All_Avi{curr_avi};
            disp(['Loading ' FileName]);
            tic;
            readerobj = VideoReader(FileName);
            toc;
            disp(['Finish loading']);
% 
            % Get basic information
            TotalFrameNum = readerobj.NumberOfFrames;
            videoWidth = readerobj.Width;
            videoHeight = readerobj.Height;

            disp('Reading all frames...');
            clear mov Mouse_Mov

            filenum = floor(TotalFrameNum/1000);
            for n = 1:filenum
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
                    curr_filename = [FileName(1:end-4) '_Seg_00' num2str(n)];
                elseif n < 100
                    curr_filename = [FileName(1:end-4) '_Seg_0' num2str(n)];
                else
                    curr_filename = [FileName(1:end-4) '_Seg_' num2str(n)];
                end
                save([TargetPath filesep curr_filename],'Mouse_mov','-v7.3');
                toc
            end
            disp(['Segmentation ' num2str(n+1)]);
            clear mov Mouse_mov
            tic
            if TotalFrameNum-filenum*1000 == 0
                continue
            end
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
                curr_filename = [FileName(1:end-4) '_Seg_00' num2str(n+1)];
            elseif n+1 < 100
                curr_filename = [FileName(1:end-4) '_Seg_0' num2str(n+1)];
            else
                curr_filename = [FileName(1:end-4) '_Seg_' num2str(n+1)];
            end
            save([TargetPath filesep curr_filename],'Mouse_mov','-v7.3');
            toc
            clear mov
            if curr_avi == 1
                save([TargetPath filesep Initial '_' Date '_' Animal '_SelectFramePer100'],'SelectFramePer100','-v7.3');
            else
                save([TargetPath filesep Initial '_' Date '_' Animal '_SelectFramePer100'],'SelectFramePer100','-append');
            end

        end
        disp('Finish reading');
    end
end

