%% This new code is based on J. Reidl et al. 2013 and discussions between Chi and Hiroshi
% Farewell 2015.

%% Set Initials and Animals
Initial = 'CR';
IN = 'VIP';
Animals = {'3438544-R'};
Events = {'CueAligned', 'MovOnsetAligned','RewardAligned'};


%% Clear
clearvars -except Initial IN Animals Events
close all
clc

%% Main Body
disp('Running ICA for all sessions...');
for curr_animal = 1:length(Animals)
    
    Animal = Animals{curr_animal};
    
    Temporal_Mean_allsession = [];
    AllPixel_allsession = [];
    
    for curr_event = 1:length(Events)
        
        Event = Events{curr_event};

        %% Data folder
        % cd(['C:/Lab/Projects/WideFieldLeverPress/Data' filesep Initial '_' Animal filesep 'EventAligned_Gap500' filesep Event]) 
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500' filesep Event]) % lab supercomputer
        All_file_list = dir(cd);
        Image_file_list = {All_file_list(cellfun(@(x) ~isempty(strfind(x, Event)), {All_file_list.name})).name};
        Image_file_list = sort(Image_file_list);
        
        %% Initialization
        FrameNumPerSession.(Event) = nan(size(Image_file_list));
        
        %% Load data session by session
        Im_Session = length(Image_file_list);
        
        for curr_session = 1:length(Image_file_list)
            tic % time it
            cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500' filesep Event]) 
            disp(['Running ' Image_file_list{curr_session} '...'])
            
            % Load image
            load(Image_file_list{curr_session}, '-mat');          

            %% Get the temporal mean of the data for reconstruction
            ValidPixel = [];
            if ~isempty(strfind(Event,'Cue'))
                ValidPixel = AlignedIm_Cue.Cue_Conc;
                FrameNumPerSession.(Event)(curr_session) = size(AlignedIm_Cue.Cue_Conc,2);
            elseif ~isempty(strfind(Event,'MovOnset'))
                ValidPixel = AlignedIm_MovOnset.MovOnset_Conc;
                FrameNumPerSession.(Event)(curr_session) = size(AlignedIm_MovOnset.MovOnset_Conc,2);
            elseif ~isempty(strfind(Event,'Reward'))
                ValidPixel = AlignedIm_Reward.Reward_Conc;
                FrameNumPerSession.(Event)(curr_session) = size(AlignedIm_Reward.Reward_Conc,2);
            end

            AllPixel_allsession = horzcat(AllPixel_allsession, ValidPixel);
            toc
        end
    end
    
    Temporal_Mean_allsession = mean(AllPixel_allsession,2);
        
    %% Check the kurtosis for all concataneted frames
    Kur_S = kurtosis(AllPixel_allsession);
    Kur_T = kurtosis(AllPixel_allsession');

    disp(['Spatial Kurtosis: ' num2str(nanmean(Kur_S)) '+/-' num2str(nanstd(Kur_S))]);    
    disp(['Temporal Kurtosis: ' num2str(nanmean(Kur_T)) '+/-' num2str(nanstd(Kur_T))]);    

    %% Step two: PCA
    % About subtracting the spatial mean: Mode pattern and sequence is same, color is reversed for most modes
    % The first latent becomes much smaller
    [COEFF_all,SCORE_all,latent_all] = pca(AllPixel_allsession'); % SCORE: Column: Component, Row: Frame
    Info = cumsum(latent_all)./sum(latent_all); % How much information
    CompNums = [80]; % Choose first XX componnets
    cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    mkdir('ICA');
    cd('ICA');
    save([Initial '_' Animal '_PCA'],'COEFF_all','SCORE_all','latent_all','-v7.3');
    
    for ii = 1:length(CompNums)
        CompNum = CompNums(ii);
        disp(num2str(CompNum));
        ModePCA = COEFF_all(:,1:CompNum); % COEFF: Row: Pixel, Column: Component
        % Plot
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']) 
        cd('ICA');
        mkdir(['ICA_' num2str(CompNum)]); cd(['ICA_' num2str(CompNum)]);
        
        figure
        set(gcf,'color','k')
        for mode = 1:20
            subaxis(4,5,mode, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
            clims = [-0.04 0.04];
            image = ModePCA(:,mode);
            imagesc(reshape(image,[128 128]),clims)
            colormap jet;
            axis square
            axis off
            title([Animal ' ModePCA' num2str(mode)]);
        end
        saveas(gcf,[Animal '_ModePCA_1-20.fig']);
        close all
        
        figure
        set(gcf,'color','k')
        for mode = 21:40
            subaxis(4,5,mode-20, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
            clims = [-0.04 0.04];
            image = ModePCA(:,mode);
            imagesc(reshape(image,[128 128]),clims)
            colormap jet;
            axis square
            axis off
            title([Animal ' ModePCA' num2str(mode)]);
        end
        saveas(gcf,[Animal '_ModePCA_21-40.fig']);
        close all
        
        if CompNum >= 60
            figure
            set(gcf,'color','k')
            for mode = 41:60
                subaxis(4,5,mode-40, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-0.04 0.04];
                image = ModePCA(:,mode);
                imagesc(reshape(image,[128 128]),clims);
                colormap jet;
                axis square
                axis off
                title([Animal ' ModePCA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModePCA_41-60.fig']);
            close all
        end
        
        if CompNum >= 80
            figure
            set(gcf,'color','k')
            for mode = 61:80
                subaxis(4,5,mode-60, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-0.04 0.04];
                image = ModePCA(:,mode);
                imagesc(reshape(image,[128 128]),clims);
                colormap jet;
                axis square
                axis off
                title([Animal ' ModePCA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModePCA_61-80.fig']);
            close all
        end
        
        if CompNum >= 100
            figure
            set(gcf,'color','k')
            for mode = 81:100
                subaxis(4,5,mode-80, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-0.04 0.04];
                image = ModePCA(:,mode);
                imagesc(reshape(image,[128 128]),clims);
                colormap jet;
                axis square
                axis off
                title([Animal ' ModePCA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModePCA_81-100.fig']);
            close all
        end


        %% Step three: ICA;
        % ICA algorithm: JADE, Cardoso, 2013, entropy, exactly same as 1997 version
        % after reconstruction
        ModeICA = [];
        SCORE_ICA_all = [];
        B = jadeR(ModePCA'); % Input: Row: Mode, Column: Pixel; Get: B: Row: Independent Component(IC), Column: Component from PCA;
        ModeICA = (B*ModePCA')'; % Get: ModeICA: Column: Independent Component(IC), Row: Pixel;
        % Get amplitude of each ICA mode
        A = inv(B)'; % column: PCA, row: IC, each column: PCA project on IC;
        SCORE_ICA_all = A*SCORE_all(:,1:CompNum)'; % Raw: ICA Component, Column: Frame
%         SCORE_ICA = mat2cell(SCORE_ICA_all', FrameNumPerSession)';
%         SCORE_ICA = cellfun(@transpose,SCORE_ICA,'UniformOutput',false);

        %% Saving
        cd(['/usr/local/lab/People/Chi/WFLP_IN/' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']) 
        cd('ICA'); cd(['ICA_' num2str(CompNum)]);
        figure
        set(gcf,'color','w')
        for mode = 1:20
            subaxis(4,5,mode, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
            clims = [-3 10];
            image = ModeICA(:,mode);
            imagesc(reshape(image,[128 128]),clims)
            colormap jet;
            axis square
            axis off
            title([Animal ' ModeICA' num2str(mode)]);
        end
        saveas(gcf,[Animal '_ModeICA_1-20.fig']);
        close all

        figure
        set(gcf,'color','w')
        for mode = 21:40
            subaxis(4,5,mode-20, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
            clims = [-3 10];
            image = ModeICA(:,mode);
            imagesc(reshape(image,[128 128]),clims)
            colormap jet;
            axis square
            axis off
            title([Animal ' ModeICA' num2str(mode)]);
        end
        saveas(gcf,[Animal '_ModeICA_21-40.fig']);
        close all
        
        if CompNum >= 60
            figure
            set(gcf,'color','w')
            for mode = 41:60
                subaxis(4,5,mode-40, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-3 10];
                image = ModeICA(:,mode);
                imagesc(reshape(image,[128 128]),clims);
                colormap jet;
                axis square
                axis off
                title([Animal ' ModeICA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModeICA_41-60.fig']);
            close all
        end
        
        if CompNum >= 80
            figure
            set(gcf,'color','w')
            for mode = 61:80
                subaxis(4,5,mode-60, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-3 10];
                image = ModeICA(:,mode);
                imagesc(reshape(image,[128 128]),clims);
                colormap jet
                axis square
                axis off
                title([Animal ' ModeICA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModeICA_61-80.fig']);
            close all
        end
        
        if CompNum >= 100
            figure
            set(gcf,'color','w')
            for mode = 81:100
                subaxis(4,5,mode-80, 'Spacing', 0.04, 'Padding', 0, 'Margin', 0.03);
                clims = [-3 10];
                image = ModeICA(:,mode);
                imagesc(reshape(image,[128 128]),clims)
                axis square
                axis off
                title([Animal ' ModeICA' num2str(mode)]);
            end
            saveas(gcf,[Animal '_ModeICA_81-100.fig']);
            close all
        end

        curr_filename = [Initial '_' Animal '_ICA_AllSession'];
        disp('Saving...')
        save(curr_filename, 'FrameNumPerSession', 'Temporal_Mean_allsession', 'Info', 'latent_all', 'ModeICA', 'SCORE_ICA_all', '-v7.3');
    end
    
    disp(['Finish ' Initial '_' Animal]);
end

    



