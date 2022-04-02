%% Clear
clear all
close all
clc


Initial = 'CR';
IN = 'VIP';
Animals = {'3438544-R'};

global goodplots

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    disp(Animal);
    
    CompNum = 80;
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500\ICA\ICA_' num2str(CompNum)]);
    load([Initial '_' Animal '_ICA_AllSession.mat']);
    
    for curr_figure = 1:CompNum./20
        figure(curr_figure); set(gcf,'color','w');
        for curr_mode = (curr_figure-1)*20+1:curr_figure*20
            h(curr_mode-(curr_figure-1)*20) = subplot(4,5,curr_mode-(curr_figure-1)*20);
            temp = reshape(ModeICA(:,curr_mode),[128 128]);
            imagesc(temp,[-3 10]); colormap jet; axis square; axis off;
        end
        
        selectplots(h)
        pause % Have to pause here to make sure Mode_select can take the right number
        Mode_Selected{1,curr_figure} = goodplots;
        close all
        clear h
    end
    disp('Selection done');
    Mode_Selected = logical(cell2mat(Mode_Selected));
    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    temp = ModeICA(:,Mode_Selected);
    for retained = 1:sum(Mode_Selected)
        temp(1:300,retained) = 0;
        Mode_Retained{retained} = reshape(temp(:,retained),[128 128]);
        [~,I] = max(Mode_Retained{retained}(:));
        [I_row(retained,1), I_col(retained,1)] = ind2sub(size(Mode_Retained{retained}),I);
    end
    [~,index] = sort(I_row);
    sortMode_Retained = Mode_Retained(index);
    for ii = 1:10
        if ii*(ii+1) >= sum(Mode_Selected)
            break
        end
    end    
    figure; set(gcf,'color','w');
    for curr_mode = 1:sum(Mode_Selected)
        subplot(ii,ii+1,curr_mode)
        imagesc(sortMode_Retained{curr_mode},[-3 10]); axis square; axis off; colormap jet;
    end
    pause
    savefig(gcf,[Initial '_' Animal '_ModeRetained_' num2str(CompNum)]);
    close all
    
    cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
    save([Initial '_' Animal '_RecICA_' num2str(CompNum)],'Mode_Selected','sortMode_Retained','-v7.3');
    clear Mode_Selected sortMode_Retained
end

%% Reconstruct
EXE = true;
if EXE

    clear all
    close all
    clc

    Initial = 'CR';
    IN = 'VIP';
    Animals = {'3438544-R'};

    for curr_animal = 1:length(Animals)
        Animal = Animals{curr_animal};
        disp(Animal);

        CompNum = 80;
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
        load([Initial '_' Animal '_RecICA_' num2str(CompNum)],'Mode_Selected','sortMode_Retained','-mat');
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal '\EventAligned_Gap500\ICA\ICA_' num2str(CompNum)]);
        load([Initial '_' Animal '_ICA_AllSession.mat']);

        RecICA_Sum = ModeICA(:,Mode_Selected)*SCORE_ICA_all(Mode_Selected,:) + repmat(Temporal_Mean_allsession,1,size(SCORE_ICA_all,2));
        fields = fieldnames(FrameNumPerSession);
        for field = 1:length(fields)
            temp_FrameNum = FrameNumPerSession.(fields{field});
            temp_FrameNum = temp_FrameNum(~isnan(temp_FrameNum));
            total_num = sum(temp_FrameNum);
            temp = RecICA_Sum(:,1:total_num);
            if strcmp('CueAligned',fields{field})
                RecICA_Cue = mat2cell(temp', temp_FrameNum)';
                RecICA_Cue = cellfun(@transpose,RecICA_Cue,'UniformOutput',false);
                RecICA_Sum(:,1:total_num) = [];
                clear temp
            elseif strcmp('MovOnsetAligned',fields{field})
                RecICA_Mov = mat2cell(temp', temp_FrameNum)';
                RecICA_Mov = cellfun(@transpose,RecICA_Mov,'UniformOutput',false);
                RecICA_Sum(:,1:total_num) = [];
                clear temp
            elseif strcmp('RewardAligned',fields{field})
                RecICA_Reward = mat2cell(temp', temp_FrameNum)';
                RecICA_Reward = cellfun(@transpose,RecICA_Reward,'UniformOutput',false);
                RecICA_Sum(:,1:total_num) = [];
                clear temp
            end
        end
        cd(['Z:\People\Chi\WFLP_IN\' IN filesep Initial '_' Animal filesep 'EventAligned_Gap500']);
        save([Initial '_' Animal '_RecICA_' num2str(CompNum)],'RecICA_Cue','RecICA_Mov','RecICA_Reward','-append');

    end

end    