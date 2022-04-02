%%
close all
clear all
clc

Filter = true;
FigureExe = false;
EyeROI = true;
MovExe = true;
HoleFill = true;

Initial = 'CR';
Animals = {'3633170-L'};
Dates = {'190411'};
Session = 'Rec_1';

for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    for curr_date = 1:length(Dates)
        Date = Dates{curr_date};
        clearvars -except Animals curr_animal Animal Dates curr_date Date Initial Session
        if ispc
            TargetPath = ['Z:\People\Chi\TwoP_IN\PupilFitting' filesep Initial '_' Animal filesep Date filesep Session];
        elseif isunix
            addpath(genpath('/usr/local/lab/People/Chi/WFLP_IN/WFIN_IN/Code/CR_Pupil_Imaging'));
            TargetPath = ['/usr/local/lab/People/Chi/TwoP_IN/PupilFitting' filesep Initial '_' Animal filesep Date filesep Session];
        end
        cd(TargetPath);
        % Load baseline information
        if exist([Initial '_' Date '_' Animal '_' Session '_Baseline_PupilDia_CR.mat'],'file')
            load([Initial '_' Date '_' Animal '_' Session '_Baseline_PupilDia_CR.mat'],'Eye_ROI_Mask','Baseline_Pupil_Dia','-mat');
            load([Initial '_' Date '_' Animal '_' Session '_Fitting_PreSet.mat'],...
                'Eye_ROI_filtered','corrThreshold','break_length','ClusterNum','col_range_max','col_range_min','row_range_max','row_range_min','params','Unwanted_BW','Seed_ROI_Mask','absMin','absMax','-mat');
                TemplateForBlink  = nanmedian(Eye_ROI_filtered,3);
                clear Eye_ROI_filtered
        else
            disp('No baseline file detected!');
            continue
        end

        AllFiles = dir('*.mat');
        AllFiles = sort({AllFiles.name})';
        AllFiles = AllFiles(cellfun(@(x) ~isempty(strfind(x, '_Seg_')), AllFiles));
        AllFiles = AllFiles(cellfun(@(x) isempty(strfind(x, 'Fitting')), AllFiles));
        SegFiles = sort(AllFiles);

        % Estimate Eye ROI and boundry for kmeans
        disp(['Total file numbers: ' num2str(length(SegFiles))]);
        for curr_file = 1:length(SegFiles)
            load([Initial '_' Date '_' Animal '_' Session '_Fitting_PreSet.mat'],'break_length','ClusterNum','col_range_max','col_range_min','row_range_max','row_range_min','params','Unwanted_BW','Seed_ROI_Mask','absMin','absMax','-mat');
            tic
            FileName = SegFiles{curr_file};
            disp(['Loading ' FileName]);
            load(FileName);

            [row,col] = find(Eye_ROI_Mask);
            Eye_ROI = Mouse_mov(min(row):max(row),min(col):max(col),:);

            for ii = 1:size(Eye_ROI,3)
                Eye_ROI_filtered(:,:,ii) = medfilt2(Eye_ROI(:,:,ii),[3 3],'symmetric');
            end
        
            AveFramePerSeg(:,:,curr_file) = nanmean(Eye_ROI_filtered,3);
            save([Initial '_' Date '_' Animal '_' Session '_AveFramePerSeg'],'AveFramePerSeg','-v7.3');

        %     figure; imagesc(AveFramePerSeg(:,:,curr_file));
            clear Eye_ROI

            disp('K-means clustering');
            laser_on = mean(mean(Eye_ROI_filtered,1),2);
            laser_on = max(squeeze(laser_on));
    %         for ii = 1:size(Eye_ROI_filtered,3)
    %             m(ii) = mean(mean(Eye_ROI_filtered(:,:,ii)));
    %         end
    %         threshholod = (prctile(m,1)+prctile(m,90))/2;
            Eye_ROI_kmeans = nan(size(Eye_ROI_filtered));
            Laser_syn = nan(size(Eye_ROI_filtered,3),1);
            for ii = 1:size(Eye_ROI_filtered,3)
                clear temp_frame temp_image temp_area level
                temp_frame = Eye_ROI_filtered(:,:,ii);
                Laser_syn(ii,:) = nanmean(temp_frame(Seed_ROI_Mask));
                corr2withTemplate(ii,:) = corr2(Eye_ROI_filtered(:,:,ii),TemplateForBlink);
                if mean(mean(Eye_ROI_filtered(:,:,ii))) < laser_on/2
                    disp(['Laser off ' num2str(ii)]);
                    corr2withTemplate(ii,:) = 0;
                    continue;
                end
                if corr2withTemplate(ii,:) < corrThreshold
                    disp(['EyeBlink: ' num2str(ii)]);
                    Eye_ROI_kmeans(:,:,ii) = false(size(TemplateForBlink));
                    continue
                end            
                
                % multithreshold = k-means clusitering, faster
                level = multithresh(temp_frame,ClusterNum);
                temp_image = imquantize(temp_frame,level);
                % Arear control
                for kk = ClusterNum+1:-1:1 % start from brightest
                    if nanmean(temp_frame(temp_image(:) == kk)) < level(round(length(level)/2)-1)
                        continue
                    end
                    if sum((temp_image(:) == kk).*Seed_ROI_Mask(:)) < 500
                        continue
                    end
                    temp_area = sum(temp_image(:) == kk);
                    if temp_area > 2000
                        Eye_ROI_kmeans(:,:,ii) = temp_image == kk;
                        clear temp_area
                        break
                    end
                end
            end
        

    %         figure; hold on;
    %         for ii = 0:30:30*floor(size(Eye_ROI_filtered,3)/30)
    %             subplot(7,5,ii/30+1);
    %             imagesc(Eye_ROI_kmeans(:,:,ii+1),[0 1]);
    %             colormap(gray); axis equal; axis off;
    %         end
    % %         
    %         figure; hold on;
    %         for ii = 0:30:30*floor(size(Eye_ROI_filtered,3)/30)
    %             subplot(7,5,ii/30+1);
    %             imagesc(Eye_ROI_filtered(:,:,ii+1)/255,[0 1]);
    %             colormap(gray); axis equal; axis off;
    %         end

            Eye_ROI_kmeans([1:row_range_min,row_range_max:end],:,:) = false;
            Eye_ROI_kmeans(:,[1:col_range_min,col_range_max:end],:) = false;

            % Refine k-means
            for ii = 1:size(Eye_ROI_filtered,3)
                clear row col Area
                if corr2withTemplate(ii,:) < corrThreshold
                    Eye_ROI_kmeans(:,:,ii) = false(size(TemplateForBlink));
                    continue
                end
                [row,col] = find(Eye_ROI_kmeans(:,:,ii));
                temp_center = [nanmean(row),nanmean(col)];
                for kk = 1:length(row)
                    temp_dist = sqrt((row(kk)-temp_center(1))^2+(col(kk)-temp_center(2))^2);
    %                 if temp_dist > Baseline_Pupil_Dia
                    if temp_dist > (col_range_max-col_range_min)/1.8
                        Eye_ROI_kmeans(row(kk),col(kk),ii) = false;
                    end
                end

                % Fill the hole
                Eye_ROI_kmeans(:,:,ii) = imfill(Eye_ROI_kmeans(:,:,ii),'holes');

                % Conserve, get rid of eyelid shadow
                [L,n] =  bwlabel(Eye_ROI_kmeans(:,:,ii),4);
                if n == 0
                    Eye_ROI_kmeans(:,:,ii) = false(size(TemplateForBlink));
                    continue
                end
                for kk = 1:n
                    Area(1,kk) = sum(L(:) == kk);
                end
                [~,index] = max(Area);
                Eye_ROI_kmeans(:,:,ii) = L == index;
            end

            disp('Edge detecting');
            tic
            for ii = 1:size(Eye_ROI_kmeans,3)
                if corr2withTemplate(ii,:) < corrThreshold
                    Eye_ROI_Edges(:,:,ii) = false(size(TemplateForBlink));
                    continue
                end
                if isempty(find(Eye_ROI_kmeans(:,:,ii)))
                    corr2withTemplate(ii,:) = 0;
                    continue;
                end
                [row,col] = find(Eye_ROI_kmeans(:,:,ii));
                temp_center = round([nanmean(row),nanmean(col)]);
                Eye_ROI_Edges(:,:,ii) = edge(Eye_ROI_kmeans(:,:,ii),'canny');
                [row,col] = find(Eye_ROI_Edges(:,:,ii));
                D = pdist([row,col],'euclidean');
                minD = round(prctile(D,15));
                Eye_ROI_Edges(temp_center(1)-minD:temp_center(1)+minD,temp_center(2)-minD:temp_center(2)+minD,ii) = false;
    %             Eye_ROI_Edges(min(row_RFL):max(row_RFL),min(col_RFL):max(col_RFL),ii) = false;

                % Get rid of relection
                [Gmag,Gdir] = imgradient(Eye_ROI_kmeans(:,:,ii),'prewitt');
                x = (1:size(Eye_ROI_kmeans(:,:,ii),2)) - temp_center(2);
                y = (1:size(Eye_ROI_kmeans(:,:,ii),1)) - temp_center(1);
                [X,Y] = meshgrid(x,y);
                angleMatrix = atan2d(-Y, X);
                angleDifference = mod(Gdir-angleMatrix + 180, 360) - 180;
                angleDifference = abs(angleDifference);
                Eye_ROI_Edges(:,:,ii) = Eye_ROI_Edges(:,:,ii).*(angleDifference>140);
    %           
                % ignore bottom and up edge
                Eye_ROI_Edges(:,temp_center(2)-break_length:temp_center(2)+break_length,ii) = false;
            end
            toc
            for ii = 1:size(Eye_ROI_Edges,3)
                temp = Eye_ROI_Edges(:,:,ii);
                temp(Unwanted_BW) = false;
                Eye_ROI_Edges(:,:,ii) = reshape(temp,size(Unwanted_BW));
            end   

    %         figure; hold on;
    %         for ii = 0:30:30*floor(size(Eye_ROI_filtered,3)/30)
    %             subplot(7,5,ii/30+1);
    %             imagesc(Eye_ROI_Edges(:,:,ii+1),[0 1]);
    %             colormap(gray); axis equal; axis off;
    %         end
        
            disp('Fitting...');
            if ~exist('params','var')
                disp('No params loaded');
                break
            end
            tic
            for ii = 1:size(Eye_ROI_Edges,3)
                if corr2withTemplate(ii,:) < corrThreshold
                    bestFits{1,ii} = [];
                    continue
                end
                [row,col] = find(Eye_ROI_Edges(:,:,ii));
                D = pdist([row,col],'euclidean');
                params.minMajorAxis = max(absMin,round(prctile(D,85)));
                params.maxMajorAxis = min(absMax,round((col_range_max-col_range_min)/0.9));
                bestFits{1,ii} = ellipseDetection(Eye_ROI_Edges(:,:,ii), params);
            end
            toc

    %         for ii = 1:30:30*floor(size(Eye_ROI_Edges,3)/30)+1
    %             figure;
    %             imagesc(Eye_ROI_filtered(:,:,ii)/255,[0 1]); colormap gray; axis equal; axis off;
    %             hold on;
    %             h = imagesc(Eye_ROI_Edges(:,:,ii)); colormap gray; axis equal; axis off;
    %             set(h,'AlphaData',Eye_ROI_Edges(:,:,ii));
    %             ellipse(bestFits{1,ii}(1,3),bestFits{1,ii}(1,4),bestFits{1,ii}(1,5)*pi/180,bestFits{1,ii}(1,1),bestFits{1,ii}(1,2),'r');
    %             ellipse(bestFits{1,ii}(2,3),bestFits{1,ii}(2,4),bestFits{1,ii}(2,5)*pi/180,bestFits{1,ii}(2,1),bestFits{1,ii}(2,2),'y');
    %             ellipse(bestFits{1,ii}(3,3),bestFits{1,ii}(3,4),bestFits{1,ii}(3,5)*pi/180,bestFits{1,ii}(3,1),bestFits{1,ii}(3,2),'w');
    %             pause;
    %         end
            close all;
        
            % Refine
            clear temp_ref
            temp = [];
            for ii = 1:10
                if ~isempty(bestFits{1,ii})
                    temp = [temp; bestFits{1,ii}(1,:)];
                end
            end
            temp_ref(1,:) = median(temp);
            clear temp;
            for ii = 1:size(Eye_ROI_Edges,3)
                if isempty(bestFits{1,ii})
                    bestFits_refine{ii,1} = [];
                    PupilDia_refine(ii,1) = nan;
                    if ii == 1
                        temp_ref(ii,:) = temp_ref(1,:);
                    else
                        temp_ref(ii,:) = temp_ref(ii-1,:);
                    end
                    continue
                end
                temp_diff_1 = abs(bestFits{1,ii}(1,:)-temp_ref(end,:));
                if temp_diff_1(3) < 2 && sqrt(temp_diff_1(1)^2+temp_diff_1(2)^2) < 10
                    bestFits_refine{ii,1} = bestFits{1,ii}(1,:);
                    PupilDia_refine(ii,1) = bestFits{1,ii}(1,3);
                    temp_ref(ii,:) = bestFits{1,ii}(1,:);
                else
                    temp_diff_2 = abs(bestFits{1,ii}(2,:)-temp_ref(end,:));
                    temp_diff_3 = abs(bestFits{1,ii}(3,:)-temp_ref(end,:));
                    temp_diff_matrix = [temp_diff_1;temp_diff_2;temp_diff_3];
                    [temp_diff_min,Imin] = sort(temp_diff_matrix(:,3),'ascend');
                    temp = nan(1,6);
                    if temp_diff_min(1) < 2
                        if sqrt(temp_diff_matrix(Imin(1),1)^2+temp_diff_matrix(Imin(1),2)^2) < 10
                            bestFits_refine{ii,1} = bestFits{1,ii}(Imin(1),:);
                            PupilDia_refine(ii,1) = bestFits{1,ii}(Imin(1),3);
                            temp = bestFits{1,ii}(Imin(1),:);
                            temp_ref(ii,:) = bestFits{1,ii}(Imin(1),:);
                        elseif temp_diff_min(2) < 2
                                if sqrt(temp_diff_matrix(Imin(2),1)^2+temp_diff_matrix(Imin(2),2)^2) < 10
                                    bestFits_refine{ii,1} = bestFits{1,ii}(Imin(2),:);
                                    PupilDia_refine(ii,1) = bestFits{1,ii}(Imin(2),3);
                                    temp = bestFits{1,ii}(Imin(2),:);
                                    temp_ref(ii,:) = bestFits{1,ii}(Imin(2),:);
                                elseif temp_diff_min(3) < 2
                                    if sqrt(temp_diff_matrix(Imin(3),1)^2+temp_diff_matrix(Imin(3),2)^2) < 10
                                        bestFits_refine{ii,1} = bestFits{1,ii}(Imin(3),:);
                                        PupilDia_refine(ii,1) = bestFits{1,ii}(Imin(3),3);
                                        temp = bestFits{1,ii}(Imin(3),:);
                                        temp_ref(ii,:) = bestFits{1,ii}(Imin(3),:);
                                    end
                                end
                        end
                    end            
                    if isnan(temp(1,1))
                        bestFits_refine{ii,1} = bestFits{1,ii}(Imin(1),:);
                        PupilDia_refine(ii,1) = bestFits{1,ii}(Imin(1),3);
                        if ii > 1
                            temp_ref(ii,:) = temp_ref(ii-1,:);
                        end
                    end
                end
            end

            cd(TargetPath);
            save([FileName(1:end-4) '_Fitting'],'Animal','Initial','Date','TargetPath','Eye_ROI_filtered','corr2withTemplate',...
                    'ClusterNum','bestFits','params','absMin','absMax','bestFits_refine','PupilDia_refine','Laser_syn','-v7.3');

            clear Eye_ROI_filtered Eye_ROI_kmeans Eye_ROI_Edges bestFits bestFits_refine PupilDia_refine params
        end
    end
end