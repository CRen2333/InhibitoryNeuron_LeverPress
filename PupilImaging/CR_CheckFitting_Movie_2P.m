clear all
close all
clc

Initial = 'WL';
Animals = {'3526642-R'};

General_Path = 'Z:\People\Chi\TwoP_IN\PupilFitting';
for curr_animal = 1:length(Animals)
    Animal = Animals{curr_animal};
    
    FittingFile_Path = [General_Path filesep Initial '_' Animal];
    Dates = dir(FittingFile_Path);
    Dates = {Dates.name};
    Dates = sort(Dates);
    Dates = Dates(3:end);
    for curr_date = [1:2] %1:length(Dates)
        Date = Dates{curr_date};
        PresetFile = [Initial '_' Date '_' Animal '_Rec_1_Fitting_Preset.mat'];   
        if ~exist([FittingFile_Path filesep filesep Date filesep 'Rec_1' filesep PresetFile])
            continue
        end
        load([FittingFile_Path filesep Date filesep 'Rec_1' filesep PresetFile],'-mat');
        for ii = 1:size(Eye_ROI_filtered,3)
            clear temp_image temp_area level
            temp_frame = Eye_ROI_filtered(:,:,ii);

            % multithreshold = k-means clusitering, faster
            level = multithresh(temp_frame,ClusterNum);% Multilevel image thresholds using Otsu's method
            temp_image = imquantize(temp_frame,level);
            % Area control
            for kk = ClusterNum+1:-1:1 % start from brightest
                if nanmean(temp_frame(temp_image(:) == kk)) < level(round(length(level)/2)) % the selected area must be bright enough
                    continue
                end
                if sum((temp_image(:) == kk).*Seed_ROI_Mask(:)) < 500 % the selected area must colocalize with seed area
                    continue
                end
                temp_area = sum(temp_image(:) == kk);
                if temp_area > 5000 % the area must be big enough to be used
                    Eye_ROI_kmeans(:,:,ii) = temp_image == kk; % create the mask of interested region and store the ROI in Eye_ROI_kmeans
                    clear temp_area
                    break
                end
            end
        end

        Eye_ROI_kmeans([1:row_range_min,row_range_max:end],:,:) = false;
        Eye_ROI_kmeans(:,[1:col_range_min,col_range_max:end],:) = false;

        for ii = 1:size(Eye_ROI_filtered,3)
            clear row col Area
            if sum(Eye_ROI_kmeans(:,:,ii)) == 0
                continue
            end
            [row,col] = find(Eye_ROI_kmeans(:,:,ii));
            temp_center = [nanmean(row),nanmean(col)];
            for kk = 1:length(row)
                temp_dist = sqrt((row(kk)-temp_center(1))^2+(col(kk)-temp_center(2))^2);
        %         if temp_dist > Baseline_Pupil_Dia
                if temp_dist > (col_range_max-col_range_min)/1.8
                    Eye_ROI_kmeans(row(kk),col(kk),ii) = false;
                end
            end
            % Fill the hole
            Eye_ROI_kmeans(:,:,ii) = imfill(Eye_ROI_kmeans(:,:,ii),'holes');

            % Conserve, get rid of eyelid shadow
            [L,n] =  bwlabel(Eye_ROI_kmeans(:,:,ii),4);
            for kk = 1:n
                Area(1,kk) = sum(L(:) == kk);
            end
            [~,index] = max(Area);
            Eye_ROI_kmeans(:,:,ii) = L == index;
        end
        
        se = strel('disk',2);
        for ii = 1:size(Eye_ROI_filtered,3)
            Eye_ROI_kmeans(:,:,ii) = imdilate(Eye_ROI_kmeans(:,:,ii),se);
        end

        break_length = 0;
        disp('Edge detecting');
        tic
        for ii = 1:size(Eye_ROI_kmeans,3)
            if sum(sum(Eye_ROI_kmeans(:,:,ii))) == 0
                Eye_ROI_Edges(:,:,ii) = false(size(Eye_ROI_kmeans(:,:,1)));
                continue
            end
            [row,col] = find(Eye_ROI_kmeans(:,:,ii));
            temp_center = round([nanmean(row),nanmean(col)]);
            Eye_ROI_Edges(:,:,ii) = edge(Eye_ROI_kmeans(:,:,ii),'canny');
            [row,col] = find(Eye_ROI_Edges(:,:,ii));
            D = pdist([row,col],'euclidean');
            minD = round(prctile(D,15));
            Eye_ROI_Edges(temp_center(1)-minD:temp_center(1)+minD,temp_center(2)-minD:temp_center(2)+minD,ii) = false;
        %     Eye_ROI_Edges(min(row_RFL):max(row_RFL),min(col_RFL):max(col_RFL),ii) = false;

            % Get rid of relection
            [Gmag,Gdir] = imgradient(Eye_ROI_kmeans(:,:,ii),'prewitt');
            x = (1:size(Eye_ROI_kmeans(:,:,ii),2)) - temp_center(2);
            y = (1:size(Eye_ROI_kmeans(:,:,ii),1)) - temp_center(1);
            [X,Y] = meshgrid(x,y);
            angleMatrix = atan2d(-Y, X);
            angleDifference = mod(Gdir-angleMatrix + 180, 360) - 180;
            angleDifference = abs(angleDifference);
            Eye_ROI_Edges(:,:,ii) = Eye_ROI_Edges(:,:,ii).*(angleDifference>140);
            Eye_ROI_Edges(:,temp_center(2)-break_length:temp_center(2)+break_length,ii) = false;
        end
        
        for ii = 1:size(Eye_ROI_Edges,3)
            temp = Eye_ROI_Edges(:,:,ii);
            temp(Unwanted_BW) = false;
            Eye_ROI_Edges(:,:,ii) = reshape(temp,size(Unwanted_BW));
        end 
        
%         for ii = 1:size(Eye_ROI_Edges,3)
%             [row,col] = find(Eye_ROI_Edges(:,:,ii));
%             D = pdist([row,col],'euclidean');
%             params.minMajorAxis = max(absMin,round(prctile(D,85)));
%         %     params.maxMajorAxis = min(330,round(prctile(D,100)));
%             params.maxMajorAxis = min(330,round((col_range_max-col_range_min)/0.9));
%             bestFits{1,ii} = ellipseDetection(Eye_ROI_Edges(:,:,ii), params);
%         end
        
        % Movie
        close all;
        % TargetPath = ['Z:\People\Ruize\PupilCheck\' IN filesep Initial '_' Animal];
        TargetPath = ['C:\Lab\Temp\PupilCheck\' filesep Initial '_' Animal];
        if ~exist(TargetPath)
            mkdir(TargetPath)
        end
        VideoName = [Initial '_' Date '_' Animal 'PreSetFittingCheck'];
        v = VideoWriter([TargetPath filesep VideoName]);
        v.FrameRate = 15; 
        open(v);
        
        figure; set(gcf,'color','k','position',[50 50 size(Eye_ROI_filtered,2) size(Eye_ROI_filtered,1)]);
        for ii = 1:size(Eye_ROI_filtered,3)
            clf;
            imagesc(Eye_ROI_filtered(:,:,ii)/255,[0 1]);
            set(gca,'position',[0,0,1,1]);
            colormap gray; axis equal; axis off;
            hold on;
            h = imagesc(~Eye_ROI_Edges(:,:,ii)); colormap gray; axis equal; axis off;
            set(h,'AlphaData',Eye_ROI_Edges(:,:,ii));
            
            h1 = ellipse(bestFits{1,ii}(1,3),bestFits{1,ii}(1,4),bestFits{1,ii}(1,5)*pi/180,bestFits{1,ii}(1,1),bestFits{1,ii}(1,2));
            set(h1,'color',[0.7 0.95 1]);
            h2 = ellipse(bestFits{1,ii}(2,3),bestFits{1,ii}(2,4),bestFits{1,ii}(2,5)*pi/180,bestFits{1,ii}(2,1),bestFits{1,ii}(2,2));
            set(h2,'color',[0.3 0.75 0.93]);
            h3 = ellipse(bestFits{1,ii}(3,3),bestFits{1,ii}(3,4),bestFits{1,ii}(3,5)*pi/180,bestFits{1,ii}(3,1),bestFits{1,ii}(3,2));
            set(h3,'color',[0 0.45 0.74]);
            
            text(10,20,[num2str(ii) '/' num2str(size(Eye_ROI_Edges,3))],'color','w');

            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        
        close(v)
        clearvars -except Initial Animals curr_animal IN Animal General_Path FittingFile_Path Dates curr_date 
    end    
end