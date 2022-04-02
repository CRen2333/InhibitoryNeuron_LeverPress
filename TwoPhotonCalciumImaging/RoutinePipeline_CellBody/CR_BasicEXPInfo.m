%% Basic info of 2P animals
% IN = 'VIP';
% Animals = {'KP_3463808_2','KP_3475729_LR','KP_3480351_1','KP_3463808_1','CR_3526643-R'};

% IN = 'SOM';
% Animals = {'KP_3459921_1','KP_3461990_1','WL_3526549-R'};

%% KP_3463808_1
clear all;
IN = 'VIP';
Animal = 'KP_3463808_1';
Imaging_Fields{1}.Date = {'180222','180224','180226','180228','180302','180304','180306','180308','180310','180312','180314'};
Imaging_Fields{1}.Piezo = [1,1,1,1,1,1,1,1,1,0,0];
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% KP_3463808_2
clear all;
IN = 'VIP';
Animal = 'KP_3463808_2';
Imaging_Fields{1}.Date = {'180317','180319','180321','180323','180325','180327','180329','180331','180402','180404','180406'};
Imaging_Fields{1}.Piezo = zeros(1,11);
Imaging_Fields{2}.Date = {'180318','180320','180322','180324','180326','180328','180330','180401','180403','180405','180407'};
Imaging_Fields{2}.Piezo = zeros(1,11);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% KP_3475729_LR
clear all;
IN = 'VIP';
Animal = 'KP_3475729_LR';
Imaging_Fields{1}.Date = {'180316','180318','180320','180322','180324','180326','180328','180330','180401','180403','180405'};
Imaging_Fields{1}.Piezo = zeros(1,11);
Imaging_Fields{2}.Date = {'180317','180319','180321','180323','180325','180327','180329','180331','180402','180404','180406'};
Imaging_Fields{2}.Piezo = zeros(1,11);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
                
%% KP_3480351_1
clear all;
IN = 'VIP';
Animal = 'KP_3480351_1';
Imaging_Fields{1}.Date = {'180316','180318','180320','180322','180324','180326','180328','180330','180401','180403','180405'};
Imaging_Fields{1}.Piezo = zeros(1,11);
Imaging_Fields{2}.Date = {'180317','180319','180321','180323','180325','180327','180329','180331','180402','180404','180406'};
Imaging_Fields{2}.Piezo = zeros(1,11);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3526643-R
clear all;
IN = 'VIP';
Animal = 'CR_3526643-R';
Imaging_Fields{1}.Date = {'180821','180823','180825','180827','180829','180831','180902','180904','180906','180908','180910'};
Imaging_Fields{1}.Piezo = zeros(1,11);
Imaging_Fields{2}.Date = {'180822','180824','180826','180828','180830','180901','180903','180905','180907','180909','180911'};
Imaging_Fields{2}.Piezo = zeros(1,11);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526642-L
clear all;
IN = 'VIP';
Animal = 'WL_3526642-L';
Imaging_Fields{1}.Date = {'180914','180916','180918','180920','180922','180924','180926','180928','180930','181002','181004'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'180914','180916','180918','180920','180922','180924','180926','180928','180930','181002','181004'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'180915','180917','180919','180921','180923','180925','180927','180929','181001','181003','181006'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'180915','180917','180919','180921','180923','180925','180927','180929','181001','181003','181006'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526641-R
clear all;
IN = 'VIP';
Animal = 'WL_3526641-R';
Imaging_Fields{1}.Date = {'180908','180910','180912','180914','180916','180918','180920','180922','180924','180926','180928'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'180908','180910','180912','180914','180916','180918','180920','180922','180924','180926','180928'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'180909','180911','180913','180915','180917','180919','180921','180923','180925','180927','180929'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'180909','180911','180913','180915','180917','180919','180921','180923','180925','180927','180929'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526642-R
clear all;
IN = 'VIP';
Animal = 'WL_3526642-R';
Imaging_Fields{1}.Date = {'180927','180929','181001','181003','181005','181007','181009','181011','181013','181015','181017'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'180927','180929','181001','181003','181005','181007','181009','181011','181013','181015','181017'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'180928','180930','181002','181004','181006','181008','181010','181012','181014','181016','181018'};
Imaging_Fields{3}.Piezo = zeros(1,11);

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526549-R
clear all;
IN = 'SOM';
Animal = 'WL_3526549-R';
Imaging_Fields{1}.Date = {'180905','180907','180909','180911','180913','180915','180917','180919','180921','180923','180925'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'180905','180907','180909','180911','180913','180915','180917','180919','180921','180923','180925'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'180906','180908','180910','180912','180914','180916','180918','180920','180922','180924','180926'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'180906','180908','180910','180912','180914','180916','180918','180920','180922','180924','180926'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% KP_3459921_1
clear all;
IN = 'SOM';
Animal = 'KP_3459921_1';
Imaging_Fields{1}.Date = {'180223','180225','180227','180301','180303','180305','180307','180309','180312','180314','180316'};
Imaging_Fields{1}.Piezo = [1,1,1,1,1,1,1,1,0,0,0];
Imaging_Fields{2}.Date = {'180223','180225','180227','180301','180303','180305','180307','180309','180312','180314','180316'};
Imaging_Fields{2}.Piezo = [2,2,2,2,2,2,2,2,0,0,0];
Imaging_Fields{3}.Date = {'180224','180226','180228','180302','180304','180306','180308','180311','180313','180315'};
Imaging_Fields{3}.Piezo = [1,1,1,1,1,1,1,0,0,0];
Imaging_Fields{4}.Date = {'180224','180226','180228','180302','180304','180306','180308','180311','180313','180315'};
Imaging_Fields{4}.Piezo = [2,2,2,2,2,2,2,0,0,0];
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% KP_3461990_1
clear all;
IN = 'SOM';
Animal = 'KP_3461990_1';
Imaging_Fields{1}.Date = {'180223','180225','180227','180301','180303','180305','180307','180309','180313','180315'};
Imaging_Fields{1}.Piezo = [1,1,1,1,1,1,1,1,0,0];
Imaging_Fields{2}.Date = {'180223','180225','180227','180301','180303','180305','180307','180309','180313','180315'};
Imaging_Fields{2}.Piezo = [2,2,2,2,2,2,2,2,0,0];
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526580-L
clear all;
IN = 'SOM';
Animal = 'WL_3526580-L';
Imaging_Fields{1}.Date = {'181008','181010','181012','181014','181016','181018','181020','181022','181024','181026','181028'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'181008','181010','181012','181014','181016','181018','181020','181022','181024','181026','181028'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'181009','181011','181013','181015','181017','181019','181021','181023','181025','181027','181029'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'181009','181011','181013','181015','181017','181019','181021','181023','181025','181027','181029'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526580-O
clear all;
IN = 'SOM';
Animal = 'WL_3526580-O';
Imaging_Fields{1}.Date = {'181010','181012','181014','181016','181019','181021','181023','181025','181027','181029','181031'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'181010','181012','181014','181016','181019','181021','181023','181025','181027','181029','181031'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'181011','181013','181015','181018','181020','181022','181024','181026','181028','181030','181101'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'181011','181013','181015','181018','181020','181022','181024','181026','181028','181030','181101'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526578-O
clear all;
IN = 'SOM';
Animal = 'WL_3526578-O';
Imaging_Fields{1}.Date = {'181012','181014','181019','181021','181023','181025','181027','181029','181031','181102','181105'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'181012','181014','181019','181021','181023','181025','181027','181029','181031','181102','181105'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'181013','181015','181020','181022','181024','181026','181028','181030','181101','181103','181106'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'181013','181015','181020','181022','181024','181026','181028','181030','181101','181103','181106'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3526578-R
clear all;
IN = 'SOM';
Animal = 'WL_3526578-R';
Imaging_Fields{1}.Date = {'181114','181116','181118','181120','181122','181124','181126','181128','181202','181204','181207'};
Imaging_Fields{1}.Piezo = zeros(1,11);
Imaging_Fields{2}.Date = {'181115','181117','181119','181121','181123','181125','181127','181129','181203','181206','181208'};
Imaging_Fields{2}.Piezo = ones(1,11);
Imaging_Fields{3}.Date = {'181115','181117','181119','181121','181123','181125','181127','181129','181203','181206','181208'};
Imaging_Fields{3}.Piezo = ones(1,11).*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3547272-O
clear all;
IN = 'SOM';
Animal = 'WL_3547272-O';
Imaging_Fields{1}.Date = {'181212','181216','181222','181230','190101'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{2}.Date = {'181212','181216','181222','181230','190101'};
Imaging_Fields{2}.Piezo = ones(1,5)*2;
Imaging_Fields{3}.Date = {'181213','181217','181223','181231','190102'};
Imaging_Fields{3}.Piezo = ones(1,5);
Imaging_Fields{4}.Date = {'181213','181217','181223','181231','190102'};
Imaging_Fields{4}.Piezo = ones(1,5)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3547272-L
clear all;
IN = 'SOM';
Animal = 'WL_3547272-L';
Imaging_Fields{1}.Date = {'181212','181214','181216','181218','181220','181222','181224','181226','181228','181230','190101'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'181212','181214','181216','181218','181220','181222','181224','181226','181228','181230','190101'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'181213','181215','181217','181219','181221','181223','181225','181227','181229','181231','190102'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'181213','181215','181217','181219','181221','181223','181225','181227','181229','181231','190102'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3547273-R
clear all;
IN = 'SOM';
Animal = 'WL_3547273-R';
Imaging_Fields{1}.Date = {'190107','190109','190111','190113','190115','190117','190121','190123','190125'};
Imaging_Fields{1}.Piezo = ones(1,9);
Imaging_Fields{2}.Date = {'190107','190109','190111','190113','190115','190117','190121','190123','190125'};
Imaging_Fields{2}.Piezo = ones(1,9)*2;
Imaging_Fields{3}.Date = {'190108','190110','190112','190114','190116','190118','190122','190124','190126'};
Imaging_Fields{3}.Piezo = ones(1,9);
Imaging_Fields{4}.Date = {'190108','190110','190112','190114','190116','190118','190122','190124','190126'};
Imaging_Fields{4}.Piezo = ones(1,9)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% WL_3547273-R
clear all;
IN = 'SOM';
Animal = 'WL_3547273-LR';
Imaging_Fields{1}.Date = {'190107','190109','190111','190113','190115','190117','190120','190122','190124','190126','190128'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190107','190109','190111','190113','190115','190117','190120','190122','190124','190126','190128'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'190108','190110','190112','190114','190116','190118','190121','190123','190125','190127','190129'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'190108','190110','190112','190114','190116','190118','190121','190123','190125','190127','190129'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3786142-O
clear all;
IN = 'SOM';
Animal = 'CR_3786142-O';
Imaging_Fields{1}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'190820','190822','190824','190826','190828','190830','190901','190903','190905','190907','190909'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'190820','190822','190824','190826','190828','190830','190901','190903','190905','190907','190909'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');


%% CR_3786142-R
clear all;
IN = 'SOM';
Animal = 'CR_3786142-R';
Imaging_Fields{1}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3786142-L
clear all;
IN = 'SOM';
Animal = 'CR_3786142-L';
Imaging_Fields{1}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190819','190821','190823','190825','190827','190829','190831','190902','190904','190906','190908'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'190820','190822','190824','190826','190828','190830','190901','190903','190905','190907','190909'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'190820','190822','190824','190826','190828','190830','190901','190903','190905','190907','190909'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887041-L
clear all;
IN = 'SOM';
Animal = 'CR_3887041-L';
Imaging_Fields{1}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887041-R
clear all;
IN = 'SOM';
Animal = 'CR_3887041-R';
Imaging_Fields{1}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3936483-O
clear all;
IN = 'SOM';
Animal = 'CR_3936483-O';
Imaging_Fields{1}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200202','200204','200206','200208','200210','200212','200214','200216','200218','200220','200222'};
Imaging_Fields{2}.Piezo = ones(1,11)*2;
Imaging_Fields{3}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200203','200205','200207','200209','200211','200213','200215','200217','200219','200221','200223'};
Imaging_Fields{4}.Piezo = ones(1,11)*2;

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3575265-O
clear all;
IN = 'DBH';
Animal = 'CR_3575265-O';
Imaging_Fields{1}.Date = {'190217'};
Imaging_Fields{1}.Piezo = zeros(1,1);

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3575265-O
clear all;
IN = 'DBH';
Animal = 'CR_3575265-LR';
Imaging_Fields{1}.Date = {'190307','190308','190311','190312','190317','190318','190323','190324','190325','190326','190327','190328'};
Imaging_Fields{1}.Piezo = zeros(1,12);

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% ********************PHARM*************************%%
%% WL_3526642-L
clear all;
IN = 'VIP';
Animal = 'WL_3526642-L';
Imaging_Fields{1}.Date = {'190115','190116'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190115','190116'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3526642-R
clear all;
IN = 'VIP';
Animal = 'WL_3526642-R';
Imaging_Fields{1}.Date = {'190116','190117'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190116','190117'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3526578-L
clear all;
IN = 'SOM';
Animal = 'CR_3526578-L';
Imaging_Fields{1}.Date = {'190118','190119'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190118','190119'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3547273-LR
clear all;
IN = 'SOM';
Animal = 'WL_3547273-LR';
Imaging_Fields{1}.Date = {'190128','190130'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190128','190130'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3547272-O
clear all;
IN = 'SOM';
Animal = 'WL_3547272-O';
Imaging_Fields{1}.Date = {'190206','190207','190208'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190206','190207','190208'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3547272-L
clear all;
IN = 'SOM';
Animal = 'WL_3547272-L';
Imaging_Fields{1}.Date = {'190208','190209','190210'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190208','190209','190210'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3526578-O
clear all;
IN = 'SOM';
Animal = 'WL_3526578-O';
Imaging_Fields{1}.Date = {'190213','190214','190215'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190213','190214','190215'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% WL_3526578-R
clear all;
IN = 'SOM';
Animal = 'WL_3526578-R';
Imaging_Fields{1}.Date = {'190222','190223','190224'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'190222','190223','190224'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3619106-LR
clear all;
IN = 'VIP';
Animal = 'CR_3619106-LR';
Imaging_Fields{1}.Date = {'190301','190302','190303'};
Imaging_Fields{1}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3619106-R
clear all;
IN = 'VIP';
Animal = 'CR_3619106-R';
Imaging_Fields{1}.Date = {'190301','190302','190303'};
Imaging_Fields{1}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702608-O
clear all;
IN = 'VIP';
Animal = 'CR_3702608-O';
% Imaging_Fields{1}.Date = {'190322','190323'};
% Imaging_Fields{1}.Piezo = ones(1,11);
% Imaging_Fields{2}.Date = {'190322','190323'};
% Imaging_Fields{2}.Piezo = ones(1,11).*2;

% Imaging_Fields{1}.Date = {'190410','190411','190412'};
% Imaging_Fields{1}.Piezo = ones(1,3);
% Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL'};
% Imaging_Fields{2}.Date = {'190410','190411','190412'};
% Imaging_Fields{2}.Piezo = ones(1,3).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL'};

Imaging_Fields{1}.Date = {'190419','190420','190421'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190419','190420','190421'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ant','CTRL'};

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702608-LR
clear all;
IN = 'VIP';
Animal = 'CR_3702608-LR';

% Drug = 'Nai_Ant';
% Imaging_Fields{1}.Date = {'190322','190323'};
% Imaging_Fields{1}.Piezo = ones(1,11);
% Imaging_Fields{2}.Date = {'190322','190323'};
% Imaging_Fields{2}.Piezo = ones(1,11).*2;

% Drug = 'Exp_Ago';
% Imaging_Fields{1}.Date = {'190410','190411','190412'};
% Imaging_Fields{1}.Piezo = ones(1,3);
% Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL'};
% Imaging_Fields{2}.Date = {'190410','190411','190412'};
% Imaging_Fields{2}.Piezo = ones(1,3).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL'};

Drug = 'Exp_Ant';
Imaging_Fields{1}.Date = {'190419','190420','190421'};
Imaging_Fields{1}.Piezo = ones(1,3);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190419','190420','190421'};
Imaging_Fields{2}.Piezo = ones(1,3).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ant','CTRL'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658844-L
clear all;
IN = 'VIP';
Animal = 'CR_3658844-L';
% Imaging_Fields{1}.Date = {'190413','190414'};
% Imaging_Fields{1}.Piezo = ones(1,11);
% Imaging_Fields{2}.Date = {'190413','190414'};
% Imaging_Fields{2}.Piezo = ones(1,11).*2;
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190501','190502','190503','190504','190505'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190501','190502','190503','190504','190505'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658844-R
clear all;
IN = 'VIP';
Animal = 'CR_3658844-R';
% Imaging_Fields{1}.Date = {'190413','190414'};
% Imaging_Fields{1}.Piezo = ones(1,11);
% Imaging_Fields{2}.Date = {'190413','190414'};
% Imaging_Fields{2}.Piezo = ones(1,11).*2;
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190501','190502','190503','190504','190505'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190501','190502','190503','190504','190505'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658845-R
clear all;
IN = 'VIP';
Animal = 'CR_3658845-R';
% Imaging_Fields{1}.Date = {'190420','190421'};
% Imaging_Fields{1}.Piezo = ones(1,11);
% Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
% Imaging_Fields{2}.Date = {'190420','190421'};
% Imaging_Fields{2}.Piezo = ones(1,11).*2;
% Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190509','190510','190512','190513','190514'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190509','190510','190512','190513','190514'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672035-L
clear all;
IN = 'SOM';
Animal = 'CR_3672035-L';
% Imaging_Fields{1}.Date = {'190520','190521'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
% Imaging_Fields{2}.Date = {'190520','190521'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190608','190609','190610','190611','190612'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190608','190609','190610','190611','190612'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672035-R
clear all;
IN = 'SOM';
Animal = 'CR_3672035-R';
% Imaging_Fields{1}.Date = {'190520','190521'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
% Imaging_Fields{2}.Date = {'190520','190521'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190608','190609','190610','190611','190612'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190608','190609','190610','190611','190612'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672031-O
clear all;
IN = 'VIP';
Animal = 'CR_3672031-O';
% Imaging_Fields{1}.Date = {'190523','190524'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
% Imaging_Fields{2}.Date = {'190523','190524'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190610','190611','190612','190613','190614'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190610','190611','190612','190613','190614'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';

% Imaging_Fields{1}.Date = {'190410','190411','190412'};
% Imaging_Fields{1}.Piezo = ones(1,3);
% Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL'};
% Imaging_Fields{2}.Date = {'190410','190411','190412'};
% Imaging_Fields{2}.Piezo = ones(1,3).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702224-O
clear all;
IN = 'SOM';
Animal = 'CR_3702224-O';
% Imaging_Fields{1}.Date = {'190604','190605'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
% Imaging_Fields{2}.Date = {'190604','190605'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190623','190624','190625','190626','190627'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190623','190624','190625','190626','190627'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702224-L
clear all;
IN = 'SOM';
Animal = 'CR_3702224-L';
% Imaging_Fields{1}.Date = {'190606','190607'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
% Imaging_Fields{2}.Date = {'190606','190607'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190625','190626','190627','190628','190629'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190625','190626','190627','190628','190629'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702224-R
clear all;
IN = 'SOM';
Animal = 'CR_3702224-R';
% Imaging_Fields{1}.Date = {'190606','190607'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
% Imaging_Fields{2}.Date = {'190606','190607'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190625','190626','190627','190628','190629'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190625','190626','190627','190628','190629'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702224-R
clear all;
IN = 'SOM';
Animal = 'CR_3702226-O';
% Imaging_Fields{1}.Date = {'190620','190621'};
% Imaging_Fields{1}.Piezo = ones(1,2);
% Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
% Imaging_Fields{2}.Date = {'190620','190621'};
% Imaging_Fields{2}.Piezo = ones(1,2).*2;
% Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
% Drug = 'Nai_Ant';

Imaging_Fields{1}.Date = {'190709','190710','190711','190712','190713'};
Imaging_Fields{1}.Piezo = ones(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Imaging_Fields{2}.Date = {'190709','190710','190711','190712','190713'};
Imaging_Fields{2}.Piezo = ones(1,5).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'Exp_Ago_Ant';
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3633192-L
clear all;
IN = 'VIP';
Animal = 'CR_3633192-L';
Imaging_Fields{1}.Date = {'190412','190413','190414'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL'};
Drug = 'EXP_Ago';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3633192-R
clear all;
IN = 'VIP';
Animal = 'CR_3633192-R';
Imaging_Fields{1}.Date = {'190402','190403','190404'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL'};
Drug = 'EXP_Ago';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');


%% CR_3633193-L
clear all;
IN = 'VIP';
Animal = 'CR_3633193-L';
Imaging_Fields{1}.Date = {'190424','190425','190426','190427','190428'};
Imaging_Fields{1}.Piezo = zeros(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'EXP_Ago_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3633193-R
clear all;
IN = 'VIP';
Animal = 'CR_3633193-R';
Imaging_Fields{1}.Date = {'190424','190425','190426','190427','190428'};
Imaging_Fields{1}.Piezo = zeros(1,5);
Imaging_Fields{1}.SessionType = {'CTRL','Exp_Ago','CTRL','Exp_Ant','CTRL'};
Drug = 'EXP_Ago_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');


%% CR_3886982-O
clear all;
IN = 'SOM';
Animal = 'CR_3886982-O';
Imaging_Fields{1}.Date = {'191019','191020'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191019','191020'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3886982-L
clear all;
IN = 'SOM';
Animal = 'CR_3886982-L';
Imaging_Fields{1}.Date = {'191019','191020'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191019','191020'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887044-O
clear all;
IN = 'VIP';
Animal = 'CR_3887044-O';
Imaging_Fields{1}.Date = {'191026','191027'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191026','191027'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887044-L
clear all;
IN = 'VIP';
Animal = 'CR_3887044-L';
Imaging_Fields{1}.Date = {'191026','191027'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191026','191027'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3786160-LR
clear all;
IN = 'VIP';
Animal = 'CR_3786160-LR';
Imaging_Fields{1}.Date = {'191028','191029'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191028','191029'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3886982-LR
clear all;
IN = 'SOM';
Animal = 'CR_3886982-LR';
Imaging_Fields{1}.Date = {'191028','191029'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191028','191029'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3886982-R
clear all;
IN = 'SOM';
Animal = 'CR_3886982-R';
Imaging_Fields{1}.Date = {'191031','191101'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191031','191101'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887040-O
clear all;
IN = 'SOM';
Animal = 'CR_3887040-O';
Imaging_Fields{1}.Date = {'191106','191107'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191106','191107'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887040-L
clear all;
IN = 'SOM';
Animal = 'CR_3887040-L';
Imaging_Fields{1}.Date = {'191106','191107'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191106','191107'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887043-O
clear all;
IN = 'VIP';
Animal = 'CR_3887043-O';
Imaging_Fields{1}.Date = {'191108','191109'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191108','191109'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887043-L
clear all;
IN = 'VIP';
Animal = 'CR_3887043-L';
Imaging_Fields{1}.Date = {'191108','191109'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_Ant'};
Imaging_Fields{2}.Date = {'191108','191109'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_Ant'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3887041-O
clear all;
IN = 'VIP';
Animal = 'CR_3887041-O';
Imaging_Fields{1}.Date = {'191109','191110'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_Ant','CTRL'};
Imaging_Fields{2}.Date = {'191109','191110'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_Ant','CTRL'};
Drug = 'Nai_Ant';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% ********************OPTO*************************%%
%% CR_3619106-LR
clear all;
IN = 'VIP';
Animal = 'CR_3619106-LR';
Imaging_Fields{1}.Date = {'190219','190221','190223'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{2}.Date = {'190220','190225','190226'};
Imaging_Fields{2}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3619106-R
clear all;
IN = 'VIP';
Animal = 'CR_3619106-R';
Imaging_Fields{1}.Date = {'190219','190221','190223'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{2}.Date = {'190224','190225','190226'};
Imaging_Fields{2}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3633192-L
clear all;
IN = 'VIP';
Animal = 'CR_3633192-L';
Imaging_Fields{1}.Date = {'190406','190408','190410'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{2}.Date = {'190407','190409','190411'};
Imaging_Fields{2}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3633192-R
clear all;
IN = 'VIP';
Animal = 'CR_3633192-R';
Imaging_Fields{1}.Date = {'190326','190328','190330'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{2}.Date = {'190327','190329','190331'};
Imaging_Fields{2}.Piezo = zeros(1,3);
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3633193-L & R
clear all;
IN = 'VIP';
Animal = 'CR_3633193-R';
Imaging_Fields{1}.Date = {'190418','190420','190422'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{2}.Date = {'190419','190421','190423'};
Imaging_Fields{2}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
Imaging_Fields{2}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702608-O
clear all;
IN = 'VIP';
Animal = 'CR_3702608-O';
Imaging_Fields{1}.Date = {'190414','190415','190416'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3702608-LR
clear all;
IN = 'VIP';
Animal = 'CR_3702608-LR';
Imaging_Fields{1}.Date = {'190414','190415','190416'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658844-L
clear all;
IN = 'VIP';
Animal = 'CR_3658844-L';
Imaging_Fields{1}.Date = {'190507','190508','190510'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658844-R
clear all;
IN = 'VIP';
Animal = 'CR_3658844-R';
Imaging_Fields{1}.Date = {'190507','190508','190510'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3658845-R
clear all;
IN = 'VIP';
Animal = 'CR_3658845-R';
Imaging_Fields{1}.Date = {'190516','190517','190518'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672031-O
clear all;
IN = 'VIP';
Animal = 'CR_3672031-O';
Imaging_Fields{1}.Date = {'190616','190617','190618'};
Imaging_Fields{1}.Piezo = zeros(1,3);
Imaging_Fields{1}.SessionType = {'Opto','Opto','Opto'};
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Opto' filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%%*************************************CHAT*************************************
%% CR_3619073-R
clear all;
IN = 'ChAT';
Animal = 'CR_3619073-R';
Imaging_Fields{1}.Date = {'190227','190228','190303','190304','190309','190310','190315','190316','190317','190318','190319','190320'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3619073-LR
clear all;
IN = 'ChAT';
Animal = 'CR_3619073-LR';
Imaging_Fields{1}.Date = {'190227','190228','190303','190304','190309','190310','190315','190316','190317','190318','190319','190320'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3619074-O
clear all;
IN = 'ChAT';
Animal = 'CR_3619074-O';
Imaging_Fields{1}.Date = {'190321','190322','190326','190327','190331','190402','190407','190408','190409','190410','190411','190412'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3633170-L
clear all;
IN = 'ChAT';
Animal = 'CR_3633170-L';
Imaging_Fields{1}.Date = {'190411','190412','190415','190416','190421','190422','190427','190428','190429','190430','190501','190502'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672020-O
clear all;
IN = 'ChAT';
Animal = 'CR_3672020-O';
Imaging_Fields{1}.Date = {'190429','190430','190503','190504','190510','190511','190516','190517','190518','190520','190521','190522'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
%% CR_3672020-L
clear all;
IN = 'ChAT';
Animal = 'CR_3672020-L';
Imaging_Fields{1}.Date = {'190429','190430','190503','190504','190510','190511','190516','190517','190518','190520','190521','190522'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3672020-R
clear all;
IN = 'ChAT';
Animal = 'CR_3672020-R';
Imaging_Fields{1}.Date = {'190511','190512','190515','190516','190522','190523','190528','190529','190530','190531','190601','190602'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4383143-L/R
clear all;
IN = 'ChAT';
Animal = 'CR_4383143-R';
Imaging_Fields{1}.Date = {'210807','210811','210818','210823','210825','210827'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
Imaging_Fields{2}.Date = {'210808','210812','210819','210824','210826','210828'};
Imaging_Fields{2}.Piezo = zeros(size(Imaging_Fields{2}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4429262-O
clear all;
IN = 'ChAT';
Animal = 'CR_4429262-O';
Imaging_Fields{1}.Date = {'210816','210820','210827','210901','210903','210905'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
Imaging_Fields{2}.Date = {'210817','210821','210828','210902','210904','210906'};
Imaging_Fields{2}.Piezo = zeros(size(Imaging_Fields{2}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4383144-R
clear all;
IN = 'ChAT';
Animal = 'CR_4383144-R';
Imaging_Fields{1}.Date = {'210817','210821','210828','210902','210904','210906'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
Imaging_Fields{2}.Date = {'210818','210822','210829','210903','210905','210907'};
Imaging_Fields{2}.Piezo = zeros(size(Imaging_Fields{2}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4412582-R
clear all;
IN = 'ChAT';
Animal = 'CR_4412583-LR';
Imaging_Fields{1}.Date = {'210912','210916','210923','210928','210930','211002'};
Imaging_Fields{1}.Piezo = zeros(size(Imaging_Fields{1}.Date));
Imaging_Fields{2}.Date = {'210911','210915','210922','210927','210929','211001'};
Imaging_Fields{2}.Piezo = zeros(size(Imaging_Fields{2}.Date));
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');


%% CR_3974574-R
clear all;
IN = 'VIP_SOM';
Animal = 'CR_3974574-R';
Imaging_Fields{1}.Date = {'200404','200406','200408','200410','200412','200414','200416','200418','200420','200422','200424'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200404','200406','200408','200410','200412','200414','200416','200418','200420','200422','200424'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200405','200407','200409','200411','200413','200415','200417','200419','200421','200423','200425'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200405','200407','200409','200411','200413','200415','200417','200419','200421','200423','200425'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_3974573-O
clear all;
IN = 'VIP_SOM';
Animal = 'CR_3974573-O';
Imaging_Fields{1}.Date = {'200404','200406','200408','200410','200412','200414','200416','200418','200420','200422','200424'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200404','200406','200408','200410','200412','200414','200416','200418','200420','200422','200424'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200405','200407','200409','200411','200413','200415','200417','200419','200421','200423','200425'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200405','200407','200409','200411','200413','200415','200417','200419','200421','200423','200425'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017386-O
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017386-O';
Imaging_Fields{1}.Date = {'200413','200415','200417','200419','200421','200423','200425','200427','200429','200501','200503'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200413','200415','200417','200419','200421','200423','200425','200427','200429','200501','200503'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200414','200416','200418','200420','200422','200424','200426','200428','200430','200502','200504'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200414','200416','200418','200420','200422','200424','200426','200428','200430','200502','200504'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017386-L
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017386-L';
Imaging_Fields{1}.Date = {'200523','200525','200527','200529','200531','200602','200604','200606','200608','200610','200612'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200523','200525','200527','200529','200531','200602','200604','200606','200608','200610','200612'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200524','200526','200528','200530','200601','200603','200605','200607','200609','200611','200613'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200524','200526','200528','200530','200601','200603','200605','200607','200609','200611','200613'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-L
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017421-L';
Imaging_Fields{1}.Date = {'200606','200608','200610','200612','200614','200616','200618','200620','200622','200625','200627'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200606','200608','200610','200612','200614','200616','200618','200620','200622','200625','200627'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200607','200609','200611','200613','200615','200617','200619','200621','200623','200626','200628'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200607','200609','200611','200613','200615','200617','200619','200621','200623','200626','200628'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-R
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017421-R';
Imaging_Fields{1}.Date = {'200606','200608','200610','200612','200614','200616','200618','200620','200622','200625','200627'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200606','200608','200610','200612','200614','200616','200618','200620','200622','200625','200627'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200607','200609','200611','200613','200615','200617','200619','200621','200623','200626','200628'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200607','200609','200611','200613','200615','200617','200619','200621','200623','200626','200628'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042831-R
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042831-R';
Imaging_Fields{1}.Date = {'200706','200708','200710','200712'};
Imaging_Fields{1}.Piezo = ones(1,4);
Imaging_Fields{2}.Date = {'200706','200708','200710','200712'};
Imaging_Fields{2}.Piezo = ones(1,4).*2;
Imaging_Fields{3}.Date = {'200707','200709','200711','200713'};
Imaging_Fields{3}.Piezo = ones(1,4);
Imaging_Fields{4}.Date = {'200707','200709','200711','200713'};
Imaging_Fields{4}.Piezo = ones(1,4).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042831-LR
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042831-LR';
Imaging_Fields{1}.Date = {'200706','200708','200710','200712'};
Imaging_Fields{1}.Piezo = ones(1,4);
Imaging_Fields{2}.Date = {'200706','200708','200710','200712'};
Imaging_Fields{2}.Piezo = ones(1,4).*2;
Imaging_Fields{3}.Date = {'200707','200709','200711','200713'};
Imaging_Fields{3}.Piezo = ones(1,4);
Imaging_Fields{4}.Date = {'200707','200709','200711','200713'};
Imaging_Fields{4}.Piezo = ones(1,4).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042830-O
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042830-O';
Imaging_Fields{1}.Date = {'200718','200720','200722','200724'};
Imaging_Fields{1}.Piezo = ones(1,4);
Imaging_Fields{2}.Date = {'200718','200720','200722','200724'};
Imaging_Fields{2}.Piezo = ones(1,4).*2;
Imaging_Fields{3}.Date = {'200719','200721','200723','200725'};
Imaging_Fields{3}.Piezo = ones(1,4);
Imaging_Fields{4}.Date = {'200719','200721','200723','200725'};
Imaging_Fields{4}.Piezo = ones(1,4).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');


%% CR_4042830-L
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042830-L';
Imaging_Fields{1}.Date = {'200718','200720','200722','200724'};
Imaging_Fields{1}.Piezo = ones(1,4);
Imaging_Fields{2}.Date = {'200718','200720','200722','200724'};
Imaging_Fields{2}.Piezo = ones(1,4).*2;
Imaging_Fields{3}.Date = {'200719','200721','200723','200725'};
Imaging_Fields{3}.Piezo = ones(1,4);
Imaging_Fields{4}.Date = {'200719','200721','200723','200725'};
Imaging_Fields{4}.Piezo = ones(1,4).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-LR
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017421-LR';
Imaging_Fields{1}.Date = {'200615','200617','200619','200621','200623','200626','200628','200630','200702','200704','200706'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200615','200617','200619','200621','200623','200626','200628','200630','200702','200704','200706'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200616','200618','200620','200622','200625','200627','200629','200701','200703','200705','200707'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200616','200618','200620','200622','200625','200627','200629','200701','200703','200705','200707'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4017421-O
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4017421-O';
Imaging_Fields{1}.Date = {'200527','200529','200531','200602','200604','200606','200608','200610','200612','200614','200616'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200527','200529','200531','200602','200604','200606','200608','200610','200612','200614','200616'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200528','200530','200601','200603','200605','200607','200609','200611','200613','200615','200617'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200528','200530','200601','200603','200605','200607','200609','200611','200613','200615','200617'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042831-O
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042831-O';
Imaging_Fields{1}.Date = {'200620','200622','200625','200627','200629','200701','200703','200705','200707','200709','200711'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200620','200622','200625','200627','200629','200701','200703','200705','200707','200709','200711'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200621','200623','200626','200628','200630','200702','200704','200706','200708','200710','200712'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200621','200623','200626','200628','200630','200702','200704','200706','200708','200710','200712'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042831-L
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042831-L';
Imaging_Fields{1}.Date = {'200620','200622','200625','200627','200629','200701','200703','200705','200707','200709','200711'};
Imaging_Fields{1}.Piezo = ones(1,11);
Imaging_Fields{2}.Date = {'200620','200622','200625','200627','200629','200701','200703','200705','200707','200709','200711'};
Imaging_Fields{2}.Piezo = ones(1,11).*2;
Imaging_Fields{3}.Date = {'200621','200623','200626','200628','200630','200702','200704','200706','200708','200710','200712'};
Imaging_Fields{3}.Piezo = ones(1,11);
Imaging_Fields{4}.Date = {'200621','200623','200626','200628','200630','200702','200704','200706','200708','200710','200712'};
Imaging_Fields{4}.Piezo = ones(1,11).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4042830-R
clear all;
IN = 'VIP_SOM';
Animal = 'CR_4042830-R';
Imaging_Fields{1}.Date = {'200720','200722','200724','200726'};
Imaging_Fields{1}.Piezo = ones(1,4);
Imaging_Fields{2}.Date = {'200720','200722','200724','200726'};
Imaging_Fields{2}.Piezo = ones(1,4).*2;
Imaging_Fields{3}.Date = {'200721','200723','200725','200727'};
Imaging_Fields{3}.Piezo = ones(1,4);
Imaging_Fields{4}.Date = {'200721','200723','200725','200727'};
Imaging_Fields{4}.Piezo = ones(1,4).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4113793-O/L
clear all;
IN = 'VIP_hM4Di';
Animal = 'CR_4113793-L';
Imaging_Fields{1}.Date = {'200911','200912','200913','200914','200915','200916','200917','200918'};
Imaging_Fields{1}.Piezo = ones(1,8);
Imaging_Fields{2}.Date = {'200911','200912','200913','200914','200915','200916','200917','200918'};
Imaging_Fields{2}.Piezo = ones(1,8).*2;
for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4259757-R/CR_4302983-LL/CR_4302983-R
clear all;
IN = 'VIP';
Animal = 'CR_4259757-R';
Imaging_Fields{1}.Date = {'210419','210420'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_nAnt','CTRL'};
Imaging_Fields{2}.Date = {'210419','210420'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_nAnt','CTRL'};
Drug = 'Nai_nAnt';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4302984-O/R/LR
clear all;
IN = 'VIP';
Animal = 'CR_4302984-LR';
Imaging_Fields{1}.Date = {'210423','210424'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_nAnt'};
Imaging_Fields{2}.Date = {'210423','210424'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_nAnt'};
Drug = 'Nai_nAnt';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4303048-L/R/LR
clear all;
IN = 'VIP';
Animal = 'CR_4303048-L';
Imaging_Fields{1}.Date = {'210524','210525'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'Nai_mAnt','CTRL'};
Imaging_Fields{2}.Date = {'210524','210525'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'Nai_mAnt','CTRL'};
Drug = 'Nai_mAnt';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4365784-O/L
clear all;
IN = 'VIP';
Animal = 'CR_4365784-L';
Imaging_Fields{1}.Date = {'210528','210529'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_mAnt'};
Imaging_Fields{2}.Date = {'210528','210529'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_mAnt'};
Drug = 'Nai_mAnt';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');

%% CR_4365807-O/L
clear all;
IN = 'VIP';
Animal = 'CR_4365807-L';
Imaging_Fields{1}.Date = {'210530','210531'};
Imaging_Fields{1}.Piezo = ones(1,2);
Imaging_Fields{1}.SessionType = {'CTRL','Nai_mAnt'};
Imaging_Fields{2}.Date = {'210530','210531'};
Imaging_Fields{2}.Piezo = ones(1,2).*2;
Imaging_Fields{2}.SessionType = {'CTRL','Nai_mAnt'};
Drug = 'Nai_mAnt';

for jj = 1:length(Imaging_Fields)
    for ii = 1:length(Imaging_Fields{jj}.Date)
        Date = Imaging_Fields{jj}.Date{ii};
        load(['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'LeverTrace' filesep Date filesep Animal '_' Date '_EphusTraces.mat'],'Frame_times');
        for curr_session = 1:length(Frame_times)
            step = floor(size(Frame_times{curr_session},1)/2);
            temp_1 = Frame_times{curr_session}(1:step);
            temp_2 = Frame_times{curr_session}(step+1:step*2);
            interval = nanmean([temp_2-temp_1])/step;
            if Imaging_Fields{jj}.Piezo(ii) == 0
                Imaging_Fields{jj}.FrameRate(ii) = 1/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 1
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            elseif Imaging_Fields{jj}.Piezo(ii) == 2
                Imaging_Fields{jj}.FrameRate(ii) = 0.5/interval;
            end
        end
    end
end
TargetPath = ['Z:\People\Chi\TwoP_IN\' IN filesep Animal filesep 'Pharm' filesep Drug filesep 'df_f'];
if ~exist(TargetPath)
    mkdir(TargetPath)
end
save([TargetPath filesep Animal '_ImagingInfo'], 'IN', 'Animal', 'Imaging_Fields','-v7.3');
