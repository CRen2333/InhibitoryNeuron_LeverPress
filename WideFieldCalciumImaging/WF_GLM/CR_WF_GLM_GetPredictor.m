%% GLM of wide-field calcium signal
% Organize behavirol data and prepare predictors acoording to image frame
% time

clear all
close all
clc

IN = 'VIP';
Initial = 'CR';
Animals = {'3183970-L','3183970-O','3183972-L','3184010-O','3184011-LL','3184012-R','3218183-O','3438483-L','3438500-L','3438544-R','3495100-R'};

% IN = 'SOM';
% Initial = 'CR';
% Animals = {'3183958-L','3183958-R','3183959-LL','3218181-L','3218181-O','3218181-R','3438521-O','3438521-L','3438521-R','3438521-LR','3453262-O','3453262-L'};

% IN = 'PV';
% Initial = 'CR';
% Animals = {'3161016-O','3161018-R','3233232-L','3233232-O','3233232-R','3233233-O','3233233-L','3233233-R','3491479-L','3491479-LR','3547207-LR'};

for curr_animal = 1:length(Animals)
    clearvars -except IN Initial Animals curr_animal
    Animal = Animals{curr_animal};
    disp(Animal);
    General_path = fullfile('Z:\People\Chi\WFLP_IN',IN,[Initial '_' Animal]);
    load(fullfile(General_path,'MovAnalysis',[Initial '_' Animal '_FrameIndex_Info.mat']),'Imaging_Day','IndexInfo');
    load(fullfile(General_path,'MovAnalysis',[Initial '_' Animal '_Traces.mat']),'Licking_Bout','Licking_Active','Lever_Active','Lever_Force_Resample','StiTime');
    Imaging_Day = Imaging_Day(1:11);
    % Change licking_Bout to frame time
    Licking_Bout = Licking_Bout(Imaging_Day);
    Lever_Active = Lever_Active(Imaging_Day);
    Licking_Active = Licking_Active(Imaging_Day);
    Lever_Force_Resample = Lever_Force_Resample(Imaging_Day);
    StiTime.CueOn = StiTime.CueOn(Imaging_Day);
    StiTime.CueOff = StiTime.CueOff(Imaging_Day);
    StiTime.RewardOn = StiTime.RewardOn(Imaging_Day);
    StiTime.RewardOff = StiTime.RewardOff(Imaging_Day);
    StiTime.PunishOn = StiTime.PunishOn(Imaging_Day);
    StiTime.PunishOff = StiTime.PunishOff(Imaging_Day);
    for curr_session = 1:length(Imaging_Day)
        for curr_block = 1:length(IndexInfo{curr_session})
            % Initialize
            FrameIndex.LickBout{curr_session}{curr_block} = [];
            FrameIndex.LickBoutOn{curr_session}{curr_block} = [];
            FrameIndex.LickBoutOff{curr_session}{curr_block} =[];
            FrameIndex.Mvm{curr_session}{curr_block} = [];
            FrameIndex.MvmOn{curr_session}{curr_block} = [];
            FrameIndex.MvmOff{curr_session}{curr_block} = [];
            FrameIndex.Cue{curr_session}{curr_block} = [];
            FrameIndex.CueOn{curr_session}{curr_block} = [];
            FrameIndex.CueOff{curr_session}{curr_block} = [];
            FrameIndex.Reward{curr_session}{curr_block} = [];
            FrameIndex.RewardOn{curr_session}{curr_block} = [];
            FrameIndex.RewardOff{curr_session}{curr_block} = [];
            FrameIndex.Punish{curr_session}{curr_block} = [];
            FrameIndex.PunishOn{curr_session}{curr_block} = [];
            FrameIndex.PunishOff{curr_session}{curr_block} = [];
            Traces.LeverPosition{curr_session}{curr_block} = [];
            Traces.LeverSpeed{curr_session}{curr_block} = [];
            Traces.LickRate{curr_session}{curr_block} = [];
            
            FrameIndex.RwdLickBout{curr_session}{curr_block} = [];
            FrameIndex.RwdLickBoutOn{curr_session}{curr_block} = [];
            FrameIndex.RwdLickBoutOff{curr_session}{curr_block} =[];
            FrameIndex.nRwdLickBout{curr_session}{curr_block} = [];
            FrameIndex.nRwdLickBoutOn{curr_session}{curr_block} = [];
            FrameIndex.nRwdLickBoutOff{curr_session}{curr_block} =[];
            FrameIndex.Outcome{curr_session}{curr_block} = [];
            FrameIndex.Outcome_ITI{curr_session}{curr_block} = [];
            
            if isempty(Licking_Active{curr_session})
                disp([Animal ' session ' num2str(curr_session) ' block ' num2str(curr_block) ': no xsg infomation'])
                continue
            end
                      
            % LickRate
            if isempty(Licking_Active{curr_session}{curr_block*2-1})
                disp([Animal ' session ' num2str(curr_session) ' block ' num2str(curr_block) ': Empty Licking']);
            else
                temp_licking_trace = Licking_Active{curr_session}{curr_block*2-1};
                temp_licking_trace = diff([0;temp_licking_trace]);
                temp_licking_trace = temp_licking_trace == 1;
                temp_licking_trace = movsum(temp_licking_trace,500)/500*1000;
                frame_stamp = [0:33.33:300000];
                frame_stamp = round(frame_stamp(2:end));
                Traces.LickRate{curr_session}{curr_block} = temp_licking_trace(frame_stamp);
                clear temp_licking_trace frame_stamp
            end
            
            % Licking bout                        
            if isempty(Licking_Bout{curr_session}{curr_block*2-1})
                disp([Animal ' session ' num2str(curr_session) ' block ' num2str(curr_block) ': Empty Licking Bout']);
            else
                temp_licking_trace = Licking_Bout{curr_session}{curr_block*2-1};
                temp_lick_on = find(diff([0;temp_licking_trace]) == 1)/1000;
                temp_lick_off = find(diff([0;temp_licking_trace]) == -1)/1000;
                temp_lick_on = round(temp_lick_on*29.98);
                temp_lick_off = round(temp_lick_off*29.98);
                temp_lick_on(temp_lick_on<1) = [];
                temp_lick_off(temp_lick_off>9000) = [];
                FrameIndex.LickBoutOn{curr_session}{curr_block} = temp_lick_on;
                FrameIndex.LickBoutOff{curr_session}{curr_block} = temp_lick_off;
                temp_bout_trace = zeros(1,9000);
                if length(temp_lick_on)==length(temp_lick_off)
                    for ii = 1:length(temp_lick_on)
                        temp_bout_trace(temp_lick_on(ii):temp_lick_off(ii)) = 1;
                    end
                elseif length(temp_lick_on)>length(temp_lick_off)
                    for ii = 1:length(temp_lick_on)-1
                        temp_bout_trace(temp_lick_on(ii):temp_lick_off(ii)) = 1;
                    end
                    temp_bout_trace(temp_lick_on(end):end) = 1;
                elseif length(temp_lick_on)<length(temp_lick_off)
                    temp_bout_trace(1:temp_lick_off(1)) = 1;
                    for ii = 1:length(temp_lick_on)
                        temp_bout_trace(temp_lick_on(ii):temp_lick_off(ii+1)) = 1;
                    end
                end
                FrameIndex.LickBout{curr_session}{curr_block} = find(temp_bout_trace==1);
                clear temp_licking_trace temp_lick_on temp_lick_off temp_bout_trace            
            end
            
            % Lever position and speed
            if isempty(Lever_Force_Resample{curr_session}{curr_block*2-1})
                disp([Animal ' session ' num2str(curr_session) ' block ' num2str(curr_block) ': Empty Lever']);
            else
                cutoff = 100; % 100 Hz
                temp_lever_trace = Lever_Force_Resample{curr_session}{curr_block*2-1};
                butterworth_stop = 1/(1000/cutoff); % fraction of nyquist 5/500 --> (cutoff = 10 Hz)
                [b a] = butter(4, butterworth_stop,'low');
                lever_force_smooth = filtfilt(b,a,temp_lever_trace);
                lever_force_smooth(1:cutoff) = lever_force_smooth(cutoff+1);
                lever_force_smooth(end-cutoff+1:end) = lever_force_smooth(end-cutoff);
                lever_force_smooth = movmean(lever_force_smooth,round(1000/29.98));
                frame_stamp = [0:33.33:300000];
                frame_stamp = round(frame_stamp(2:end));
                Traces.LeverPosition{curr_session}{curr_block} = lever_force_smooth(frame_stamp);
                lever_speed = abs(diff([lever_force_smooth([1,frame_stamp])]));          
                Traces.LeverSpeed{curr_session}{curr_block} = lever_speed;
                clear temp_lever_trace lever_speed frame_stamp lever_force_smooth
            end
            
            % Lever on, off, and binary
            temp_mvm_trace = Lever_Active{curr_session}{curr_block*2-1};
            temp_mvm_on = find(diff([0;temp_mvm_trace]) == 1)/1000;
            temp_mvm_off = find(diff([0;temp_mvm_trace]) == -1)/1000;
            temp_mvm_on = round(temp_mvm_on*29.98);
            temp_mvm_off = round(temp_mvm_off*29.98);
            FrameIndex.MvmOn{curr_session}{curr_block} = temp_mvm_on;
            FrameIndex.MvmOff{curr_session}{curr_block} = temp_mvm_off;
            temp_mvm_trace = zeros(1,9000);
            if length(temp_mvm_on)==length(temp_mvm_off)
                for ii = 1:length(temp_mvm_on)
                    temp_mvm_trace(temp_mvm_on(ii):temp_mvm_off(ii)) = 1;
                end
            elseif length(temp_mvm_on)>length(temp_mvm_off)
                for ii = 1:length(temp_mvm_on)-1
                    temp_mvm_trace(temp_mvm_on(ii):temp_mvm_off(ii)) = 1;
                end
                temp_mvm_trace(temp_mvm_on(end):end) = 1;
            elseif length(temp_mvm_on)<length(temp_mvm_off)
                temp_mvm_trace(1:temp_mvm_off(1)) = 1;
                for ii = 1:length(temp_mvm_on)
                    temp_mvm_trace(temp_mvm_on(ii):temp_mvm_off(ii+1)) = 1;
                end
            end
            FrameIndex.Mvm{curr_session}{curr_block} = find(temp_mvm_trace==1);
            clear temp_mvm_trace temp_mvm_on temp_mvm_off
            
            % Stimuli time
            FrameIndex.CueOn{curr_session}{curr_block} = round(StiTime.CueOn{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.CueOn{curr_session}{curr_block}(FrameIndex.CueOn{curr_session}{curr_block}<1) = [];
            FrameIndex.CueOn{curr_session}{curr_block}(FrameIndex.CueOn{curr_session}{curr_block}>9000) = [];
            FrameIndex.CueOn{curr_session}{curr_block}(isnan(FrameIndex.CueOn{curr_session}{curr_block})) = [];
            FrameIndex.CueOff{curr_session}{curr_block} = round(StiTime.CueOff{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.CueOff{curr_session}{curr_block}(FrameIndex.CueOff{curr_session}{curr_block}<1) = [];
            FrameIndex.CueOff{curr_session}{curr_block}(FrameIndex.CueOff{curr_session}{curr_block}>9000) = [];
            FrameIndex.CueOff{curr_session}{curr_block}(isnan(FrameIndex.CueOff{curr_session}{curr_block})) = [];
            temp_on = FrameIndex.CueOn{curr_session}{curr_block};
            temp_off = FrameIndex.CueOff{curr_session}{curr_block};
            temp_trace = zeros(1,9000);
            if length(temp_on)==length(temp_off)
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
            elseif length(temp_on)>length(temp_off)
                for ii = 1:length(temp_on)-1
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
                temp_trace(temp_on(end):end) = 1;
            elseif length(temp_on)<length(temp_off)
                temp_trace(1:temp_off(1)) = 1;
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii+1)) = 1;
                end
            end
            FrameIndex.Cue{curr_session}{curr_block} = find(temp_trace==1);
            clear temp_trace temp_on temp_off
            
            FrameIndex.RewardOn{curr_session}{curr_block} = round(StiTime.RewardOn{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.RewardOn{curr_session}{curr_block}(FrameIndex.RewardOn{curr_session}{curr_block}<1) = [];
            FrameIndex.RewardOn{curr_session}{curr_block}(FrameIndex.RewardOn{curr_session}{curr_block}>9000) = [];
            FrameIndex.RewardOn{curr_session}{curr_block}(isnan(FrameIndex.RewardOn{curr_session}{curr_block})) = [];
            FrameIndex.RewardOff{curr_session}{curr_block} = round(StiTime.RewardOff{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.RewardOff{curr_session}{curr_block}(FrameIndex.RewardOff{curr_session}{curr_block}<1) = [];
            FrameIndex.RewardOff{curr_session}{curr_block}(FrameIndex.RewardOff{curr_session}{curr_block}>9000) = [];
            FrameIndex.RewardOff{curr_session}{curr_block}(isnan(FrameIndex.RewardOff{curr_session}{curr_block})) = [];
            temp_on = FrameIndex.RewardOn{curr_session}{curr_block};
            temp_off = FrameIndex.RewardOff{curr_session}{curr_block};
            temp_trace = zeros(1,9000);
            if length(temp_on)==length(temp_off)
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
            elseif length(temp_on)>length(temp_off)
                for ii = 1:length(temp_on)-1
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
                temp_trace(temp_on(end):end) = 1;
            elseif length(temp_on)<length(temp_off)
                temp_trace(1:temp_off(1)) = 1;
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii+1)) = 1;
                end
            end
            FrameIndex.Reward{curr_session}{curr_block} = find(temp_trace==1);
            clear temp_trace temp_on temp_off
            
            FrameIndex.PunishOn{curr_session}{curr_block} = round(StiTime.PunishOn{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.PunishOn{curr_session}{curr_block}(FrameIndex.PunishOn{curr_session}{curr_block}<1) = [];
            FrameIndex.PunishOn{curr_session}{curr_block}(FrameIndex.PunishOn{curr_session}{curr_block}>9000) = [];
            FrameIndex.PunishOn{curr_session}{curr_block}(isnan(FrameIndex.PunishOn{curr_session}{curr_block})) = [];
            FrameIndex.PunishOff{curr_session}{curr_block} = round(StiTime.PunishOff{curr_session}{curr_block*2-1}*29.98);
            FrameIndex.PunishOff{curr_session}{curr_block}(FrameIndex.PunishOff{curr_session}{curr_block}<1) = [];
            FrameIndex.PunishOff{curr_session}{curr_block}(FrameIndex.PunishOff{curr_session}{curr_block}>9000) = [];
            FrameIndex.PunishOff{curr_session}{curr_block}(isnan(FrameIndex.PunishOff{curr_session}{curr_block})) = [];
            temp_on = FrameIndex.PunishOn{curr_session}{curr_block};
            temp_off = FrameIndex.PunishOff{curr_session}{curr_block};
            temp_trace = zeros(1,9000);
            if length(temp_on)==length(temp_off)
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
            elseif length(temp_on)>length(temp_off)
                for ii = 1:length(temp_on)-1
                    temp_trace(temp_on(ii):temp_off(ii)) = 1;
                end
                temp_trace(temp_on(end):end) = 1;
            elseif length(temp_on)<length(temp_off)
                temp_trace(1:temp_off(1)) = 1;
                for ii = 1:length(temp_on)
                    temp_trace(temp_on(ii):temp_off(ii+1)) = 1;
                end
            end
            FrameIndex.Punish{curr_session}{curr_block} = find(temp_trace==1);
            clear temp_trace temp_on temp_off
            
            % Rewarded related licking or not
            if isempty(FrameIndex.LickBoutOn{curr_session}{curr_block})
                disp([Animal ' session ' num2str(curr_session) ' block ' num2str(curr_block) ': Empty Licking Bout']);
            else          
                temp_lick_on = FrameIndex.LickBoutOn{curr_session}{curr_block};
                temp_lick_off = FrameIndex.LickBoutOff{curr_session}{curr_block};
                temp_bout_trace = zeros(1,9000); 
                temp_bout_trace_2 = zeros(1,9000); 
                if length(temp_lick_on)>length(temp_lick_off)
                    temp_lick_on = temp_lick_on(1:end-1);
                elseif length(temp_lick_on)<length(temp_lick_off)
                    temp_lick_off = temp_lick_off(2:end);
                end
                for curr_lick = 1:length(temp_lick_on)
                    if sum(ismember([temp_lick_on(curr_lick):temp_lick_off(curr_lick)],FrameIndex.Reward{curr_session}{curr_block}))
                        temp_bout_trace(temp_lick_on(curr_lick):temp_lick_off(curr_lick)) = 1;
                        FrameIndex.RwdLickBoutOn{curr_session}{curr_block} = [FrameIndex.RwdLickBoutOn{curr_session}{curr_block},temp_lick_on(curr_lick)];
                        FrameIndex.RwdLickBoutOff{curr_session}{curr_block} = [FrameIndex.RwdLickBoutOff{curr_session}{curr_block},temp_lick_off(curr_lick)];
                    else
                        temp_bout_trace_2(temp_lick_on(curr_lick):temp_lick_off(curr_lick)) = 1;
                        FrameIndex.nRwdLickBoutOn{curr_session}{curr_block} = [FrameIndex.nRwdLickBoutOn{curr_session}{curr_block},temp_lick_on(curr_lick)];
                        FrameIndex.nRwdLickBoutOff{curr_session}{curr_block} = [FrameIndex.nRwdLickBoutOff{curr_session}{curr_block},temp_lick_off(curr_lick)];
                    end
                end
                FrameIndex.RwdLickBout{curr_session}{curr_block} = find(temp_bout_trace==1);
                FrameIndex.nRwdLickBout{curr_session}{curr_block} = find(temp_bout_trace_2==1);

                clear temp_licking_trace temp_lick_on temp_lick_off temp_bout_trace temp_bout_trace_2           
            end
            
            % Outcome
            temp_trace = zeros(1,9000);
            for curr_trial = 1:length(StiTime.CueOn{curr_session}{curr_block*2-1})
                if ~isnan(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial))
                    temp_trace(round(StiTime.CueOn{curr_session}{curr_block*2-1}(curr_trial)*29.98):round(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial)*29.98)) = 1;
                end
            end
            FrameIndex.Outcome{curr_session}{curr_block} = find(temp_trace==1);
            clear temp_trace
            
            % Outcome ITI
            temp_trace = zeros(1,9000);
            for curr_trial = 1:length(StiTime.RewardOn{curr_session}{curr_block*2-1})-1
                if ~isnan(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial))
                    temp_trace(round(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial)*29.98):round(8*29.98)) = 1;
                end
            end
            if ~isempty(StiTime.RewardOn{curr_session}{curr_block*2-1})
                curr_trial = length(StiTime.RewardOn{curr_session}{curr_block*2-1});
                if ~isnan(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial))
                    if round(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial)*29.98)+round(8*29.98) > 9000
                        temp_trace(round(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial)*29.98):end) = 1;
                    else
                        temp_trace(round(StiTime.RewardOn{curr_session}{curr_block*2-1}(curr_trial)*29.98):round(8*29.98)) = 1;
                    end
                end
            end
            FrameIndex.Outcome_ITI{curr_session}{curr_block} = find(temp_trace==1);
            clear temp_trace
            
            Temp = IndexInfo{curr_session}{curr_block};
            if ~isempty(Temp)
                Temp = Temp(~isnan(Temp(:,6)),:); % Only look at rewarded, doesn't matter whether cued or no
                IndexInfo_forGLM{curr_session}{curr_block,1} = Temp;
            else
                IndexInfo_forGLM{curr_session}{curr_block,1} = [];
            end
            clear Temp
    
        end
    end
    curr_filename = [General_path filesep 'EventAligned_Gap500' filesep Initial '_' Animal '_FrameIndexForGLM.mat'];
    save(curr_filename,'FrameIndex','IndexInfo_forGLM','Traces','-v7.3')
end
