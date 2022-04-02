function [Licking_active, Licking_bout, Licking_bout_switch] = CR_GetLickingState(xsg_data,threshold)

if isstruct(xsg_data)
    % If xsg_data is structure from looped data, pull out sample rate
    xsg_sample_rate = xsg_data.header.acquirer.acquirer.sampleRate;
    Lick_raw = xsg_data.data.acquirer.trace_3;
else
   % If not, assume sample rate is 10 kHz (not recorded in continuous mode)
   xsg_sample_rate = 10000;
   Lick_raw = xsg_data;
end

% Get licking state
Licking_xsg_resample = resample(Lick_raw,1,10);
if threshold == 3
    Licking_active = Licking_xsg_resample < threshold;
else
    Licking_active = Licking_xsg_resample > threshold;
end
Licking_active([1,end],1) = false;
% Defind licking bout
gap_allowance = 300; % 7Hz or 300ms
Licking_allowance_max = 10000; % length limit for each licking
Licking_allowance_min = 10; % length limit for each licking
Licking_bout_min = 500; % length limit for each licking bout
Licking_active_switch = diff([0;Licking_active;0]);
Licking_active_starts = find(Licking_active_switch == 1);
Licking_active_stops = find(Licking_active_switch == -1)-1;

Licking_active_times = Licking_active_stops - Licking_active_starts;
Licking_interactive_times = Licking_active_starts(2:end) - Licking_active_stops(1:end-1);

Licking_active_erase = logical((Licking_active_times > Licking_allowance_max)+(Licking_active_times < Licking_allowance_min));
for i = find(Licking_active_erase)'
    Licking_active(Licking_active_starts(i): ...
        Licking_active_stops(i)) = 0;
end

% Update
Licking_active_switch = diff([0;Licking_active;0]);
Licking_active_starts = find(Licking_active_switch == 1);
Licking_active_stops = find(Licking_active_switch == -1)-1;

% find bout
Licking_active_times = Licking_active_stops - Licking_active_starts;
Licking_interactive_times = Licking_active_starts(2:end) - Licking_active_stops(1:end-1);
Licking_active_fill = Licking_interactive_times < gap_allowance;
Licking_bout = Licking_active;
for i = find(Licking_active_fill)'
    Licking_bout(Licking_active_stops(i): ...
        Licking_active_starts(i+1)) = 1;
end

Licking_bout_switch = diff([0;Licking_bout]);
Licking_bout_starts = find(Licking_bout_switch == 1);
Licking_bout_stops = find(Licking_bout_switch == -1)-1;

% for i = find(Licking_bout_starts)'
%     if length(find(Licking_active_switch(Licking_bout_starts(i):Licking_bout_stops(i))==1)) < 3
%         Licking_bout(Licking_bout_starts(i):Licking_bout_stops(i)) = 0;
%     end
% end

Licking_bout_times = Licking_bout_stops - Licking_bout_starts;
Licking_interbouttimes = Licking_bout_starts(2:end) - Licking_bout_stops(1:end-1);

Licking_bout_erase = logical((Licking_bout_times < Licking_bout_min));
for i = find(Licking_bout_erase)'
    Licking_bout(Licking_bout_starts(i): ...
        Licking_bout_stops(i)) = 0;
end

Licking_bout_switch = diff([0;Licking_bout]);
Licking_bout_starts = find(Licking_bout_switch == 1);
Licking_bout_stops = find(Licking_bout_switch == -1)-1;

