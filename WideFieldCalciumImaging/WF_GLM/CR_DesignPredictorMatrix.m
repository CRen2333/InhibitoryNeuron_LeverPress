function [OutputMatrix] = CR_DesignPredictorMatrix(frame_index,n_frame,shift_start,shift_end)

OutputMatrix = zeros(n_frame,shift_start+shift_end+1);
for ff=1:shift_start+shift_end+1
    curr_frame = frame_index-shift_start-1+ff;
    curr_frame(curr_frame<1) = [];
    curr_frame(curr_frame>n_frame) = [];
    OutputMatrix(curr_frame,ff) = 1;
end