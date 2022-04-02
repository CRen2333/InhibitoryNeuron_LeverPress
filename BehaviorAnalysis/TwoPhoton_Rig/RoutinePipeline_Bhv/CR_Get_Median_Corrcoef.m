function median_corrcoef = CR_Get_Median_Corrcoef(V)
% Get the median of correlation coefficient of lever traces
% V is the lever trace matrix, each column corresponds to each trace

if size(V,1) == 1 % In certain cases only one trace
    median_corrcoef = [];
end

Corr_all_lever_trace_matrix = corrcoef(V);
Template = tril(true(size(Corr_all_lever_trace_matrix)),-1);
concat_corrcoef = Corr_all_lever_trace_matrix(Template);
median_corrcoef = nanmedian(concat_corrcoef);