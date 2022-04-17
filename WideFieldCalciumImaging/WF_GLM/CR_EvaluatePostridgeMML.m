function [MSE_train,MSE_test,VarExp_train,VarExp_test,Corr_train,Corr_test] = CR_EvaluatePostridgeMML(true_trace, fit_trace, training_index)
    % trace: row: obs    
    MeanSuqare = (fit_trace-true_trace).^2;
    TSS_train = sum((true_trace(training_index,:)-nanmean(true_trace(training_index,:))).^2,1);
    TSS_test = sum((true_trace(~training_index,:)-nanmean(true_trace(~training_index,:))).^2,1);
    MSE_train = sum(MeanSuqare(training_index,:),1)./sum(training_index);
    MSE_test = sum(MeanSuqare(~training_index,:),1)./sum(~training_index);
    VarExp_train = 1-sum(MeanSuqare(training_index,:),1)./TSS_train;
    VarExp_test = 1-sum(MeanSuqare(~training_index,:),1)./TSS_test;
    tempcorr = corrcoef([true_trace(training_index,:),fit_trace(training_index,:)]);
    tempcorr = diag(tempcorr,size(true_trace,2));
    Corr_train = tempcorr.^2;
    tempcorr_2 = corrcoef([true_trace(~training_index,:),fit_trace(~training_index,:)]);
    tempcorr_2 = diag(tempcorr_2,size(true_trace,2));
    Corr_test = tempcorr_2.^2;
end