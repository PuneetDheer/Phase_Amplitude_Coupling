 function [MVL] = MVL_canolty(LF_Phase,HF_Amplitude)
                
    CV_CA_sig = HF_Amplitude.*exp(1i*LF_Phase); % composite complex-valued signal time series or vector time series
    MVL = mean(CV_CA_sig); %Average vector time series ACROSS TIME
    MVL = abs(MVL); %length of the average vector time series (raw mean vector length)  
    
 end