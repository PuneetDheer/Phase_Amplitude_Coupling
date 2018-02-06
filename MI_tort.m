 function [MI] = MI_tort(LF_Phase,HF_Amplitude,nbin)
        
    steps=((pi-(-pi))/nbin);
    bins=-pi:steps:pi;
    %bins=linspace(-pi,+pi,n_binss);
    
    for i=1:length(bins)-1
        n_bins(i,1:2)=bins(i:i+1);
    end
        
    % calculate mean amplitude in each phase
    for i=1:length(n_bins)
        selected_ampl = find(LF_Phase>= n_bins(i,1) & LF_Phase< n_bins(i,2));
        mean_selected_ampl(i,1) = sum(HF_Amplitude(selected_ampl));
    end    
        mean_selected_ampl=mean_selected_ampl/sum(mean_selected_ampl); %[0 1] PDF      
      
        
     % normalized entropy
        
        MI=-sum(mean_selected_ampl.*log(mean_selected_ampl));
        MI = (log(nbin)-MI)/log(nbin); % log(nbin) is maximum entropy (Normalized [0 1])
        
    end