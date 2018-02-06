% CODED BY : PUNEET DHEER (RF)
% DATE : 05-02-2018
% PHASE AMPLITUDE COUPLING (TORT and CANOLTY)

% INPUT:
% low_freq_data = low frequency containing data channel
% high_freq_data = high frequency containing data channel  
% phase_range =  lower frequency range of interenst for phase (for example:[8 13] in hz)
% phase_step =  increasing in Hz steps (for example:1 in hz, then 8,9,10...13)
% phase_Bandwidth = bandwidth to cover central frequency (in hz) be carefull
% amp_range = higher frequency range of interenst for amplitude (for example:[30 80] in hz)
% amp_step =  increasing in Hz steps (for example:2 in hz, then 30,32,34...80)
% amp_Bandwidth = bandwidth to cover central frequency (in hz) be carefull
% showplot = to show plot 'yes' or 'no'
% surrogates = for significance 'yes' or 'no'
% No_surrogate = Number of surrogates
% Canolty_surrogate = 'yes' or 'no'
%
% Method = 'Tort' or 'Canolty'
% OUTPUT:
% raw_MI_matrix = phase amplitude comodulogram (raw)
% surr_MI_matrix =  phase amplitude comodulogram (zscore)
%
% References:
% [1] Tort et al. 10.1073/pnas.0810524105
% [2] High Gamma Power Is Phase-Locked to Theta Oscillations in Human Neocortex
% [3] Towards a proper estimation of phase-amplitude coupling in neural oscillations
% [4] Untangling cross-frequency coupling in neuroscience


function [raw_MI_matrix,surr_MI_matrix] = Modulation_Index(srate,low_freq_data, high_freq_data, phase_range, phase_step, phase_Bandwidth, amp_range, amp_step, amp_Bandwidth, showplot, surrogates,No_surrogate,Canolty_surrogate, Method)
tic
%number of bins used for Tort
nbin = 18;

phase_length = length(phase_range(1):phase_step:phase_range(2)); %Hz steps
amp_length = length(amp_range(1):amp_step:amp_range(2)); %Hz steps

% For comodulogram plot
raw_MI_matrix = zeros(amp_length,phase_length);
surr_MI_matrix = zeros(amp_length,phase_length);

Amp_row = 1;
Phase_column = 1;

for phase_freq = phase_range(1):phase_step:phase_range(2) %with 1Hz step
    clear Phase_f1 Phase_f2
    
    Phase_f1 = phase_freq -(phase_Bandwidth/2);
    Phase_f2 = phase_freq +(phase_Bandwidth/2);
    [time_series_for_phase] = eegfilt(low_freq_data,srate,Phase_f1,[]);
    [time_series_for_phase] = eegfilt(time_series_for_phase,srate,[],Phase_f2);
    Phase_time_series=angle(hilbert(time_series_for_phase));
     
    for amp_freq = amp_range(1):amp_step:amp_range(2) %with 2Hz step
        clear Amp_f1 Amp_f2
        
        % Specifiy bandwidth
        Amp_f1 = (amp_freq -(amp_Bandwidth/2));  %choose bandwidth in Hz to capture the sufficient amount of Amplitude (be carefull)
        Amp_f2 = (amp_freq +(amp_Bandwidth/2)); %choose bandwidth in Hz to capture the sufficient amount of Amplitude (be carefull)
                
        % Filter data at amp frequency
        
         [time_series_for_amp] = eegfilt(high_freq_data,srate,Amp_f1,[]);
         [time_series_for_amp] = eegfilt(time_series_for_amp,srate,[],Amp_f2);
         Amp_time_series= abs(hilbert(time_series_for_amp));
            
        switch Method
           
            case 'Tort'
                [MI] = MI_tort(Phase_time_series,Amp_time_series,nbin);

            case 'Canolty'
                [MI] = MVL_canolty(Phase_time_series,Amp_time_series);
            
        end
        
        %Raw MI values in matrix
        raw_MI_matrix(Amp_row,Phase_column) = MI;
            
       if strcmp(surrogates, 'yes')
            
                   if strcmp(Canolty_surrogate,'yes')
                        %Canolty Method

                       numpoints = length(low_freq_data); %% number of sample points in raw signal
                       numsurrogate = No_surrogate;
                       minskip = srate; %time lag must be at least this big
                       maxskip = numpoints-srate; %time lag must be smaller than this
                       skip = ceil(numpoints.*rand(numsurrogate*3,1));
                       skip ( find ( skip > maxskip ) ) = [ ] ;
                       skip ( find ( skip < minskip ) ) = [ ] ;
                       skip = skip ( 1 : numsurrogate, 1 ) ;
                       surrogate_m = zeros(numsurrogate,1) ;

                       for s=1:numsurrogate
                           surrogate_amplitude = [ Amp_time_series( skip(s) : end )  Amp_time_series(1 : skip(s) - 1) ];
                           %         surrogate_phase = [ LF_Phase( skip(s) : end )  LF_Phase(1 : skip(s) - 1) ];

                           surrogate_m(s)=MVL_canolty(Phase_time_series,surrogate_amplitude);
                           %          surrogate_m(s)=MVL_canolty(Amp_time_series,surrogate_phase);
                          
                           %disp(numsurrogate-s)

                       end
                   
           
                       [surrogate_mean,surrogate_std] = normfit(surrogate_m);
                       Norm_MVL=(MI-surrogate_mean)/surrogate_std;  %Z score Actual output
                       
                       surr_MI_matrix(Amp_row,Phase_column) = Norm_MVL;
                       %Norm_phase=angle(mo_raw);
                       %mvl_norm=Norm_MVL*exp(1i*Norm_phase);
                        
           
           
                   else
                        for surr = 1:No_surrogate
                                                       
                            % shuffle phase time series without altering the spectral power characteristics of the individual time series
                             Phase_time_series=Phase_time_series(randperm(length(Phase_time_series)));
                            
                             switch Method
                                
                                case 'Tort'
                                    [MI] = MI_tort(Phase_time_series,Amp_time_series,nbin);
                                    
                                    
                                case 'Canolty'
                                    [MI] = MVL_canolty(Phase_time_series,Amp_time_series);
                                    
                            end
                                                        
                            MI_surr(surr) = MI;
                        
                         end
                      
            
                        %Z-scored MI values
                        MI_zscore = (raw_MI_matrix(Amp_row,Phase_column)-mean(MI_surr))/std(MI_surr);
                        surr_MI_matrix(Amp_row,Phase_column) = MI_zscore;
            
                  end
        end
       
       
        if strcmp(showplot, 'yes')
            subplot(121)
            title('Raw Phase Amplitude Coupling')
            pcolor(phase_range(1):phase_step:phase_range(2),amp_range(1):amp_step:amp_range(2),raw_MI_matrix)
            shading interp; 
            colormap(jet)
            hold on; 
            %set(gca,'FontSize',20);
            %ylabel('Amplitude (Hz)','FontSize',15)
            %xlabel('Phase (Hz)','FontSize',15)
            ylabel('Amplitude (Hz)')
            xlabel('Phase (Hz)')
            colorbar
            drawnow
            %set(gcf, 'Color', 'w');
            set(gca,'XTick',[phase_range(1):phase_step:phase_range(2)]);
            subplot(122)
            title('Zscored Phase Amplitude Coupling')
            pcolor(phase_range(1):phase_step:phase_range(2),amp_range(1):amp_step:amp_range(2),surr_MI_matrix)
            shading interp; 
            colormap(jet)
            hold on; 
            %set(gca,'FontSize',20);
            %ylabel('Amplitude (Hz)','FontSize',15)
            %xlabel('Phase (Hz)','FontSize',15)
            ylabel('Amplitude (Hz)')
            xlabel('Phase (Hz)')
            colorbar
            drawnow
            %set(gcf, 'Color', 'w');
            set(gca,'XTick',[phase_range(1):phase_step:phase_range(2)]);
        end
            

       
        
        fprintf('Phase at: %d Hz && Amplitude at: %d Hz --> MI: %d \n',phase_freq,amp_freq, raw_MI_matrix(Amp_row,Phase_column));
        Amp_row = Amp_row + 1;
         
    end

  
    Amp_row = 1; %reset to first Amplitude range
    Phase_column = Phase_column + 1;  % next Phase_range
    
       
    toc
end
surr_MI_matrix(find(surr_MI_matrix<1.96))=0;%p<0.05
