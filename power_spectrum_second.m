clear all
close all
clc

%change

path='D:\Monash\Semester 2\Internship\BlueLightData\';
% path='C:\Users\Dan\Dropbox\Monash\PhD\Blue Light\Data\';

subject_folder = {'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10','BL11','BL12','BL13','BL14','BL15','BL16','BL17','BL18','BL19','BL20','BL21','BL22','BL23','BL24'};
allsubj = {'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10','BL11','BL12','BL13','BL14','BL15','BL16','BL17','BL18','BL19','BL20','BL21','BL22','BL23','BL4'};

%Low Light
subject_folder_L = {'BL1\BL1_1','BL2\BL2_L','BL3\BL3_L','BL4\BL4_L','BL5\BL5_L','BL6\BL6_L','BL7\BL7_L','BL8\BL8_L','BL9\BL9_L','BL10\BL10_L','BL11\BL11_L','BL12\BL12_L','BL13\BL13_L','BL14\BL14_L','BL15\BL15_L','BL16\BL16_L','BL17\BL17_L','BL18\BL18_L','BL19\BL19_L','BL20\BL20_L','BL21\BL21_L','BL22\BL22_L','BL23\BL23_L','BL24\BL24_L'};
allsubj_L = {'BL1_1','BL2_L','BL3_L','BL4_L','BL5_L','BL6_L','BL7_L','BL8_L','BL9_L','BL10_L','BL11_L','BL12_L','BL13_L','BL14_L','BL15_L','BL16_L','BL17_L','BL18_L','BL19_L','BL20_M','BL21_L','BL22_L','BL23_L','BL4_L'};
%FFT_L = {'FFT_L1','FFT_L2', 'FFT_L3', 'FFT_L4', 'FFT_L5', 'FFT_L6', 'FFT_L7', 'FFT_L8', 'FFT_L9', 'FFT_L10', 'FFT_L11', 'FFT_L12', 'FFT_L13', 'FFT_L14', 'FFT_L15','FFT_L16','FFT_L17','FFT_L18','FFT_L19','FFT_L20','FFT_L21','FFT_L22','FFT_L23','FFT_L24'}; 
%Medium Light
subject_folder_M = {'BL1\BL1_M','BL2\BL2_M','BL3\BL3_M','BL4\BL4_M','BL5\BL5_M','BL6\BL6_M','BL7\BL7_M','BL8\BL8_M','BL9\BL9_M','BL10\BL10_M','BL11\BL11_M','BL12\BL12_M','BL13\BL13_M','BL14\BL14_M','BL15\BL15_M','BL16\BL16_M','BL17\BL17_M','BL18\BL18_M','BL19\BL19_M','BL20\BL20_M','BL21\BL21_M','BL22\BL22_M','BL23\BL23_M','BL24\BL24_M'};
allsubj_M = {'BL1_M','BL2_M','BL3_M','BL4_M','BL5_M','BL6_M','BL7_M','BL8_M','BL9_M','BL10_Meeg','BL11_M','BL13_M','BL13_M','BL14_M','BL15_M','BL16_M','BL17_M','BL18_M','BL19_M','BL20_Mm','BL21_M','BL22_M','BL23_M','BL24_M'};
%High Light
subject_folder_H = {'BL1\BL1_H','BL2\BL2_H','BL3\BL3_H','BL4\BL4_H','BL5\BL5_H','BL6\BL6_H','BL7\BL7_H','BL8\BL8_H','BL9\BL9_H','BL10\BL10_H','BL11\BL11_H','BL12\BL12_H','BL13\BL13_H','BL14\BL14_H','BL15\BL15_H','BL16\BL16_H','BL17\BL17_H','BL18\BL18_H','BL19\BL19_H','BL20\BL20_H','BL21\BL21_H','BL22\BL22_H','BL23\BL23_H','BL24\BL24_H'};
allsubj_H = { 'BL1_H','BL2_H','BL3_H','BL4_H','BL5_H','BL6_H','BL7_H','BL8_H','BL9_H','BL10_H','BL11_H','BL12_H','BL13_H','BL14_H','BL15_H','BL16_H','BL17_H','BL18_H','BL19_H','BL20_H','BL21_H','BL22_H','BL23_H','BL24_H'};


duds = []; %Participants who you want to leave out of the analysis
single_participants = []; %If you want to run the analysis on single Participants 
file_start = 1; %The participant number you want to start the analysis at 

if ~isempty(duds) && isempty(single_participants)
    subject_folder([duds]) = [];
    allsubj([duds]) = [];
    subject_folder_L([duds]) = [];
    allsubj_L([duds]) = [];
    subject_folder_M([duds]) = [];
    allsubj_M([duds]) = [];
    subject_folder_H([duds]) = [];
    allsubj_H([duds]) = [];
    
end

if ~isempty(single_participants)
    subject_folder = subject_folder(single_participants);
    allsubj = allsubj(single_participants);
    subject_folder_L = subject_folder_L(single_participants);
    allsubj_L = allsubj_L(single_participants);
    subject_folder_M = subject_folder_M(single_participants);
    allsubj_M = allsubj_M(single_participants);
    subject_folder_H = subject_folder_H(single_participants);
    allsubj_H = allsubj_H(single_participants);
end


targcodes = [101:124]; % Continuous Dots

old_fs = 500; % old sample rate
fs = 500; %new sample rate

nchan = 65;

LPFcutoff=35;       % Low Pass Filter cutoff

LPF = 1;    % 1 = low-pass filter the data, 0=don't.

% how much of the spectrum to use?
speclims = [1 LPFcutoff];  % Limits in Hz

chanlocs = readlocs ('actiCAP65_ThetaPhi.elp','filetype','besa'); %DN for actiCAP with reference channel 'FCz' included - hence 65 chans
chanlocs = chanlocs(1:nchan)';


for s=1:length(allsubj)
    disp(['Subject: ',num2str(s)])
    disp(['Subject: ',allsubj{s}])
      
    for LC=1:3 %LC= light condition (low medium high)
        clear chanVar
        
        if LC==1
           load([path subject_folder_L{s} '\' allsubj_L{s} '_8_to_13Hz_neg1000_to_1000_ARchans1to65_35HzLPF_point0HzHPF_ET.mat']); 
          
       
        elseif LC==2
          load([path subject_folder_M{s} '\' allsubj_M{s} '_8_to_13Hz_neg1000_to_1000_ARchans1to65_35HzLPF_point0HzHPF_ET.mat']);
          
        elseif LC==3
          load([path subject_folder_H{s} '\' allsubj_H{s} '_8_to_13Hz_neg1000_to_1000_ARchans1to65_35HzLPF_point0HzHPF_ET.mat']);
          
       
        end
  
    
    %% Reshaping and FFT analysis
%     So in the EEG data from each electrode (1:65), on each participant (1:24), 
%     we need to extract the EEG power spectrum (1Hz - 35Hz) seperately for 
%     each of the 3 light conditions (LC). We are only interested on the
%     relevant data where the participant did the right procedure.
 new_erp = erp(:,find(t_crop==-500):find(t_crop==0),rejected_trial_n~=1 & artifact_neg1000_to_0ms_n~=1); % Includes only relevant data
 erp_2 = reshape(new_erp,nchan,size(new_erp,2)*size(new_erp,3)); 
 FFT_amplitude_spectrum = abs(fft(erp_2'))'/(size(erp_2,2)'/2); % FFT amplitude spectrum
 Frequency_scale = [0:size(FFT_amplitude_spectrum,2)-1]*fs/size(FFT_amplitude_spectrum,2); % Frequency scale

    %% Butterworth Filter for smoothier FFT      
period = 0.2; % Frequency period empiracally determined according to the desired number of points for statistical analysis

Frequency_Limit = 35; %Desired frequency limit for the plot

[b,a] = butter(2,0.004); % Parameters of the Butterworth filter 
% % %  The cutoff frequency Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding
% % %  to half the sample rate defined as 500 Hz. This way the cutoff frequency
% % %  was set by Wn = 2*w/500. W is the non normalized cutoff frequency and it
% % %  was empirically defined as 1 Hz.
  for i = 1:nchan    
      
    Filtered_FFT = filtfilt(b,a,FFT_amplitude_spectrum(i,:)); % This function applys the butterworth parameters to the desired signal which in this case is done for each channel
    limit_vectorPos = Frequency_Limit/(Frequency_scale(2)-Frequency_scale(1)); % Position of the last data in the vector according to the chosen frequency limit
    step = round((limit_vectorPos+1)/Frequency_Limit*period); %  Sample Period across the vector
    newFreq_scale = zeros(1,Frequency_Limit/period); % Allocation of variables
    rawFFT_amplitude_spectrum = zeros(1,Frequency_Limit/period);
    newFFT_amplitude_spectrum2 = zeros(1,Frequency_Limit/period);
    newfilter = zeros(1,Frequency_Limit/period);

    for j = 1:Frequency_Limit/period         % Performs a selection of the data for its shorter version
       
        newFreq_scale(j) = Frequency_scale(j*step); % Adjusts the frequency scale to match the reduced version (X axis)
        newFFT_amplitude_spectrum2(j)=Filtered_FFT(j*step); % Takes especific data according to the defined step
    end
  
  FULL_filtered_fft(s,i,:,LC) = newFFT_amplitude_spectrum2; % stores the filtered fft for each channel and light condition of each participant (full data set)
  
  end
  
    end

end


%% Procedure to take the standard deviation error and difference of intensity between high and low ligh condition

light_tag = {'Low Light','Med Light','High Light'}; 
chanlocs = readlocs ('actiCAP65_ThetaPhi.elp','filetype','besa');
FFT_reduced_participants = squeeze(mean(FULL_filtered_fft(:,:,10:5:175,:),1)); %average accross participants
FFT_reduced_channels = squeeze(mean(FULL_filtered_fft(:,[63,64],10:5:175,:),2)); %avarage accross Channels. For the purpose of this experiment the channels 63 and 64 were selected.
for participant = 1:size(FFT_reduced_channels,1)
    for Frequency = 1:size(FFT_reduced_channels,2)
        for LC = 1:size(FFT_reduced_channels,3)
            ste(Frequency,LC) = std(FFT_reduced_channels(:,Frequency,LC))/sqrt(size(FFT_reduced_channels,1)); %These loops calculate the standard deviation error accross the selected channels.
        end
    end
end
FFT_LC_DIFF = FFT_reduced_participants(:,:,3) - FFT_reduced_participants(:,:,1); %Difference in magnitude between High light and Low light condition accross the selected channels.
set(gca,'NextPlot','replaceChildren');
for i = 1:size(FFT_LC_DIFF,2)
    figure;
    topoplot(FFT_LC_DIFF(:,i),chanlocs,'electrodes','numbers','plotchans',1:65); %Plots a map of scalp data for each frequency
    title(['Power: High minus Low Light: ' num2str(i) ' Hz']);
end






%% Alternantive method to smooth FFT using iFFT

% Y = fft(FFT_amplitude_spectrum(25,:),size(FFT_amplitude_spectrum,2));
% r = 150; % range of frequencies we want to preserve
% rectangle = zeros(size(Y));
% rectangle(1:r+1) = 1;               % preserve low +ve frequencies
% y_half = ifft(Y.*rectangle,size(FFT_amplitude_spectrum,2));   % +ve low-pass filtered signal
% rectangle(end-r+1:end) = 1;         % preserve low -ve frequencies
% y_rect = ifft(Y.*rectangle,size(FFT_amplitude_spectrum,2));   % full low-pass filtered signal
% 
% plot(Frequency_scale,y_half,'b','LineWidth',2); 
% hold on
% plot(Frequency_scale,y_rect,'r','LineWidth',2);
% hold on
% plot(Frequency_scale,FFT_amplitude_spectrum(25,:),'.m','MarkerSize',0.02)
% hold on
% plot(Frequency_scale, filtfilt(b,a,FFT_amplitude_spectrum(25,:)),'g')
% legend('iFFT filter 1','iFFT filter 2','Raw data','Butterworth filter')



%% Stft for epoched data (i.e. the erp matrix)


% %% STFT for epoch around target:
% nchan=65;
% fs=500;
% clear stftC
% TC = [-1500:100:700];%100ms sliding window
% fftlen = 300;
% F = [0:20]*fs/fftlen; %Frequencies
% for tt=1:length(TC)
%     temp = abs(fft(erp(:,find(t_crop>=TC(tt),1)-fftlen/2+[1:fftlen],:),[],2))./(fftlen/2);
%     stftC(:,:,tt,:) = reshape(temp(:,1:length(F),:),[nchan length(F) 1 size(erp,3)]);
% end %stftC(Electrode,Frequency,Time,Trial)
% 
% %Isolate time-range and collapse accross it
% Trange = find(TC>0 & TC<800);
% spec =squeeze(mean(stftC(:,:,Trange,:),3)); %spec(Electrode,Frequency,Trial)
% %Isolate frequency band within that time range and collapse accross it:
% band = find(F>8 & F<14);
% 
% PreAlpha=squeeze(mean(spec(:,band,:),2)); %PreAlpha(electrode trial)
% 
% 
% Alpha_simon=squeeze(mean(stftC(:,band,:,:),2)); %Alpha_simon(Electrode,Time,Trial)
% % 





 