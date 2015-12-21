clear all
close all
clc

%change

path='D:\Monash\Semester 2\Internship\BlueLightData\';
% path='C:\Users\Dan\Dropbox\Monash\PhD\Blue Light\Data\';

subject_folder = {'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10','BL11','BL12','BL13','BL14','BL15','BL16','BL17','BL18','BL19','BL20','BL21','BL22','BL23','BL24'};
allsubj=         {'BL1','BL2','BL3','BL4','BL5','BL6','BL7','BL8','BL9','BL10','BL11','BL12','BL13','BL14','BL15','BL16','BL17','BL18','BL19','BL20','BL21','BL22','BL23','BL4'};

%Low Light
subject_folder_L = {'BL1\BL1_1','BL2\BL2_L','BL3\BL3_L','BL4\BL4_L','BL5\BL5_L','BL6\BL6_L','BL7\BL7_L','BL8\BL8_L','BL9\BL9_L','BL10\BL10_L','BL11\BL11_L','BL12\BL12_L','BL13\BL13_L','BL14\BL14_L','BL15\BL15_L','BL16\BL16_L','BL17\BL17_L','BL18\BL18_L','BL19\BL19_L','BL20\BL20_L','BL21\BL21_L','BL22\BL22_L','BL23\BL23_L','BL24\BL24_L'};
allsubj_L = {'BL1_1','BL2_L','BL3_L','BL4_L','BL5_L','BL6_L','BL7_L','BL8_L','BL9_L','BL10_L','BL11_L','BL12_L','BL13_L','BL14_L','BL15_L','BL16_L','BL17_L','BL18_L','BL19_L','BL20_M','BL21_L','BL22_L','BL23_L','BL4_L'};
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
          EEG = pop_loadbv([path subject_folder_L{s} '\'], [allsubj_L{s} '.vhdr']); 
          
           allbadchans = {[5,59]%BL1\BL1_1
                      [59]%BL2\BL2_L
                      []%BL3\BL3_L
                      [41,46,22,16,28]%BL4\BL4_L
                      [56,19,33,62]%BL5\BL5_L
                      [41,17,46,22,37,58,19,34,29]%BL6\BL6_L 
                      [1,8,12]%BL7\BL7_L 
                      [37,41,45,50]%BL8\BL8_L 
                      [7,36,37,41,46,45]%BL9\BL9_L 
                      [17,22,28,37]%BL10\BL10_L 
                      [45,37,6]%BL11\BL11_L 
                      [32,37,45]%BL12\BL12_L
                      [7,40,11,16,45,37]%{'BL13\BL13_L'
                      [12,37,29]%{'BL14\BL14_L'
                      [48,53,17,22,28,41,60]%{'BL15\BL15_L'
                      [8,13,42]%'BL16\BL16_L'
                      [45,12,60,64]%'BL17\BL17_L'
                      [12,17,37,64,55,32,16]%%'BL18\BL18_L'
                      [37,45,7,12,46]%%'BL19\BL19_L'
                      [35,45,16]%%'BL20\BL20_L'
                      [37,45,64,42,41]%%'BL21\BL21_L'
                      [45,61,5]%%'BL22\BL22_L'
                      [45,34,13,15,58,16,56]%%'BL23\BL23_L'
                      [12,20,42,51,56,41,7,46,30,31]};%%'BL24\BL24_L'
        elseif LC==2
          EEG = pop_loadbv([path subject_folder_M{s} '\'], [allsubj_M{s} '.vhdr']);
          
          allbadchans = {[5,16,55]%BL1\BL1_M
                      [17,28,49,64]%BL2\BL2_M
                      [17,36]%BL3\BL3_M
                      [1,56,12,33]%BL4\BL4_M
                      [41,46,56]%BL5\BL5_M
                      [17,22,41,46,23,28,60,64,27]%BL6\BL6_M
                      [37,45]%,BL7\BL7_M
                      [37,45,50]%BL8\BL8_M
                      [41,46,27,59,45,8]%BL9\BL9_M
                      [37,17,12,60,33]%BL10\BL10_M
                      [28,17,41]%BL11\BL11_M 
                      [37,45,32,51,64,10,12]%BL12\BL12_M 
                      [45,34,35,39,11,38,46]%{'BL13\BL13_M'
                      [12,28,33]%{'BL14\BL14_M'
                      [43,28,32]%{'BL15\BL15_M'
                      [1,23,40]%BL16\BL16_M
                      [45,37,64,61,60]%'BL17\BL17_M'
                      [45,64,60,31,29]%'BL18\BL18_M'
                      [56,62,32,16,46]%%'BL19\BL19_M'
                      [17,22,37,49]%%'BL20\BL20_M'
                      [45,61]%%'BL21\BL21_M'
                      [37,45,16]%%'BL22\BL22_M'
                      [37,45,16,28]%%'BL23\BL23_M'
                      [11,12,17,23,60,46,36,51,41,28]};%%'BL24\BL24_M'
        elseif LC==3
          EEG =pop_loadbv([path subject_folder_H{s} '\'], [allsubj_H{s} '.vhdr']);
          
          allbadchans = {[55,12]%BL1\BL1_H 
                      [28,17,60,41]%BL2\BL2_H %AC
                      [1,3,22,55,31,32]%BL3\BL3_H %AC
                      [41,46,37,3,17,22]%BL4\BL4_H %AC
                      [46,41,22,60,43,37,3,17,16,64]%BL5\BL5_H %AC (4 SWITCHED OFF)
                      [17,22,41,46,23,28,60,64,27]%BL6\BL6_H %AC (4 SWITCHED OFF)
                      [45]%BL7\BL7_H %AC
                      [45,11,24,53,46,7]%BL8\BL8_H %AC_DEFINITELY_CHECK_THIS_ONE!!!
                      [2,41,46,45,37,59,3]%BL9\BL9_H %AC
                      [45,46,37,12,7,17,28]%BL10\BL10_H %AC
                      [45,44,22]%BL11\BL11_H
                      [37,32,28,41,42,8,3,12]%BL12\BL12_H
                      [10,40]%'BL13\BL13_H'
                      [1,7,47,36,37,45,12,48,24,29,43]%'BL14\BL14_H'
                      [33,37,40,42,28]%'BL15\BL15_H'
                      [45,22,50]%BL16\BL16_H
                      [1,37,45,32]%'BL17\BL17_H'
                      [45]%'BL18\BL18_H'
                      [37,8,61,16]%%'BL19\BL19_H'
                      [12,16,22,42,28]%%'BL20\BL20_H'
                      [37,45,64]%%'BL21\BL21_H'
                      [37,45,17,22,32]%%'BL22\BL22_H'
                      [45,20,16]%%'BL23\BL23_H'
                      [2,41,33,42,12,36,46]};%%'BL24\BL24_H'
        end
        
   
    
   
    loadbvSK_DN % From loadbvSK that simon kelly wrote. This cleans and detrends 
                % the data around the DCC events  
    
    EEG = letterkilla_old(EEG); %DN: removes the letters that Brain Products appends to the triggers
    
    % First LP Filter
%     if LPF, EEG.data = eegfilt(EEG.data,old_fs,0,LPFcutoff); end
    if LPF, EEG = pop_eegfiltnew(EEG, 0, LPFcutoff); end %new FIR filter
    
     
    % interpolate bad channels
    badchans = allbadchans{s}; %for interpolating bad channels
    if ~isempty(badchans)
        EEG.chanlocs = chanlocs;
        EEG=eeg_interp(EEG,[badchans],'spherical');
    end
        
    EEG.data = double(EEG.data);
       
    % average-reference the whole continuous data (safe to do this now after interpolation):
    EEG.data = EEG.data - repmat(mean(EEG.data([1:nchan],:),1),[nchan,1]);
    

    %%
    numev = length(EEG.event);
    
    % Fish out the event triggers and times
    clear trigs stimes RT motion_on
    for i=1:numev
        trigs(i)=EEG.event(i).type;
        stimes(i)=round(EEG.event(i).latency);
    end
    
    targtrigs = [];
    for i=1:length(trigs)
        if any(targcodes(:)==trigs(i))
            targtrigs = [targtrigs i];
        end
    end
    
    %% Michel and Ralph!
%     So this is where I am up to with this analysis. And this is where I pass
%     it over for you guys to have some fun with.... :-D
%     WHAT WE NEED:
%     So in the EEG data from each electrode (1:65), on each participant (1:24), 
%     we need to extract the EEG power spectrum (1Hz - 35Hz) seperately for 
%     each of the 3 light conditions (LC). We are only interested on the EEG 
%     data recorded after the first target trigger and before the last target 
%     trigger, so from: stimes(targtrigs(1)):stimes(targtrigs(end))
%     You'll see I've just started playing around with fft, but I will now pass it
%     over to you guys to have some fun with:
    
        FFT_amplitude_spectrum = abs(fft(EEG.data(:,stimes(targtrigs(1)):stimes(targtrigs(end)))'))'; % FFT amplitude spectrum
        Frequency_scale = [0:size(FFT_amplitude_spectrum,2)-1]*EEG.srate/size(FFT_amplitude_spectrum,2); % Frequency scale
        %chanVar = mean(FFT_amplitude_spectrum(:,find(Frequency_scale>speclims(1) & Frequency_scale<speclims(2))),2);       % ROW of average variances for each channel 
         test=(FFT_amplitude_spectrum(:,find(Frequency_scale>speclims(1) & Frequency_scale<speclims(2))));
         %plot(test(25,:))%just choose channel 25 (Pz)
         plot(test)
         
         
         
   %%      
period=0.2; % Frequency period
Frequency_Limit=35;
[b,a]=butter(2,0.004); % Parameters of the Butterworth filter 
%  The cutoff frequency Wn must be 0.0 < Wn < 1.0, with 1.0 corresponding
%  to half the sample rate defined as 500 Hz. This way the cutoff frequency
%  was set by Wn = 2*w/500. W is the non normalized cutoff frequency and it
%  was empirically defined as 1 Hz.
Filtered_FFT=filtfilt(b,a,FFT_amplitude_spectrum(25,:)); % This function applys the butterworth parameters to the desired signal
%figure(1)
%plot(Frequency_scale,abs(Filtered_FFT),'b','LineWidth',3)
% hold on
% plot(Frequency_scale,FFT_amplitude_spectrum(25,:),'r.','MarkerSize',0.01)
% legend('Filtered FFT','Raw FFT')
% xlabel('Frequency(Hz)')
% ylabel('MicroVolts')

limit_vectorPos=Frequency_Limit/(Frequency_scale(2)-Frequency_scale(1)); % Position of the last data in the vector according to the chosen frequency limit
step=round((limit_vectorPos+1)/Frequency_Limit*period); % Period across the vector
newFreq_scale=zeros(Frequency_Limit/period); % Allocation of variables
newFFT_amplitude_spectrum=zeros(Frequency_Limit/period);
newFFT_amplitude_spectrum2=zeros(Frequency_Limit/period);
newfiltro=zeros(Frequency_Limit/period);
for i=1:Frequency_Limit/period         
    newFreq_scale(i)=Frequency_scale(i*step); % Assigning data to the new vectors for a limited and compressed FFT version.
    newFFT_amplitude_spectrum(i)=FFT_amplitude_spectrum(25,i*step);
    newFFT_amplitude_spectrum2(i)=Filtered_FFT(i*step);
    newfiltro(i)=Filtered_FFT(i*step);
end


figure(2)
plot(newFreq_scale,newFFT_amplitude_spectrum,'r')
hold on
plot(newFreq_scale,newFFT_amplitude_spectrum2,'g+')
legend('Limited filtered FFT', 'Limited raw FFT')
 %%        

test=0;
         
         
% By the way Michel and Ralph, for now we will pull out at the data from 
%     all 65 channels, but if you want to just look at one channel to see 
%     if the plots look reasonable, then channel number [25] is a good one
%     to look at (25 is channel Pz).....
% Also don't forget that you'll have to download the EEGLAB toolbox for matlab
% Because the script above uses functions from there. You can download EEGLAB
% for free from here: http://sccn.ucsd.edu/eeglab/install.html
               
    end    
end

%save([path subject_folder{s} '\' allsubj{s} '_PowerSpectrum'],...
%    'speclims','Frequency_scale','FFT_amplitude_spectrum')










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



