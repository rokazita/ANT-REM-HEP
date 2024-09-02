%% this part needs to be done on every subject separately 
% Selecting and segmenting 4 sec epochs of EEG, EOG, Thalamic, ECG data for later preprocessing
ft_defaults
clc;
clear all;
close all;
% outfolder_p = 'd:\ANT_HEP\phasic_tonic_HEP_segments\phasic';
% outfolder_t = 'd:\ANT_HEP\phasic_tonic_HEP_segments\tonic';
% outfolder_w = 'd:\ANT_HEP\phasic_tonic_HEP_segments\wake';
outfolder_n = 'd:\ANT_HEP\phasic_tonic_HEP_segments\nrem';

% This is the output folder for the figure and xlsx files
patients = {'BA', 'EÉ', 'FSI','HaJu', 'KB', 'KEA', 'MaFe', 'PiRi','PJ', 'TI', 'ToZa','TöTa'};
eeg_kea={'EEG Fp1-G2','EEG F7-G2','EEG T3-G2','EEG T5-G2','EEG Fp2-G2','EEG F8-G2','EEG T4-G2','EEG T6-G2',...
    'EEG F3-G2','EEG C3-G2','EEG P3-G2','EEG O1-G2','EEG F4-G2','EEG C4-G2','EEG P4-G2','EEG O2-G2',...    
    'EEG Fpz-G2','EEG Oz-G2','EEG Ft10-G2','EEG Tp10-G2','EEG Po10-G2',...
    'EEG Ft9-G2','EEG Tp9-G2','EEG Po9-G2'};
eeg_kb={'EEG Fp1-G2', 'EEG F7-G2', 'EEG T3-G2', 'EEG T5-G2', 'EEG Fp2-G2', 'EEG F8-G2',...
    'EEG T4-G2', 'EEG T6-G2', 'EEG F3-G2', 'EEG C3-G2', 'EEG P3-G2', 'EEG O1-G2', 'EEG F4-G2',...
    'EEG C4-G2', 'EEG P4-G2', 'EEG O2-G2', 'EEG Fpz-G2', 'EEG Oz-G2',...
    'EEG Ft10-G2', 'EEG Tp10-G2', 'EEG Po10-G2', 'EEG Ft9-G2', 'EEG Tp9-G2', 'EEG Po9-G2'};
eeg_fsi={'EEG Fp1-G2', 'EEG Ft9-G2', 'EEG Fp2-G2', 'EEG F7-G2', 'EEG F3-G2', 'EEG Ft10-G2',...
    'EEG F4-G2', 'EEG F8-G2', 'EEG T3-G2', 'EEG C3-G2', 'EEG Tp10-G2', 'EEG C4-G2', 'EEG T4-G2',...
    'EEG T5-G2', 'EEG P3-G2', 'EEG Tp9-G2', 'EEG P4-G2', 'EEG T6-G2', 'EEG O1-G2', 'EEG Oz-G2',...
    'EEG O2-G2', 'EEG Zyg2-G2', 'EEG Zyg1-G2'}; 

eeg_piri={ 'EEG Fp2-G2', 'EEG F7-G2', 'EEG F8-G2', 'EEG T3-G2',...
     'EEG T4-G2', 'EEG T5-G2', 'EEG Pz-G2', 'EEG T6-G2', 'EEG O1-G2', 'EEG O2-G2',...
     'EEG Fpz-G2', 'EEG Oz-G2', 'EEG Ft10-G2', 'EEG Ft9-G2', 'EEG Tp10-G2', 'EEG Tp9-G2',...
     'EEG Po10-G2', 'EEG Po9-G2'};

 % read and preprocess edf files 
 % select patient (i)
i=12;

% cd('d:\ANT_HEP\Phasic_Tonic_segments')
% [phasic_codes]=xlsread(strcat(char(patients(i)),'.xlsx'),1);
% [tonic_codes]=xlsread(strcat(char(patients(i)),'.xlsx'),2);
% [wake_codes]=xlsread(strcat(char(patients(i)),'.xlsx'),3);

cd('C:\Users\zital\OneDrive\Asztali gép\ELTE KUTATÁS\data\nrem_excels')
[nrem_codes]=xlsread(strcat(char(patients(i)),'.xlsx'),1);


% change codes into sample number multiplied by sr and artefact window length (2048 * 4)
% phasic_samples=phasic_codes.*8192; %except for PiRi and FSI where it is 4096!!!!!!!!!! otherwise: 8192
% tonic_samples=tonic_codes.*8192;
% wake_samples=wake_codes.*8192;
nrem_samples=nrem_codes.*8192;

cd('C:\Users\zital\OneDrive\Asztali gép\MENTÉS\Asztal Mentés\ELTE KUTATÁS\data')
display(strcat('We are analizing patient:',char(patients(i))));

cfg=[];
cfg.dataset=strcat(char(patients(i)),'.edf');
cfg.dftfilter  = 'yes';
cfg.dftfreq    = [50 100];
cfg.bpfilter = 'yes';
cfg.bpfreq  = [0.7 70]; %0.3
cfg.continuous = 'yes';

%CHOOSE THE MONTAGE for the EEG channels
cfg.channel=eeg_kea;  % change if necessary
cfg.reref='yes';
cfg.refchannel={'EEG Tp10-G2','EEG Tp9-G2'};
% 
% % read the specific segments defined by the trials
% cfg.trl=phasic_samples;
% phasic_eeg=ft_preprocessing(cfg); 
% % select tonic samples
% cfg.trl=tonic_samples;
% tonic_eeg=ft_preprocessing(cfg); 
% % select wake samples
% cfg.trl = wake_samples;
% wake_eeg=ft_preprocessing(cfg); 
cfg.trl = nrem_samples;
nrem_eeg=ft_preprocessing(cfg); 


% select the ANT thalamic channel CHECK CORRECT CHANNELS!!!
cfg.channel={'EEG JTh10-JTh11'};   

cfg.reref='no';  % if no rereferencing is necessary
cfg.reref='yes';
cfg.refchannel={'EEG Th11-G2'};
% 
% cfg.trl=phasic_samples;
% phasic_ant=ft_preprocessing(cfg); 
% 
% cfg.trl=tonic_samples;
% tonic_ant=ft_preprocessing(cfg); 
% 
% cfg.trl=wake_samples;
% wake_ant=ft_preprocessing(cfg);

cfg.trl=nrem_samples;
nrem_ant=ft_preprocessing(cfg);


%remove double reference if necessary
cfg.channel={'EEG Th10-G2'};
% phasic_ant=ft_selectdata(cfg,phasic_ant);
% tonic_ant=ft_selectdata(cfg,tonic_ant);
% wake_ant=ft_selectdata(cfg,wake_ant);
nrem_ant=ft_selectdata(cfg,nrem_ant);


% select EOG
cfg.channel={'EEG Zyg2-Zyg1'};   

cfg.reref='no';

cfg.reref='yes';
cfg.refchannel={'EEG ZYG1-G2'};

% cfg.trl=phasic_samples;
% phasic_eog=ft_preprocessing(cfg);
% 
% cfg.trl=tonic_samples;
% tonic_eog=ft_preprocessing(cfg);
% 
% cfg.trl=wake_samples;
% wake_eog=ft_preprocessing(cfg);

cfg.trl=nrem_samples;
nrem_eog=ft_preprocessing(cfg);

% remove mastoids from the eog data
cfg.channel={'EEG ZYG2-G2'};
% phasic_eog=ft_selectdata(cfg,phasic_eog);
% tonic_eog=ft_selectdata(cfg,tonic_eog);
% wake_eog=ft_selectdata(cfg,wake_eog);
nrem_eog=ft_selectdata(cfg,nrem_eog);


cfg.channel={'all' '-EEG Zyg2-Zyg1'};
% phasic_eeg=ft_selectdata(cfg,phasic_eeg);
% tonic_eeg=ft_selectdata(cfg,tonic_eeg);
% wake_eeg=ft_selectdata(cfg,wake_eeg);
nrem_eeg=ft_selectdata(cfg,nrem_eeg);

% select ECG
cfg.channel={'EEG ECG-G2', 'EEG Tp10-G2', 'EEG Tp9-G2'};  
cfg.reref='yes';
cfg.refchannel={'EEG Tp10-G2', 'EEG Tp9-G2'};

% cfg.trl=phasic_samples;
% phasic_ecg=ft_preprocessing(cfg);
% 
% cfg.trl=tonic_samples;
% tonic_ecg=ft_preprocessing(cfg);
% 
% cfg.trl=wake_samples;
% wake_ecg=ft_preprocessing(cfg);

cfg.trl=nrem_samples;
nrem_ecg=ft_preprocessing(cfg);

% remove mastoids from the ecg data
cfg.channel={'EEG ECG-G2'};
% phasic_ecg=ft_selectdata(cfg,phasic_ecg);
% tonic_ecg=ft_selectdata(cfg,tonic_ecg);
% wake_ecg=ft_selectdata(cfg,wake_ecg);
nrem_ecg=ft_selectdata(cfg,nrem_ecg);


% append data
cfg=[];
cfg.keepsampleinfo  = 'yes';

% eeg_ant_ecg_phasic= ft_appenddata(cfg,phasic_eeg,phasic_ant,phasic_ecg); %,phasic_eog);
% 
% eeg_ant_ecg_tonic  = ft_appenddata(cfg,tonic_eeg,tonic_ant,tonic_ecg); %,tonic_eog);
% 
% eeg_ant_ecg_wake= ft_appenddata(cfg,wake_eeg,wake_ant,wake_ecg); %,wake_eog);

eeg_ant_ecg_nrem= ft_appenddata(cfg,nrem_eeg,nrem_ant,nrem_ecg); %,wake_eog);


% segment data to non-overlapping and overlapping data
cfg=[];
cfg.length=4;  % in sec
cfg.overlap=0;
% phasic_segmented = ft_redefinetrial(cfg, eeg_ant_ecg_phasic);
% tonic_segmented = ft_redefinetrial(cfg, eeg_ant_ecg_tonic);
% wake_segmented =ft_redefinetrial(cfg,eeg_ant_ecg_wake);
nrem_segmented =ft_redefinetrial(cfg,eeg_ant_ecg_nrem);



%% DOWNSAMPLING
cfg=[];
cfg.resamplefs      = 512;
% tonic_resampled = ft_resampledata(cfg, tonic_segmented);
% phasic_resampled = ft_resampledata(cfg, phasic_segmented);
% wake_resampled = ft_resampledata(cfg, wake_segmented);
nrem_resampled = ft_resampledata(cfg, nrem_segmented);


outfolder=('d:\ANT_HEP\resampled_512_70hz'); % PATH!
% save(strcat(char(outfolder),'/',char(patients(i)),'_clean_tonic_512.mat'),'tonic_resampled'); 
% save(strcat(char(outfolder),'/',char(patients(i)),'_clean_phasic_512.mat'),'phasic_resampled'); 
% save(strcat(char(outfolder),'/',char(patients(i)),'_clean_wake_512.mat'),'wake_resampled'); 
save(strcat(char(outfolder),'/',char(patients(i)),'_clean_nrem_512.mat'),'nrem_resampled'); 

      
%% SEGMENTING
cfg=[];
cfg.method          = 'trial';
cfg.metric          = 'var';
cfg.xlim            = 100; 
cfg.ylim(1) = -100;
cfg.ylim(2) = 100;
% phasic_clean=ft_rejectvisual(cfg,phasic_resampled); 
% tonic_clean=ft_rejectvisual(cfg,tonic_resampled);
% wake_clean =ft_rejectvisual(cfg,wake_resampled);
nrem_clean =ft_rejectvisual(cfg,nrem_resampled);
% outfolder_p_512 = 'f:\ANT_HEP\HEP_segments_512_70hz\phasic';
% outfolder_t_512 = 'f:\ANT_HEP\HEP_segments_512_70hz\tonic';
% outfolder_w_512 = 'f:\ANT_HEP\HEP_segments_512_70hz\wake';
outfolder_n_512 = 'd:\ANT_HEP\HEP_segments_512_70hz\nrem';

% save(strcat(char(outfolder_p_512),'\',char(patients(i)),'_clean_phasic.mat'),'phasic_clean');
% 
% save(strcat(char(outfolder_t_512),'\',char(patients(i)),'_clean_tonic.mat'),'tonic_clean');
% 
% save(strcat(char(outfolder_w_512),'\',char(patients(i)),'_clean_wake.mat'),'wake_clean');
save(strcat(char(outfolder_n_512),'\',char(patients(i)),'_clean_nrem.mat'),'nrem_clean');


%% preprocessing: detection of R-peaks in ECG signals


% select data file 
%select the files that were saved in the previous part
[file, path] = uigetfile ('*.mat');

% get file name
point = find(file == '.');
filename = file(1:point-1);

% load the data
load ([path,file]);

% detect phasic or tonic and remoredundant channels 
if file(end-6) == 's' 
    cond = 'phasic_clean';

elseif file(end-6) == 'n' 
    cond = 'tonic_clean';
    
elseif file(end-6) == 'a' 
    cond = 'wake_clean';

elseif file(end-6) == 'r' 
    cond = 'nrem_clean';

else 
    cond="selection error";
end

% create trials variable and get number of trials
trials = eval([cond '.trial']);
ntrials = length(trials);

channels= eval([cond '.label']);
nchannels=length(channels);

% find the ECG channel
str = eval([cond '.label']);
ecgchn = find(strcmpi(str,'EEG ECG-G2')); %PIRI'S IS DIFFERENT


% plot 1st trial ecg and decide if necessary to invert
ecg = eval(['trials{1,' int2str(30) '}(' int2str(ecgchn) ',:)']);
h = figure; 
plot(ecg);
% input any number to invert or press enter to skip
invert = input ('invert signal? ');
close (h);

% counter for valid epochs
epoch_count = 1;

%loop ntrials
for i = 1:ntrials
    % data matrix only with ecg and eeg
    data = trials{1,i};
    
    % detect Rpeaks trial-by-trial
    % get ecg signal
    ecg = data(ecgchn,:);
    
    % invert if necessary
    if ~isempty(invert)
        ecg = ecg*(-1);
    end
    
    % plot the ecg trial
    h = figure ('Name',['ECG_trial_' int2str(i)]);
    plot (ecg)
    
    % manually detect R peaks
    x = ginput;
    close (h)
    
    % counter for detected epochs
    count = 1;
    
    % find maximum 100 samples around each selection point and decide if
    % they are still within the 4 sec window
    if ~isempty(x)
        for j = 1:length(x)
            t = floor(x(j));
            if t-100>0 && t+100<2048           % segments are 4 sec long max datapoints: 4 x 2048: 8192 
                [m,ind_max] = max(ecg(t-100:t+100));
                r(count) = ind_max + (t - 101);
                count = count+1;
            end
        end
        
      
        % number of epochs per trial
        nepochs = length(r);
        
        % loop nepochs
        for j = 1:nepochs
            % check if epoch is within matrix boundaries
            if r(j)-101>0 && r(j)+410<2048                                  
%                  if r(j)-410>0 && r(j)+1640<8192                                   % -200 and +800 ms (sr:2048 Hz)!!!

                epochdata(:,:,epoch_count) = data(:,r(j)-101:r(j)+410);
                epoch_count = epoch_count+1;
                [mx,imax] = max(data(1,:));
            end
        end
    end
    if exist('r')==1
        save_r{i}=r;
    else
        save_r{i}=0;
    end
    clear r;
end

cd('d:\ANT_HEP\HEP_epochdata_512_70hz')  % PATH!
save ([filename '_hep_epochs'], 'epochdata');

save ([filename '_hep_epochs_rpeaks'], 'save_r');
