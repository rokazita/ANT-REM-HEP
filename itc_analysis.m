%% Inter Trial Coherence Analysis on ANT channels
% Perform itc analysis for specific time windows and frequencies

% Clear the command window, workspace, and close all figuresclc;
clc;
clear all;
close all;

% Set up the FieldTrip toolbox
ft_defaults

% List of patient names and sleep-wake phases
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
% sleep_wake_phases = {'phasic', 'tonic','wake'};
sleep_wake_phases = {'nrem'};

% Some initial settings
fs=512;

% Iterate through sleep-wake phases
for i = 1:length(sleep_wake_phases)
    count=0;

    % Iterate through patients
    for j=1:length(patients)
    cd 'd:\ANT_HEP\HEP_epochdata_512_70hz';
    load(strcat(char(patients(j)), '_clean_', char(sleep_wake_phases(i)), '_hep_epochs'));
    
    % Load segment data based on the sleep-wake phase   
    if strcmpi(sleep_wake_phases(i),'tonic') == 1 
           cd 'f:\ANT_HEP\HEP_segments_512_70hz\tonic';
        elseif  strcmpi(sleep_wake_phases(i),'phasic') == 1
            cd 'f:\ANT_HEP\HEP_segments_512_70hz\phasic';
    elseif strcmpi(sleep_wake_phases(i),'wake') == 1
            cd 'f:\ANT_HEP\HEP_segments_512_70hz\wake';
    else
            cd 'd:\ANT_HEP\HEP_segments_512_70hz\nrem';
    end
    segment=load(strcat(char(patients(j)), '_clean_', char(sleep_wake_phases(i))));
    fns=fieldnames(segment);

    % Create a structure for the analysis using epoch data
    file.label=segment.(fns{1}).label;
    file.fsample=segment.(fns{1}).fsample;
    for k = 1:size(epochdata,3)
        file.trial{1,k}=epochdata(:,:,k);
        file.time{1,k}=linspace(-0.2,0.8,fs);
    end
    for d = 0:size(epochdata,3)-1 
        file.sampleinfo(d+1,1)=d*fs+1;
        file.sampleinfo(d+1,2)=(d+1)*fs;
    end
    file.sampleinfo=segment.(fns{1}).sampleinfo;
    file.hdr=segment.(fns{1}).hdr;
    
    % Specify channel configurations based on the patients
   if strcmpi(patients(j),'BA') == 1 
            channel1='EEG JTH9-JTH10';
            channel2='EEG JTH10-JTH11';
            channel3=zeros(5);
            channel4=zeros(5);
            numberofchannels=2;
        elseif strcmpi(patients(j),'HaJu') == 1 
            channel1='EEG JTh10-JTh11';
            channel2=NaN(5);
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=1;
        elseif strcmpi(patients(j),'KB') == 1 
            channel1='EEG BTh0-BTh1';
            channel2='EEG BTh1-BTh2';
            channel3='EEG BTh2-BTh3';
            channel4=NaN(5);
            numberofchannels=3;
        elseif strcmpi(patients(j),'KEA') == 1 
            channel1='EEG BTh0-BTh1';
            channel2='EEG JTh9-JTh10';
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=2;
        elseif strcmpi(patients(j),'MaFe') == 1 
            channel1='EEG BTH1-BTH2';
            channel2='EEG BTH2-BTH3';
            channel3='EEG JTH10-JTH11';
            channel4=NaN(5);
            numberofchannels=3;
        elseif strcmpi(patients(j),'PiRi') == 1 
            channel1='EEG Th10-G2';
            channel2=NaN(5);
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=1;
        elseif strcmpi(patients(j),'PJ') == 1 
            channel1='EEG BTh1-BTh2';
            channel2='EEG BTh2-BTh3';
            channel3='EEG JTh10-JTh11';
            channel4=NaN(5);
            numberofchannels=3;
        elseif strcmpi(patients(j),'TI') == 1 
            channel1='EEG BTh2-BTh3';
            channel2=NaN(5);
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=1;
        elseif strcmpi(patients(j),'TöTa') == 1 
            channel1='EEG JTh10-JTh11';
            channel2=NaN(5);
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=1;
        elseif strcmpi(patients(j),'ToZa') == 1 
            channel1='EEG Jth8-Jth9';
            channel2='EEG Jth9-Jth10';
            channel3='EEG Bth1-Bth2';
            channel4='EEG Bth2-Bth3';
            numberofchannels=4;
        elseif strcmpi(patients(j), 'FSI') == 1
            channel1='EEG JTh3-G2';
            channel2=NaN(5);
            channel3=NaN(5);
            channel4=NaN(5);
            numberofchannels=1;
        else
            disp('Not known');
            numberofchannels=NaN;
        end

    if numberofchannels>=1
      channels{1}=channel1;
    end
    if numberofchannels>=2
      channels{2}=channel2;
    end
    if numberofchannels>=3
      channels{3}=channel3;
    end
    if numberofchannels>=4
      channels{4}=channel4; 
    end   
    
    % Perform ITC analysis
    cfg        = [];
    cfg.method = 'mtmconvol';
    cfg.channel      = channels; 
    cfg.taper='hanning';
    cfg.toi          = -0.125:0.01:0.725; 
    cfg.foi= 7:1:45;
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.15;  
    cfg.output = 'fourier';
    cfg.pad          = 1.2;
    freq       = ft_freqanalysis(cfg, file);

    % Save the result into a new structure based on the sleep-wake phase
    if strcmpi(sleep_wake_phases(i),'phasic') == 1 
            phasic_itc{j}=freq;
    elseif strcmpi(sleep_wake_phases(i),'tonic') == 1 
            tonic_itc{j}=freq;
    elseif strcmpi(sleep_wake_phases(i),'wake') == 1 
           wake_itc{j}=freq;
    else
           nrem_itc{j}=freq;
    end

    % Clear temporary variables
    clear file fns numberofchannels p k d f s segment channel1 channel2 channel3 channel4 channels;
    end
end

% Save the ITC analysis results
cd 'd:\ANT_HEP\HEP_itc_512_70hz' ;
save('nrem_itc_7_45', 'nrem_itc', '-v7.3');

save('tonic_itc_7_45', 'tonic_itc', '-v7.3');
save('phasic_itc_7_45', 'phasic_itc', '-v7.3');
save('wake_itc_7_45', 'wake_itc', '-v7.3');


%%%%% make a new FieldTrip-style data structure containing the ITC
% copy the descriptive fields over from the frequency decomposition

%%% Phasic
load 'f:\ANT_HEP\HEP_itc_512_70hz\phasic_itc_7_45' 
p_freq=phasic_itc;

clear i
for i=1:size(phasic_itc,2)
    p_itc           = [];
    p_itc.label     = p_freq{1,i}.label;
    p_itc.freq      = p_freq{1,i}.freq;
    p_itc.time      = p_freq{1,i}.time;
    p_itc.dimord    = 'chan_freq_time';
    
    F = p_freq{1, i}.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials
    
    % compute inter-trial phase coherence (itpc)
    p_itc.itpc      = F./abs(F);         % divide by amplitude
    p_itc.itpc      = sum(p_itc.itpc,1);   % sum angles
    p_itc.itpc      = abs(p_itc.itpc)/N;   % take the absolute value and normalize
    p_itc.itpc      = squeeze(p_itc.itpc); % remove the first singleton dimension
    itpc_phasic_all{1, i}=p_itc.itpc;
    
    % compute inter-trial linear coherence (itlc)
    p_itc.itlc      = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
    p_itc.itlc      = abs(p_itc.itlc);     % take the absolute value, i.e. ignore phase
    p_itc.itlc      = squeeze(p_itc.itlc); % remove the first singleton dimension
    itlc_phasic_all{1, i}=p_itc.itlc;
end

%%% Tonic
load 'f:\ANT_HEP\HEP_itc_512_70hz\tonic_itc_7_45' 
t_freq=tonic_itc;
clear i
for i=1:size(tonic_itc,2)
    t_itc           = [];
    t_itc.label     = t_freq{1,i}.label;
    t_itc.freq      = t_freq{1,i}.freq;
    t_itc.time      = t_freq{1,i}.time;
    t_itc.dimord    = 'chan_freq_time';
    
    F = t_freq{1, i}.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials
    
    % compute inter-trial phase coherence (itpc)
    t_itc.itpc    = F./abs(F);         % divide by amplitude
    t_itc.itpc      = sum(t_itc.itpc,1);   % sum angles
    t_itc.itpc   = abs(t_itc.itpc)/N;   % take the absolute value and normalize
    t_itc.itpc     = squeeze(t_itc.itpc); % remove the first singleton dimension
    itpc_tonic_all{1, i}=t_itc.itpc;
    % compute inter-trial linear coherence (itlc)
    t_itc.itlc    = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
    t_itc.itlc  = abs(t_itc.itlc);     % take the absolute value, i.e. ignore phase
    t_itc.itlc    = squeeze(t_itc.itlc); % remove the first singleton dimension
    itlc_tonic_all{1, i}=t_itc.itlc;
end

%%% Wake
load 'f:\ANT_HEP\HEP_itc_512_70hz\wake_itc_7_45' 
w_freq=wake_itc;
clear i
for i=1:size(wake_itc,2)
    w_itc           = [];
    w_itc.label     = w_freq{1,i}.label;
    w_itc.freq      = w_freq{1,i}.freq;
    w_itc.time      = w_freq{1,i}.time;
    w_itc.dimord    = 'chan_freq_time';
    
    F = w_freq{1, i}.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials
    
    % compute inter-trial phase coherence (itpc)
    w_itc.itpc    = F./abs(F);         % divide by amplitude
    w_itc.itpc      = sum(w_itc.itpc,1);   % sum angles
    w_itc.itpc   = abs(w_itc.itpc)/N;   % take the absolute value and normalize
    w_itc.itpc     = squeeze(w_itc.itpc); % remove the first singleton dimension
    itpc_wake_all{1, i}=w_itc.itpc;
    % compute inter-trial linear coherence (itlc)
    w_itc.itlc    = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
    w_itc.itlc  = abs(w_itc.itlc);     % take the absolute value, i.e. ignore phase
    w_itc.itlc    = squeeze(w_itc.itlc); % remove the first singleton dimension
    itlc_wake_all{1, i}=w_itc.itlc;
end    
%%%%%

%%% NREM
load 'd:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_7_45' 
n_freq=nrem_itc;
clear i
for i=1:size(nrem_itc,2)
    n_itc           = [];
    n_itc.label     = n_freq{1,i}.label;
    n_itc.freq      = n_freq{1,i}.freq;
    n_itc.time      = n_freq{1,i}.time;
    n_itc.dimord    = 'chan_freq_time';
    
    F = n_freq{1, i}.fourierspctrm;   % copy the Fourier spectrum
    N = size(F,1);           % number of trials
    
    % compute inter-trial phase coherence (itpc)
    n_itc.itpc    = F./abs(F);         % divide by amplitude
    n_itc.itpc      = sum(n_itc.itpc,1);   % sum angles
    n_itc.itpc   = abs(n_itc.itpc)/N;   % take the absolute value and normalize
    n_itc.itpc     = squeeze(n_itc.itpc); % remove the first singleton dimension
    itpc_nrem_all{1, i}=n_itc.itpc;
    % compute inter-trial linear coherence (itlc)
    n_itc.itlc    = sum(F) ./ (sqrt(N*sum(abs(F).^2)));
    n_itc.itlc  = abs(n_itc.itlc);     % take the absolute value, i.e. ignore phase
    n_itc.itlc    = squeeze(n_itc.itlc); % remove the first singleton dimension
    itlc_nrem_all{1, i}=n_itc.itlc;
end    


% Resize for baseline correction (it needs 3 dimensional data)
clear i
for i=1:size(itlc_nrem_all,2)
    if (i==2) || (i==6) || (i==8) || (i==9) || (i==11) %for those patients who have only 1 ant channel
        itlc_nrem_all{1, i}= reshape(itlc_nrem_all{1, i},1,39,86);
        itpc_nrem_all{1, i}=reshape(itpc_nrem_all{1, i},1,39,86);
%         itlc_wake_all{1, i}= reshape(itlc_wake_all{1, i},1,39,86);
%         itpc_wake_all{1, i}=reshape(itpc_wake_all{1, i},1,39,86);
%         itlc_phasic_all{1, i}=reshape(itlc_phasic_all{1, i},1,39,86);
%         itpc_phasic_all{1, i}=reshape(itpc_phasic_all{1, i},1,39,86);
%         itlc_tonic_all{1, i}=reshape(itlc_tonic_all{1, i},1,39,86);
%         itpc_tonic_all{1, i}=reshape(itpc_tonic_all{1, i},1,39,86);
    end
end

% Create a structure that contains all the results [people x phase]
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
sleep_wake_phases = {'nrem'};
% sleep_wake_phases = {'phasic', 'tonic','wake'};
clear z t
for z=1:length(patients)
    for t=1:length(sleep_wake_phases)
        if strcmpi(sleep_wake_phases(t),'phasic') == 1 
                patients_itc{z,t}=phasic_itc{z};
                patients_itc{z,t}.itlc=itlc_phasic_all{1,z};
                patients_itc{z,t}.itpc=itpc_phasic_all{1,z};
        elseif strcmpi(sleep_wake_phases(t),'tonic') == 1 
                patients_itc{z,t}=tonic_itc{z};
                patients_itc{z,t}.itlc=itlc_tonic_all{z};
                patients_itc{z,t}.itpc=itpc_tonic_all{z};
        elseif strcmpi(sleep_wake_phases(t),'wake') == 1 
                patients_itc{z,t}=wake_itc{z};
                patients_itc{z,t}.itpc=itpc_wake_all{z};
                patients_itc{z,t}.itlc=itlc_wake_all{z};
        elseif strcmpi(sleep_wake_phases(t),'nrem') == 1 
                nrem_itc{t,z}=nrem_itc{1,z};
                nrem_itc{t,z}.itpc=itpc_nrem_all{1,z};
                nrem_itc{t,z}.itlc=itlc_nrem_all{1,z};
        end    
    end
end

save('nrem_itc_without_baseline_7_45', 'nrem_itc', '-v7.3');


%% Baseline Correction

clear all

% load('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_without_baseline_7_45.mat');
load('d:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_without_baseline_7_45.mat');

%baselined_patients_itc=patients_itc; %the baselined data will be saved into this variable
baselined_patients_itc=nrem_itc; %the baselined data will be saved into this variable

% Perform baseline correction on the ITPC data
for row=1:size(nrem_itc,1)
    for column=1:size(nrem_itc,2)
        new_patients_itc{row,column}.powspctrm= nrem_itc{row,column}.itpc; %here itpc is chosen
        new_patients_itc{row,column}.dimord= 'chan_freq_time';
        new_patients_itc{row,column}.freq= nrem_itc{row,column}.freq;
        new_patients_itc{row,column}.time= nrem_itc{row,column}.time;
        new_patients_itc{row,column}.label= nrem_itc{row,column}.label;

        cfg=[];
        cfg.baseline     = [-0.2 -0.05];
        cfg.baselinetype = 'relative';
        cfg.parameter = 'powspctrm';
        [freq] = ft_freqbaseline(cfg, new_patients_itc{row,column});
        baselined_patients_itc{row,column}=freq;
 
        clear freq
    end
end

% Save the baseline-corrected itpc data
cd 'd:\ANT_HEP\HEP_itc_512_70hz';
save('nrem_itc_itpc_with_baseline_7_45', 'baselined_patients_itc', '-v7.3');


% Perform the same baseline correction on the ITLC data
clear baselined_patients_itc

baselined_patients_itc=nrem_itc;

for row=1:size(nrem_itc,1)
    for column=1:size(nrem_itc,2)
        new_patients_itc{row,column}.powspctrm= nrem_itc{row,column}.itlc; %here itlc is chosen
        new_patients_itc{row,column}.dimord= 'chan_freq_time';
        new_patients_itc{row,column}.freq= nrem_itc{row,column}.freq;
        new_patients_itc{row,column}.time= nrem_itc{row,column}.time;
        new_patients_itc{row,column}.label= nrem_itc{row,column}.label;

        cfg=[];
        cfg.baseline     = [-0.2 -0.05];
        cfg.baselinetype = 'relative';
        cfg.parameter = 'powspctrm';
        [freq] = ft_freqbaseline(cfg, new_patients_itc{row,column});
        baselined_patients_itc{row,column}=freq;
       
        clear freq
    end
end

% Save the baseline-corrected itlc data
save('nrem_itc_itlc_with_baseline_7_45', 'baselined_patients_itc', '-v7.3');

%% Manual grand average computation

%%%% Run only one of the sections below
%open files WITHOUT baseline correction
clear all
load ('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_without_baseline_7_45');
    %want to work with itpc values?
    name='inter-trial phase coherence';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itpc;
        end
    end
    %or itlc values?
    name='inter-trial linear coherence';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itlc;
        end
    end
%%%% Or
% open files WITH baseline correction 
    %want to work with itpc values?
    clear all
%     load('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itpc_with_baseline_7_45');
    load('d:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_itpc_with_baseline_7_45');
    name='inter-trial phase coherence';
    patients_itc=baselined_patients_itc;
   %or itlc values?
    clear all
%     load('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itlc_with_baseline_7_45');
    load('d:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_itlc_with_baseline_7_45');
    name='inter-trial linear coherence';
    patients_itc=baselined_patients_itc;
%%%%%


%%%% Working with the loaded variables
% Averaging across ANT channels
for row=1:size(patients_itc,1)
    for column=1:size(patients_itc,2)
        patients_itc{row,column}.powspctrm = mean(patients_itc{row,column}.powspctrm,1);
        patients_itc{row,column}.label={'General-ANT'};
    end
end
clear column row

for newrow=1:size(patients_itc,2)
    nrem{newrow}=patients_itc{1,newrow};    
%     phasic{newrow}=patients_itc{newrow,1};
%     tonic{newrow}=patients_itc{newrow,2};
%     wake{newrow}=patients_itc{newrow,3};
end
clear newrow

% Averaging across patients
grandaverage_phasic = zeros(1, size(phasic{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(phasic,2)
    grandaverage_phasic = grandaverage_phasic + phasic{1, i}.powspctrm;
end
grandaverage_phasic = squeeze(grandaverage_phasic / size(phasic, 2)); 
clear i

grandaverage_tonic = zeros(1, size(tonic{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(tonic,2)
    grandaverage_tonic = grandaverage_tonic + tonic{1, i}.powspctrm;
end
grandaverage_tonic = squeeze(grandaverage_tonic / size(tonic, 2)); 
clear i

grandaverage_wake = zeros(1, size(wake{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(wake,2)
    grandaverage_wake = grandaverage_wake + wake{1, i}.powspctrm;
end
grandaverage_wake = squeeze(grandaverage_wake / size(wake, 2)); 
clear i


grandaverage_nrem = zeros(1, size(nrem{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(nrem,2)
    grandaverage_nrem = grandaverage_nrem + nrem{1, i}.powspctrm;
end
grandaverage_nrem = squeeze(grandaverage_nrem / size(nrem, 2)); 
clear i



%%%% Creating a structure for each phase for plotting later on
% this makes it easier to use the plotting functions (imagesc)
phasic_itc_avg=patients_itc{1,1};
phasic_itc_avg.powspctrm = grandaverage_phasic;

tonic_itc_avg=patients_itc{1,2};
tonic_itc_avg.powspctrm = grandaverage_tonic;

wake_itc_avg=patients_itc{1,3};
wake_itc_avg.powspctrm = grandaverage_wake;

nrem_itc_avg=patients_itc{1,1};
nrem_itc_avg.powspctrm = grandaverage_nrem;



% Save the grand averaged data 
% check the name
cd 'd:\ANT_HEP\HEP_itc_512_70hz';
save('grand_average_nrem_itpc_baselined_7_45', 'nrem_itc_avg', '-v7.3');

save('grand_average_phasic_itlc_baselined_7_45', 'phasic_itc_avg', '-v7.3');
save('grand_average_tonic_itlc_baselined_7_45', 'tonic_itc_avg', '-v7.3');
save('grand_average_wake_itlc_baselined_7_45', 'wake_itc_avg', '-v7.3');



%%%%% Plot

% Initialize some settings
itc_s=[nrem_itc_avg];
titles=["NREM "];
font_size=10;
x_label='Time (sec)';
y_label='Frequency (Hz)';

% Set up figure standards
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

% Loop through phases and create subplots
for i=1:length(titles)
    %subplot(3, 1, i);
    fig1_comps.p1=imagesc(itc_s(i).time, itc_s(i).freq, itc_s(i).powspctrm);
    axis xy
    xline(0, '--w', 'LineWidth',2)
    caxis([0.8 2])
    c = colorbar;
    c.Label.String = 'Relative ITC';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(strcat(titles(i), name),'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end

%% Statistical Analysis

%%%% Run only one of the sections below
%open files WITHOUT baseline correction
clear all
load ('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_without_baseline_7_45');
    %want to work with itpc values?
    name='inter-trial phase coherence';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itpc;
        end
    end
    %or itlc values?
    name='inter-trial linear coherence';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itlc;
        end
    end
%%%%or
% open files WITH baseline correction 
    clear all
    load('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itpc_with_baseline_7_45');
    name='inter-trial phase coherence';
    patients_itc=baselined_patients_itc;
    %or!
    clear all
    rem= load('d:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itpc_with_baseline_7_45');
    nrem= load('d:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_itpc_with_baseline_7_45');
    name='inter-trial phase coherence';
    name='inter-trial linear coherence';
   % patients_itc=baselined_patients_itc;
%%%%%

patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

for i=1:length(patients)
    rem.baselined_patients_itc(i,4)=nrem.baselined_patients_itc(1,i);
end
loop_version = {'phasic', 'tonic','wake','nrem'};
for j=1:length(loop_version)
    for i=1:length(patients)
        baselined_patients_itc(i,j)=rem.baselined_patients_itc(i,j);
    end
end



%%%% Continue with loaded variables
% Averaging the ANT channels for itpc or itlc (named powspctrm in this case)
for row=1:size(baselined_patients_itc,1)
    for column=1:size(baselined_patients_itc,2)
        baselined_patients_itc{row,column}.powspctrm=mean(baselined_patients_itc{row,column}.powspctrm,1);
        baselined_patients_itc{row,column}.label={'General-ANT'};
    end
end
clear column row

clear nrem
% Divide the data into new structures based on the sleep-wake phases
for newrow=1:size(baselined_patients_itc,1)
    phasic{newrow}=baselined_patients_itc{newrow,1};
    tonic{newrow}=baselined_patients_itc{newrow,2};
    wake{newrow}=baselined_patients_itc{newrow,3};
    nrem{newrow}=baselined_patients_itc{newrow,4};
end
clear newrow


%%%% Statistical Analysis

patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

% Design matrix
subj = length(patients); 
design = zeros(2,2*subj); 
for i = 1:subj 
    design(1,i) = i; 
end 
for i = 1:subj 
design(1,subj+i) = i; 
end 
design(2,1:subj) = 1; 
design(2,subj+1:2*subj) = 2;


% Parameters for the analysis
cfg = []; 
cfg.design = design;
cfg.parameter='powspctrm'; %=itpc!!! just renamed. could be itlc too.
cfg.method = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability 
cfg.statistic = 'ft_statfun_depsamplesT';  %paired samples t test
cfg.avgoverfreq='no';
cfg.latency=[0.05 0.65];
cfg.avgovertime='no';
cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
cfg.alpha = 0.05; 
cfg.uvar = 1; 
cfg.ivar = 2; 
cfg.numrandomization = 'all'; 

% Perform frequency statistics
[stat_p_n_itpc] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
[stat_t_n_itpc] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
[stat_w_n_itpc] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
% [stat_p_n_itlc] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
% [stat_t_n_itlc] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
% [stat_w_n_itlc] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
% 

% Save the results
% check the naming 
cd 'd:\ANT_HEP\HEP_itc_512_70hz'
save('stat_p_n_itpc_with_baseline_7_45', 'stat_p_n_itpc', '-v7.3');
save('stat_t_n_itpc_with_baseline_7_45', 'stat_t_n_itpc', '-v7.3');
save('stat_w_n_itpc_with_baseline_7_45', 'stat_w_n_itpc', '-v7.3');
% save('stat_p_t_itlc_without_baseline_7_45', 'stat_p_t_itlc', '-v7.3');
% save('stat_p_w_itlc_without_baseline_7_45', 'stat_p_w_itlc', '-v7.3');
% save('stat_t_w_itlc_without_baseline_7_45', 'stat_t_w_itlc', '-v7.3');


%%%% Plot

% Initialize some settings
font_size = 10;
x_label = 'Time (sec)';
y_label = 'Frequency (Hz)';
stat=[stat_p_n_itpc stat_t_n_itpc stat_w_n_itpc];
titles=[" Phasic-NREM", " Tonic- REM", " Wake-NREM"];


%p Plot prob, p value
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

for i=1:3
    subplot(3, 1,i)
    imagesc(squeeze(stat(i).time), squeeze(stat(i).freq), squeeze(stat(i).prob));
    set(gca,'YDir','normal')
    colormap(bone)
    xline(0, '--k', 'LineWidth',2)
    caxis([0.001 0.05])
    c = colorbar;
    c.Label.String = 'p value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(strcat(name,' -' , titles(i)), 'fontsize', 16);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end 
clear fig1_comps

%p Plot stat, T value
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

for i=1:3
    subplot(3,1,i)
    imagesc(squeeze(stat(i).time), squeeze(stat(i).freq), squeeze(stat(i).stat));
    set(gca,'YDir','normal')
    caxis([-1 2])
    xline(0, '--w', 'LineWidth',2)
    c = colorbar;
    c.Label.String = 'T value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(strcat(name, ' -', titles(i)), 'fontsize', 16);
    fig1_comps.plotXLabel = xlabel('Time (sec)') ;
    fig1_comps.plotYLabel =ylabel('Frequency (Hz)');
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end


%% Cluster Analysis

%%%%% Run only one of the sections below
% Open files WITHOUT baseline correction
clear all
load ('d:\ANT_HEP\HEP_itc_512_70hz\patients_itc_without_baseline_7_45');
    %want to work with itpc values?
    name='inter-trial phase coherence';
    baselined='no';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itpc;
        end
    end

    %or itlc values?
    name='inter-trial linear coherence';
    baselined='no';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=patients_itc{i,j}.itlc;
        end
    end
%%%Or
% Open files WITH baseline correction
    %want to work with itpc values?
    clear all
    load('d:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itpc_with_baseline_7_45');
    name='inter-trial phase coherence';
    baselined='yes';
    patients_itc=baselined_patients_itc;
    %or itlc values?
    clear all
    load('f:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itlc_with_baseline_7_45');
    name='inter-trial linear coherence';
    baselined='yes';
    patients_itc=baselined_patients_itc;
%%%%%
clear all
    rem= load('d:\ANT_HEP\HEP_itc_512_70hz\patients_itc_itlc_with_baseline_7_45');
    nrem= load('d:\ANT_HEP\HEP_itc_512_70hz\nrem_itc_itlc_with_baseline_7_45');
    name='inter-trial phase coherence';
    name='inter-trial linear coherence';
   % patients_itc=baselined_patients_itc;
%%%%%

patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

for i=1:length(patients)
    rem.baselined_patients_itc(i,4)=nrem.baselined_patients_itc(1,i);
end
loop_version = {'phasic', 'tonic','wake','nrem'};
for j=1:length(loop_version)
    for i=1:length(patients)
        baselined_patients_itc(i,j)=rem.baselined_patients_itc(i,j);
    end
end





%%%% Continue with loaded variables

% Average the ANT channel for itpc or itlc (named powspctrm in this case)
for row=1:size(baselined_patients_itc,1)
    for column=1:size(baselined_patients_itc,2)
        baselined_patients_itc{row,column}.powspctrm=mean(baselined_patients_itc{row,column}.powspctrm,1);
        baselined_patients_itc{row,column}.label={'General-ANT'};
    end
end
clear column row

clear nrem
% Divide the data into new structures based on the sleep-wake phases
for newrow=1:size(baselined_patients_itc,1)
    phasic{newrow}=baselined_patients_itc{newrow,1};
    tonic{newrow}=baselined_patients_itc{newrow,2};
    wake{newrow}=baselined_patients_itc{newrow,3};
   % nrem{newrow}=baselined_patients_itc{newrow,4};
end
clear newrow

%%%%% Perform cluster analysis

patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

% Create the design matrix
subj = length(patients); 
design = zeros(2,2*subj); 
for i = 1:subj 
    design(1,i) = i; 
end 
for i = 1:subj 
design(1,subj+i) = i; 
end 
design(2,1:subj) = 1; 
design(2,subj+1:2*subj) = 2;

% Set analysis parameters
cfg = []; 
cfg.design = design;
cfg.parameter='powspctrm'; %=itpc!!! it is just renamed. but it could be itlc too.
cfg.method = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability 
cfg.statistic = 'ft_statfun_depsamplesT';  
cfg.avgoverfreq='no';
cfg.avgovertime='no';
cfg.correctm = 'cluster'; 
cfg.clusterstatistic = 'maxsum'; 
cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
cfg.clusteralpha = 0.05; 
cfg.latency=[0.05 0.65];
cfg.uvar = 1; 
cfg.ivar = 2; 
cfg.numrandomization = 'all'; 

% Perform cluster-based statistics
% [stat_p_n_itlc_cluster] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); %tonic vs phasic
% [stat_t_n_itlc_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
% [stat_w_n_itlc_cluster] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
[stat_p_t_itpc_cluster] = ft_freqstatistics(cfg, phasic{1,:}, tonic{1,:}); %tonic vs phasic
[stat_p_w_itpc_cluster] = ft_freqstatistics(cfg,  phasic{1,:}, wake{1,:}); 
[stat_t_w_itpc_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, wake{1,:}); 


% Save the results
% check the naming 
cd 'd:\ANT_HEP\HEP_itc_512_70hz'
save('stat_p_n_itpc_cluster_without_baseline_7_45', 'stat_p_n_itpc_cluster', '-v7.3');
save('stat_t_n_itpc_cluster_without_baseline_7_45', 'stat_t_n_itpc_cluster', '-v7.3');
save('stat_w_n_itpc_cluster_without_baseline_7_45', 'stat_w_n_itpc_cluster', '-v7.3');
% save('stat_p_t_itpc_cluster_without_baseline_7_45', 'stat_p_t_itpc_cluster', '-v7.3');
% save('stat_p_w_itpc_cluster_without_baseline_7_45', 'stat_p_w_itpc_cluster', '-v7.3');
% save('stat_t_w_itpc_cluster_without_baseline_7_45', 'stat_t_w_itpc_cluster', '-v7.3');



%%%%% Select clusters for the plot

% Select clusters for itpc p_t: baselined: positive cluster with p=0.02
% itpc non baselined: positive cluster with p=0 0.0049 0.0146 0.0342 
p_t_pos_cluster_pvals = [stat_p_t_itpc_cluster.posclusters(:).prob];
p_t_pos_clust         = find(p_t_pos_cluster_pvals < 0.05);
stat_p_t_itpc_cluster.posclusterslabelmat=squeeze(stat_p_t_itpc_cluster.posclusterslabelmat);
p_t_clustermatrix= stat_p_t_itpc_cluster.posclusterslabelmat == p_t_pos_clust;

% Select clusters for p_w: baselined: positive cluster with p=0.03
%non baselined: positive cluster with p=0.0439
p_w_pos_cluster_pvals = [stat_p_w_itpc_cluster.posclusters(:).prob];
p_w_pos_clust         = find(p_w_pos_cluster_pvals < 0.05);
stat_p_w_itpc_cluster.posclusterslabelmat=squeeze(stat_p_w_itpc_cluster.posclusterslabelmat);
p_w_clustermatrix = stat_p_w_itpc_cluster.posclusterslabelmat == p_w_pos_clust;

% Select clusters for t_w: baselined: no significant cluster
%non baselined: negative cluster with p=0.0557
t_w_neg_cluster_pvals = [stat_t_w_itpc_cluster.negclusters(:).prob];
t_w_neg_clust         = find(t_w_neg_cluster_pvals < 0.05);
stat_t_w_itpc_cluster.negclusterslabelmat=squeeze(stat_t_w_itpc_cluster.negclusterslabelmat);
t_w_clustermatrix = stat_t_w_itpc_cluster.negclusterslabelmat == t_w_neg_clust;


%%%%% Perform operations on selected clusters
% Subtract phasic and tonic grand averaged data
% Subtract phasic and wake grand averaged data
% Subtract tonic and wake grand averaged data
cd 'f:\ANT_HEP\HEP_itc_512_70hz'
if baselined=="yes"
    load 'grand_average_phasic_itpc_baselined_7_45.mat'
    load 'grand_average_tonic_itpc_baselined_7_45.mat'
    load 'grand_average_wake_itpc_baselined_7_45.mat'
else
    load 'grand_average_phasic_itpc_nonbaselined_7_45.mat'
    load 'grand_average_tonic_itpc_nonbaselined_7_45.mat'
    load 'grand_average_wake_itpc_nonbaselined_7_45.mat'
end
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
% For phasic-tonic
grandaverage_phasic_vs_tonic = ft_math(cfg,phasic_itc_avg,tonic_itc_avg);
stat_p_t_itpc_cluster.raw = (grandaverage_phasic_vs_tonic.powspctrm);
% For phasic-wake
grandaverage_phasic_vs_wake = ft_math(cfg,phasic_itc_avg,wake_itc_avg);
stat_p_w_itpc_cluster.raw = (grandaverage_phasic_vs_wake.powspctrm);
% For tonic-wake
grandaverage_tonic_vs_wake = ft_math(cfg,tonic_itc_avg,wake_itc_avg);
stat_t_w_itpc_cluster.raw = (grandaverage_tonic_vs_wake.powspctrm);



% Load significant differences from previous analysis
cd 'E:\ANT_HEP\HEP_itc_512_70hz'
if baselined=="yes"
    load('stat_p_t_itpc_with_baseline_7_45');
    load('stat_p_w_itpc_with_baseline_7_45');
    load('stat_t_w_itpc_with_baseline_7_45');
else
    load('stat_p_t_itpc_without_baseline_7_45');
    load('stat_p_w_itpc_without_baseline_7_45');
    load('stat_t_w_itpc_without_baseline_7_45');
end
% For phasic-tonic
pt_pvals = squeeze(stat_p_t_itpc.prob);
pt_loc         = pt_pvals < 0.05 ;
simple_boundary_pt = bwboundaries(pt_loc);
% For phasic-wake
pw_pvals = squeeze(stat_p_w_itpc.prob);
pw_loc         = pw_pvals < 0.05 ;
simple_boundary_pw = bwboundaries(pw_loc);
% For tonic-wake
tw_pvals = squeeze(stat_t_w_itpc.prob);
tw_loc         = tw_pvals < 0.05 ;
simple_boundary_tw = bwboundaries(tw_loc);


%%%%% Plot the results

% Get outlines of clusters
boundary_pt = bwboundaries(p_t_clustermatrix);
boundary_pw = bwboundaries(p_w_clustermatrix);
boundary_tw = bwboundaries(t_w_clustermatrix); 
%boundary_tw = {0,0};  %if no clusters were found

% Initialize some settings
font_size=10;
x_label='Time (sec)';
y_label='Frequency (Hz)';

averages=[grandaverage_phasic_vs_tonic grandaverage_phasic_vs_wake grandaverage_tonic_vs_wake];
titles=["Phasic-Tonic ", "Phasic-Wake ", "Tonic-Wake "];
cluster_boundaries={boundary_pt boundary_pw boundary_tw};
simple_boundaries={simple_boundary_pt simple_boundary_pw simple_boundary_tw};


PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

% Loop through phases and create subplots
for i=1:length(titles)
    subplot(3, 1, i);
    fig1_comps.p1=imagesc(averages(i).time, averages(i).freq, averages(i).powspctrm);
    axis xy
    xline(0, '--w', 'LineWidth',2)
    caxis([-1 1.5])
    c = colorbar;
    c.Label.String = 'Relative ITC';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(strcat(titles(i), name),'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    
    hold on; % Allow the plot to overlay on the image
    
    % Plot the most significant clusters boundary as a line
    b = cluster_boundaries{1,i};
    for j=1:length(cluster_boundaries{1,i}) 
       bd=cluster_boundaries{1,i}{j,1};
       plot(averages(i).time(bd(:,2)), averages(i).freq(bd(:,1)), 'w', 'LineWidth', 2);
    end
    hold on

    % Plot the significant differences as a shaded area
    % Loop through boundaries and shade the areas
    for k=1:length(simple_boundaries{1,i})
        d=simple_boundaries{1,i}{k,1};
        fill(averages(i).time(d(:,2)), averages(i).freq(d(:,1)), 'black', 'FaceAlpha', 0.3);
    end
    hold off
end



%% Segmented statistical analysis

%open files WITHOUT baseline correction
clear all
load ('d:\ANT_HEP\HEP_itc_512_70hz\patients_itc_without_baseline_7_45');
    %want to work with itpc values?
    name='inter-trial linear coherence';
    %create powspctrm and average the ANT channels
    for i=1:size(nrem_itc,1)
        for j=1:size(nrem_itc,2)
            nrem_itc{i,j}.powspctrm=mean(nrem_itc{i,j}.itlc,1);
            nrem_itc{i,j}.label={'General-ANT'};
        end
    end
        name='inter-trial phase coherence';
    for i=1:size(patients_itc,1)
        for j=1:size(patients_itc,2)
            patients_itc{i,j}.powspctrm=mean(patients_itc{i,j}.itpc,1);
            patients_itc{i,j}.label={'General-ANT'};
        end
    end

%nrem_itc=nrem_itc';
patients_itc=patients_itc';
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 
sleep_wake_phases = {'phasic', 'tonic', 'wake'};


% Segment into 75 ms segments
baseline_indices=1:16; % baseline: -125 ms - -50ms
start=19; % first segment starts at 50ms

% Creating the segments
for i=1:floor((length(patients_itc{1,1}.time)-start)/length(baseline_indices))
    timepoints=patients_itc{1,1}.time;
    segments{i}=timepoints(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
    segment_indices{i}=(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
end
baseline_window=timepoints(baseline_indices);

baseline_indices=1:8; % baseline: -125 ms - -50ms

% Run the statistical analysis on each segment
for i=1:length(segments)

    % Create the structure for each segment(i)

    for j=1:length(sleep_wake_phases)
        % Cut every patient's data into segments
        % (new variable: phase x patient)
        for k=1:length(patients)
            act_segment{j,k}.label=patients_itc{j,k}.label;
            act_segment{j,k}.dimord=patients_itc{j,k}.dimord;
            act_segment{j,k}.freq=patients_itc{j,k}.freq;
            act_segment{j,k}.cfg=patients_itc{j,k}.cfg;
            act_segment{j,k}.time=patients_itc{j,k}.time(1,segment_indices{i});
            act_segment{j,k}.powspctrm=patients_itc{j,k}.powspctrm(:,:,segment_indices{i});
        
            bl_segment{j,k}.label=patients_itc{j,k}.label;
            bl_segment{j,k}.dimord=patients_itc{j,k}.dimord;
            bl_segment{j,k}.freq=patients_itc{j,k}.freq;
            bl_segment{j,k}.cfg=patients_itc{j,k}.cfg;
            bl_segment{j,k}.time=patients_itc{j,k}.time(1,segment_indices{i}); %using this time interval, because the statistical analysis runs only with the same time points
            c=patients_itc{j,k}.powspctrm(:,:,baseline_indices);
            bl_segment{j,k}.powspctrm=repelem(c,1,1,2);
        end
    end
    clear j k

    % Statistical analysis
    
    % Design matrix
    subj = length(patients); 
    design = zeros(2,2*subj); 
    for m = 1:subj 
        design(1,m) = m; 
    end 
    clear m
    for m = 1:subj 
        design(1,subj+m) = m; 
    end 
    design(2,1:subj) = 1; 
    design(2,subj+1:2*subj) = 2;

    % Set the parameters for the analysis
    cfg = []; 
    cfg.design = design;
    cfg.parameter='powspctrm';
    cfg.method = 'montecarlo'; 
    cfg.statistic = 'ft_statfun_depsamplesT';  % ft_statfun_actvsblt ???
    cfg.correctm = 'cluster'; 
    cfg.clusterstatistic = 'maxsum'; 
    cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
    cfg.clusteralpha = 0.05; 
    cfg.avgoverfreq='no';
    cfg.avgovertime='no';
    cfg.numrandomization = 'all'; 
    cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
    cfg.alpha = 0.05; 
    cfg.uvar = 1; 
    cfg.ivar = 2; 

    % Run the analysis
%     stat_n_itc_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
    stat_p_itc_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
    stat_t_itc_seg{i} = ft_freqstatistics(cfg, act_segment{2,:}, bl_segment{2,:}); 
    stat_w_itc_seg{i} = ft_freqstatistics(cfg, act_segment{3,:}, bl_segment{3,:}); 

end

%Concatenate the segments 
stats=[stat_p_itc_seg; stat_t_itc_seg; stat_w_itc_seg];
%stats=[stat_n_itc_seg];

unified{1,1}=stats{1,1};
unified{2,1}=stats{2,1};
unified{3,1}=stats{3,1};
% 
% for i = 1:size(stats,1) 
%    unified{i,1}.time=cat(2, stats{i,1}.time, stats{i,2}.time, stats{i,3}.time, stats{i,4}.time, stats{i,5}.time, stats{i,6}.time, stats{i,7}.time, stats{i,8}.time);
%    unified{i,1}.prob=cat(3, stats{i,1}.prob, stats{i,2}.prob, stats{i,3}.prob, stats{i,4}.prob, stats{i,5}.prob, stats{i,6}.prob, stats{i,7}.prob, stats{i,8}.prob);
%    unified{i,1}.stat=cat(3, stats{i,1}.stat, stats{i,2}.stat, stats{i,3}.stat, stats{i,4}.stat, stats{i,5}.stat, stats{i,6}.stat, stats{i,7}.stat, stats{i,8}.stat);
%    unified{i,1}.ref=cat(3, stats{i,1}.ref, stats{i,2}.ref, stats{i,3}.ref, stats{i,4}.ref, stats{i,5}.ref, stats{i,6}.ref, stats{i,7}.ref, stats{i,8}.ref);
%    unified{i,1}.cirange=cat(3, stats{i,1}.cirange, stats{i,2}.cirange, stats{i,3}.cirange, stats{i,4}.cirange, stats{i,5}.cirange, stats{i,6}.cirange, stats{i,7}.cirange, stats{i,8}.cirange);
%    unified{i,1}.mask=cat(3, stats{i,1}.mask, stats{i,2}.mask, stats{i,3}.mask, stats{i,4}.mask, stats{i,5}.mask, stats{i,6}.mask, stats{i,7}.mask, stats{i,8}.mask);
% end

for i = 1:size(stats,1) 
   unified{i,1}.time=cat(2, stats{i,1}.time, stats{i,2}.time, stats{i,3}.time, stats{i,4}.time);
   unified{i,1}.prob=cat(3, stats{i,1}.prob, stats{i,2}.prob, stats{i,3}.prob, stats{i,4}.prob);
   unified{i,1}.stat=cat(3, stats{i,1}.stat, stats{i,2}.stat, stats{i,3}.stat, stats{i,4}.stat);
   unified{i,1}.ref=cat(3, stats{i,1}.ref, stats{i,2}.ref, stats{i,3}.ref, stats{i,4}.ref);
   unified{i,1}.cirange=cat(3, stats{i,1}.cirange, stats{i,2}.cirange, stats{i,3}.cirange, stats{i,4}.cirange);
   unified{i,1}.mask=cat(3, stats{i,1}.mask, stats{i,2}.mask, stats{i,3}.mask, stats{i,4}.mask);
end
%this structure can be now plotted

% Plot the results

% Initialize some settings
font_size = 10;
x_label = 'Time (sec)';
y_label = 'Frequency (Hz)';
stat=[unified{1,1} unified{2,1} unified{3,1}];
%titles=["NREM - baseline"];
 titles=["Phasic - baseline vs post-R period", "Tonic -  baseline vs post-R period", "Wake - baseline vs post-R period"];

% Adjusted Code to Create a 3x2 Plot with Columns Swapped
PS = PLOT_STANDARDS();
figure;
for i=1:length(titles)
    % First column: Plot stat, T values
    ax1 = subplot(3, 2, (i-1)*2 + 1); % Creates subplots in a 3x2 grid, first column
    imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).stat));
    hold on
    set(gca,'YDir','normal')
    xline(0, '--w', 'LineWidth',2)
    colormap(ax1, parula)  % Apply 'parula' colormap to this subplot only
    caxis([-1 2])
    c = colorbar;
    c.Label.String = 'T value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    title(titles(i), 'fontsize', font_size, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel(x_label, 'FontName', 'Times New Roman');
    ylabel(y_label, 'FontName', 'Times New Roman');

    % Second column: Plot prob, p values
    ax2 = subplot(3, 2, (i-1)*2 + 2); % Creates subplots in a 3x2 grid, second column
    imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).prob));
    set(gca,'YDir','normal')
    xline(0, '--k', 'LineWidth', 2)
    colormap(ax2, bone)  % Apply 'bone' colormap to this subplot only
    caxis([0.001 0.05])
    c = colorbar;
    c.Label.String = 'p value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    title(titles(i), 'fontsize', font_size, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    xlabel(x_label, 'FontName', 'Times New Roman');
    ylabel(y_label, 'FontName', 'Times New Roman');
end


%% Plot that contains every result from ITC analysis 

clear all

% Choose from the following options:
    baselined="yes";
    %baselined="no";
% and
    measure="itpc";
  %  measure="itlc";

% open the grand averaged data
% then the statistics
% then the cluster based statistics
cd 'd:\ANT_HEP\HEP_itc_512_70hz'
if baselined=="yes"
    %load(strcat('grand_average_nrem_', measure, '_baselined_7_45.mat'));
    load(strcat('grand_average_phasic_', measure, '_baselined_7_45.mat'));
    load(strcat('grand_average_tonic_', measure,'_baselined_7_45.mat'));
    load(strcat('grand_average_wake_', measure,'_baselined_7_45.mat'));
    stat_p_t_itc=load(strcat('stat_p_t_', measure, '_with_baseline_7_45'));
   % stat_p_n_itc=load(strcat('stat_p_n_', measure', '_with_baseline_7_45'));
    %stat_t_n_itc=load(strcat('stat_t_n_', measure', '_with_baseline_7_45'));
    %stat_w_n_itc=load(strcat('stat_w_n_', measure', '_with_baseline_7_45'));
    %stat_p_n_itc_cluster=load(strcat('stat_p_n_', measure', '_cluster_with_baseline_7_45'));
    %stat_t_n_itc_cluster=load(strcat('stat_t_n_', measure', '_cluster_with_baseline_7_45'));
    %stat_w_n_itc_cluster=load(strcat('stat_w_n_', measure', '_cluster_with_baseline_7_45'));
    
    stat_p_w_itc=load(strcat('stat_p_w_', measure', '_with_baseline_7_45'));
    stat_t_w_itc=load(strcat('stat_t_w_', measure', '_with_baseline_7_45'));
    stat_p_t_itc_cluster=load(strcat('stat_p_t_', measure', '_cluster_with_baseline_7_45'));
    stat_p_w_itc_cluster=load(strcat('stat_p_w_', measure', '_cluster_with_baseline_7_45'));
    stat_t_w_itc_cluster=load(strcat('stat_t_w_', measure', '_cluster_with_baseline_7_45'));
else
    load(strcat('grand_average_phasic_', measure, '_nonbaselined_7_45.mat'));
    load(strcat('grand_average_tonic_', measure, '_nonbaselined_7_45.mat'));
    load(strcat('grand_average_wake_', measure, '_nonbaselined_7_45.mat'));
    stat_p_t_itc=load(strcat('stat_p_t_', measure, '_without_baseline_7_45'));
    stat_p_w_itc=load(strcat('stat_p_w_', measure, '_without_baseline_7_45'));
    stat_t_w_itc=load(strcat('stat_t_w_', measure, '_without_baseline_7_45'));
    stat_p_t_itc_cluster=load(strcat('stat_p_t_', measure, '_cluster_without_baseline_7_45'));
    stat_p_w_itc_cluster=load(strcat('stat_p_w_', measure, '_cluster_without_baseline_7_45'));
    stat_t_w_itc_cluster=load(strcat('stat_t_w_', measure, '_cluster_without_baseline_7_45'));
end


% Plot
figure
tcl = tiledlayout(3,3);  
ax = gobjects(1,9);
index=1;
for i = 1:3
     ax(index) = nexttile;
    
    % Plot the grand averaged data in the first column
    % Initialize some settings
    phases=[phasic_itc_avg tonic_itc_avg wake_itc_avg];
    titles=["Phasic", "Tonic", "Wake"];
    z_lim= [0.8 2];
    x_lim=[-0.13 0.6];
    font_size=10;
    x_label='Time (sec)';
    y_label='Frequency (Hz)';
    
    % Set up figure standards
    PS = PLOT_STANDARDS();
    fig1_comps.fig = gcf;
    cfg = [];
    % Pay attention to cfg.baseline, this is optional to use
    % cfg.baseline     = [-0.2 -0.05];
    %  cfg.baselinetype = 'relative';
    cfg.zlim=z_lim;
    cfg.xlim= x_lim;
    %cfg.title = 'NREM';
    cfg.figure = 'gcf';
    phases(i).powspctrm= reshape(phases(i).powspctrm,1,39,86); %this way it can be plotted
    fig1_comps.p1=ft_singleplotTFR(cfg, phases(i));
    xline(0, '--k', 'LineWidth',2);
    c = colorbar;
    c.Label.String = 'Relative ITPC';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(titles(i), 'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');

    % Step onto the next tile
    index=index+1;
    ax(index) = nexttile;

    % Plot the T values in the second column

    % Initialize some settings
    if i==1
       stat=stat_p_t_itc;
       name=fieldnames(stat);
       toplot=stat.(name{1});
    elseif i==2
       stat=stat_p_w_itc; 
       name=fieldnames(stat);
       toplot=stat.(name{1});
    else
       stat=stat_t_w_itc;
       name=fieldnames(stat);
       toplot=stat.(name{1});
    end
    titles=["Phasic-Tonic", "Phasic-Wake", "Tonic-Wake"];

    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(squeeze(toplot.time), squeeze(toplot.freq), squeeze(toplot.stat));
    set(gca,'YDir','normal')
    ax=gca;
    ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xline(0, '--k', 'LineWidth',2)
    caxis([-2 2])
    xticks([0 0.2 0.4 0.6])
    c = colorbar;
    c.Label.String = 'T value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title(titles(i),'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');

    % Step onto the next tile
    index=index+1;
    ax(index) = nexttile;

    % Plot the p values and significant clusters in the third column
   % clusterstats=[stat_p_t_itc_cluster stat_p_w_itc_cluster stat_t_w_itc_cluster];

    % Initialize some settings
    if i==1
       clusterstat=stat_p_t_itc_cluster;
       name=fieldnames(clusterstat);
       toplotcluster=clusterstat.(name{1});
    elseif i==2
       clusterstat=stat_p_w_itc_cluster; 
       name=fieldnames(clusterstat);
       toplotcluster=clusterstat.(name{1});
    else
       clusterstat=stat_t_w_itc_cluster;
       name=fieldnames(clusterstat);
       toplotcluster=clusterstat.(name{1});
    end


    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(toplot.time, squeeze(toplot.freq), squeeze(toplot.prob));
    hold on
    colormap(ax(index),"bone")
    set(gca,'YDir','normal')
    ax=gca;
    ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xline(0, '--k', 'LineWidth', 2)
    caxis([0.001 0.05])
    xticks([0 0.2 0.4 0.6])
    c = colorbar;
    c.Label.String = 'p value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title( titles(i), 'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    clear fig1_comps
     hold on

    % increase the value of index
    index=index+1;
end
