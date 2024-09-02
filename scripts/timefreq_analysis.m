%% Time-frequency analysis on ANT channels
% Perform time-frequency analysis for specific time windows and frequencies

% Clear the command window, workspace, and close all figuresclc;
clc;
clear all;
close all;

% Set up the FieldTrip toolbox
ft_defaults

% List of patient names and sleep-wake phases
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

%sleep_wake_phases = {'phasic', 'tonic','wake'};
sleep_wake_phases = {'nrem'};


% Some initial settings
fs=512;

% Iterate through sleep-wake phases
for i = 1:length(sleep_wake_phases)
    count=0;
    PS = PLOT_STANDARDS();

    % Create a new figure for each phase
    figure;
    fig1_comps.fig = gcf;

    % Iterate through patients
    for j=1:length(patients)
        % Load the epoch data for the current patient and sleep-wake phase
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

        % Perform time-frequency analysis
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = channels; 
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.pad          = 1.2;         % only works with this value
        cfg.foi          = 7:1:45;     % analysis 7 to 100 Hz in steps of 1 Hz
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.15;  
        cfg.toi          = -0.125:0.01:0.725;   % time window "slides" from -0.125 to 0.725 sec in steps of 0.01 sec (10 ms)
        TFRhann = ft_freqanalysis(cfg, file);

        % Save the result into a new structure based on the sleep-wake phase
        if strcmpi(sleep_wake_phases(i),'phasic') == 1 
                phasic_timefreq{j}=TFRhann;
        elseif strcmpi(sleep_wake_phases(i),'tonic') == 1 
                tonic_timefreq{j}=TFRhann;
        elseif strcmpi(sleep_wake_phases(i),'wake') == 1 
               wake_timefreq{j}=TFRhann;
       elseif strcmpi(sleep_wake_phases(i),'nrem') == 1 
            nrem_timefreq{j}=TFRhann;
        end
         
        % Iterate through channels and plot time-frequency analysis results     
        % Plot
        for p=1:numberofchannels
            count=count+1;
            cfg = [];
            cfg.baseline     = [-0.2 -0.05];
            cfg.baselinetype = 'relative';
            cfg.figure = 'gcf';
            cfg.maskstyle    = 'saturation';
            cfg.zlim         = [0.8 1.2];
            cfg.xlim= [-0.2 0.6];
            cfg.ylim         = [7 45];
            cfg.channel      = char(channels{p});
            subplot(6,5,count)
            fig1_comps.p1=ft_singleplotTFR(cfg, TFRhann);
            fig1_comps.plotTitle=title(strcat(char(patients(j)), '-', char(sleep_wake_phases(i)),'-', channels{p}));
            fig1_comps.plotXLabel = xlabel('Time (sec)') ;
            fig1_comps.plotYLabel =ylabel({'Frequency', '(Hz)'});
            fig1_comps.plotZLabel =zlabel({'Power', '(dB)'});
            set([fig1_comps.plotXLabel, fig1_comps.plotYLabel, fig1_comps.plotZLabel], 'FontName', 'Times New Roman', 'FontSize', 8);
            set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 10, 'FontWeight' , 'bold');
            c = colorbar;
            c.Label.String = 'Power';
            c.Label.FontName= 'Times New Roman';
            c.Label.FontSize=8;
        end
     
        % Clear temporary variables
        clear file fns numberofchannels p k d f s segment channel1 channel2 channel3 channel4 channels;
    end
end

% Save the time-frequency analysis results
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz'
% save('tonic_timefreq_without_baseline_7_45', 'tonic_timefreq', '-v7.3');
% save('phasic_timefreq_without_baseline_7_45', 'phasic_timefreq', '-v7.3');
% save('wake_timefreq_without_baseline_7_45', 'wake_timefreq', '-v7.3');
save('nrem_timefreq_without_baseline_7_45', 'nrem_timefreq', '-v7.3');


% Create a structure that contains all the results
for z=1:length(patients)
    for t=1:length(sleep_wake_phases)
        if strcmpi(sleep_wake_phases(t),'phasic') == 1 
            patients_timefreq{z,t}=phasic_timefreq{z};
        elseif strcmpi(sleep_wake_phases(t),'tonic') == 1 
            patients_timefreq{z,t}=tonic_timefreq{z};
        elseif strcmpi(sleep_wake_phases(t),'wake') == 1 
            patients_timefreq{z,t}=wake_timefreq{z};
        end    
    end
end
% Save the time-frequency data
save('patients_timefreq_without_baseline_7_45', 'patients_timefreq', '-v7.3');

%% Baseline Correction

clear all

% load('f:\ANT_HEP\HEP_timefreq_512_70hz\patients_timefreq_without_baseline_7_45');
load('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_without_baseline_7_45');

%baselined_patients_timefreq=patients_timefreq; %the baselined data will be saved into this variable
baselined_patients_timefreq=nrem_timefreq; %the baselined data will be saved into this variable

% Perform baseline correction on the time-frequency data
for row=1:size(patients_timefreq,1)
    for column=1:size(patients_timefreq,2)
        cfg=[];
        cfg.baseline     = [-0.2 -0.05];
        cfg.baselinetype = 'relative';
        cfg.parameter = 'powspctrm';

        [freq] = ft_freqbaseline(cfg, patients_timefreq{row,column});
        baselined_patients_timefreq{row,column}=freq;
        clear freq
    end
end

% Save the baseline-corrected time-frequency data
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz';
save('nrem_timefreq_with_baseline_7_45', 'baselined_patients_timefreq', '-v7.3');


%% Manual grand average computation

%%%% Run only one of the sections below
% Open files WITHOUT baseline correction 
clear all
load('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_without_baseline_7_45');
% load('f:\ANT_HEP\HEP_timefreq_512_70hz\tonic_timefreq_without_baseline_7_45');
% load('f:\ANT_HEP\HEP_timefreq_512_70hz\phasic_timefreq_without_baseline_7_45');
% load('f:\ANT_HEP\HEP_timefreq_512_70hz\wake_timefreq_without_baseline_7_45');
%%%% Or
% Open files WITH baseline correction 
clear all
load('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_with_baseline_7_45');
% Divide the loaded data into new structures based on the sleep-wake phases
for row=1:size(baselined_patients_timefreq,1)
    for column=1:size(baselined_patients_timefreq,2)
        if column==1
            phasic_timefreq{row}=baselined_patients_timefreq{row,column};
        elseif column==2
            tonic_timefreq{row}=baselined_patients_timefreq{row,column};
        else
            wake_timefreq{row}=baselined_patients_timefreq{row,column};
        end
    end
end
clear column row
%%%%%


%%%% Working with the loaded variables
% Averaging across ANT channels
for column=1:size(baselined_patients_timefreq,2)
      baselined_patients_timefreq{1,column}.powspctrm = mean(baselined_patients_timefreq{1, column}.powspctrm,1);
%     phasic_timefreq{1,column}.powspctrm = mean(phasic_timefreq{1, column}.powspctrm,1);
%     tonic_timefreq{1,column}.powspctrm = mean(tonic_timefreq{1, column}.powspctrm,1);
%     wake_timefreq{1,column}.powspctrm = mean(wake_timefreq{1, column}.powspctrm,1);
end


% Averaging across patients
grandaverage_phasic_timefreq = zeros(1, size(phasic_timefreq{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(phasic_timefreq,2)
    grandaverage_phasic_timefreq = grandaverage_phasic_timefreq + phasic_timefreq{1, i}.powspctrm;
end
grandaverage_phasic_timefreq = grandaverage_phasic_timefreq / size(phasic_timefreq, 2); 
clear i

grandaverage_tonic_timefreq = zeros(1, size(tonic_timefreq{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(tonic_timefreq,2)
    grandaverage_tonic_timefreq = grandaverage_tonic_timefreq + tonic_timefreq{1, i}.powspctrm;
end
grandaverage_tonic_timefreq = grandaverage_tonic_timefreq / size(tonic_timefreq, 2); 
clear i

grandaverage_wake_timefreq = zeros(1, size(wake_timefreq{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(wake_timefreq,2)
    grandaverage_wake_timefreq = grandaverage_wake_timefreq + wake_timefreq{1, i}.powspctrm;
end
grandaverage_wake_timefreq = grandaverage_wake_timefreq / size(wake_timefreq, 2); 
clear i

grandaverage_nrem_timefreq = zeros(1, size(baselined_patients_timefreq{1, 1}.powspctrm, 2)); % pre-allocate
for i = 1:size(baselined_patients_timefreq,2)
    grandaverage_nrem_timefreq = grandaverage_nrem_timefreq + baselined_patients_timefreq{1, i}.powspctrm;
end
grandaverage_nrem_timefreq = grandaverage_nrem_timefreq / size(baselined_patients_timefreq, 2); 
clear i



%%%% Creating a structure for each phase for plotting later on
% this makes it easier to use the plotting functions (ft_singleplotTFR)

phasic=phasic_timefreq{1,1};
phasic.powspctrm=grandaverage_phasic_timefreq;
%label needs to be a 1x1 cell
temp=phasic_timefreq{1,1}.label;
phasic.label='';
phasic.label{1}=temp{1,1};

tonic=tonic_timefreq{1,1};
tonic.powspctrm=grandaverage_tonic_timefreq;
%label needs to be a 1x1 cell
temp=tonic_timefreq{1,1}.label;
tonic.label='';
tonic.label{1}=temp{1,1};

wake=wake_timefreq{1,1};
wake.powspctrm=grandaverage_wake_timefreq;
%label needs to be a 1x1 cell
temp=wake_timefreq{1,1}.label;
wake.label='';
wake.label{1}=temp{1,1};

nrem=baselined_patients_timefreq{1,1};
nrem.powspctrm=grandaverage_nrem_timefreq;
%label needs to be a 1x1 cell
temp=baselined_patients_timefreq{1,1}.label;
nrem.label='';
nrem.label{1}=temp{1,1};
% Save the grand averaged data
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz';
save('grand_average_nrem_timefreq_baselined_7_45', 'nrem', '-v7.3');
save('grand_average_phasic_timefreq_baselined_7_45', 'phasic', '-v7.3');
save('grand_average_tonic_timefreq_baselined_7_45', 'tonic', '-v7.3');
save('grand_average_wake_timefreq_baselined_7_45', 'wake', '-v7.3');


load('grand_average_phasic_timefreq_baselined_7_45')

%%%%% Plot

% Initialize some settings
%phases=[phasic tonic wake];
phases=[nrem];
titles=["NREM"];
% titles=["Phasic", "Tonic", "Wake"];
z_lim= [0.8 1.2];
x_lim=[-0.2 0.6];
font_size=10;
x_label='Time (sec)';
y_label='Frequency (Hz)';

% Set up figure standards
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

% Loop through phases and create subplots
for i=1:length(phases)
    subplot(1,1,i)
    cfg = [];
    % Pay attention to cfg.baseline, this is optional to use
    % cfg.baseline     = [-0.2 -0.05];
    %  cfg.baselinetype = 'relative';
    cfg.zlim=z_lim;
    cfg.xlim= x_lim;
    cfg.title = 'Phasic';
    cfg.colormap ='jet';
    cfg.figure = 'gcf';
    fig1_comps.p1=ft_singleplotTFR(cfg, phases(i));
    xline(0, '--w', 'LineWidth',3);
    c = colorbar;
    c.Label.String = 'Power';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(titles(i), 'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    hold on
end

    
%% Statistical Analysis

%%%% Run only one of the sections below
% Load the data without baseline correction
clear all
load 'd:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_without_baseline_7_45';
load 'd:\ANT_HEP\HEP_timefreq_512_70hz\patients_timefreq_without_baseline_7_45';

%%%% Or
% Load the data with baseline correction
clear all
rem=load ('d:\ANT_HEP\HEP_timefreq_512_70hz\patients_timefreq_with_baseline_7_45');
nrem=load ('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_with_baseline_7_45');
%patients_timefreq=baselined_patients_timefreq;
%%%%
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

for i=1:length(patients)
    rem.baselined_patients_timefreq(i,4)=nrem.baselined_patients_timefreq(1,i);
end
loop_version = {'phasic', 'tonic','wake','nrem'};
for j=1:length(loop_version)
    for i=1:length(patients)
        baselined_patients_timefreq(i,j)=rem.baselined_patients_timefreq(i,j);
    end
end


%%%% Continue with loaded variables
% Average the ANT channels
for row=1:size(baselined_patients_timefreq,1)
    for column=1:size(baselined_patients_timefreq,2)
        baselined_patients_timefreq{row,column}.powspctrm = mean(baselined_patients_timefreq{row,column}.powspctrm,1);
        baselined_patients_timefreq{row,column}.label={'General-ANT'};
    end
end
clear column row

clear nrem
% Divide the data into new structures based on the sleep-wake phases
for newrow=1:size(baselined_patients_timefreq,1)
    phasic{newrow}=baselined_patients_timefreq{newrow,1};
    tonic{newrow}=baselined_patients_timefreq{newrow,2};
    wake{newrow}=baselined_patients_timefreq{newrow,3};
    nrem{newrow}=baselined_patients_timefreq{newrow,4};
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
cfg.parameter='powspctrm';
cfg.method = 'montecarlo'; 
cfg.statistic = 'ft_statfun_depsamplesT';  
cfg.avgoverfreq='no';
cfg.latency=[0.05 0.65];
cfg.avgovertime='no';
cfg.numrandomization = 'all'; 
cfg.tail = 0; 
cfg.alpha = 0.05; 
cfg.uvar = 1; 
cfg.ivar = 2; 

% Perform frequency statistics
[stat_p_n_timefreq] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
[stat_t_n_timefreq] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
[stat_w_n_timefreq] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
% [stat_p_t_timefreq] = ft_freqstatistics(cfg, phasic{1,:}, tonic{1,:}); 
% [stat_p_w_timefreq] = ft_freqstatistics(cfg,  phasic{1,:}, wake{1,:}); 
% [stat_t_w_timefreq] = ft_freqstatistics(cfg,  tonic{1,:}, wake{1,:}); 

%Save the results
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz'
save('stat_p_n_timefreq_with_baseline_7_45', 'stat_p_n_timefreq', '-v7.3');
save('stat_t_n_timefreq_with_baseline_7_45', 'stat_t_n_timefreq', '-v7.3');
save('stat_w_n_timefreq_with_baseline_7_45', 'stat_w_n_timefreq', '-v7.3');
% save('stat_p_t_timefreq_without_baseline_7_45', 'stat_p_t_timefreq', '-v7.3');
% save('stat_p_w_timefreq_without_baseline_7_45', 'stat_p_w_timefreq', '-v7.3');
% save('stat_t_w_timefreq_without_baseline_7_45', 'stat_t_w_timefreq', '-v7.3');


%%%% Plot

% Initialize some settings
font_size = 10;
x_label = 'Time (sec)';
y_label = 'Frequency (Hz)';
stat=[stat_p_n_timefreq stat_t_n_timefreq stat_w_n_timefreq];
% stat=[stat_p_t_timefreq stat_p_w_timefreq stat_t_w_timefreq];
titles=["Phasic-NREM", "Tonic-NREM", "Wake-NREM"];

% Plot prob, p values
PS = PLOT_STANDARDS();
figure;
for i=1:3
    fig1_comps.fig = gcf;
    ax1=subplot(3, 1, i);
    fig1_comps.p1=imagesc(squeeze(stat(i).time), squeeze(stat(i).freq), squeeze(stat(i).prob));
    set(gca,'YDir','normal')
    xline(0, '--k', 'LineWidth', 2)
    colormap(ax1,bone)
    caxis([0.001 0.05])
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
end

hold off

% Plot stat, T values
PS = PLOT_STANDARDS();
figure;
for i=1:3
    fig1_comps.fig = gcf;
    ax2=subplot(3,1,i);
    fig1_comps.p1=imagesc(squeeze(stat(i).time), squeeze(stat(i).freq), squeeze(stat(i).stat));
    set(gca,'YDir','normal')
    xline(0, '--w', 'LineWidth',2)
    caxis([-2 2])
    colormap(ax2,"jet")
    c = colorbar;
    c.Label.String = 'T value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title(titles(i),'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end
hold off


%% Cluster Analysis

%%%%% Run only one of the sections below
% Open files WITHOUT baseline correction
clear all
load 'f:\ANT_HEP\HEP_timefreq_512_70hz\patients_timefreq_without_baseline_7_45';
%%% Or
% Load the data with baseline correction
clear all
rem=load ('d:\ANT_HEP\HEP_timefreq_512_70hz\patients_timefreq_with_baseline_7_45');
nrem=load ('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_with_baseline_7_45');
%patients_timefreq=baselined_patients_timefreq;
%%%%
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 

for i=1:length(patients)
    rem.baselined_patients_timefreq(i,4)=nrem.baselined_patients_timefreq(1,i);
end
loop_version = {'phasic', 'tonic','wake','nrem'};
for j=1:length(loop_version)
    for i=1:length(patients)
        baselined_patients_timefreq(i,j)=rem.baselined_patients_timefreq(i,j);
    end
end



%%%% Continue with loaded variables
% Average the ANT channels
for row=1:size(baselined_patients_timefreq,1)
    for column=1:size(baselined_patients_timefreq,2)
        baselined_patients_timefreq{row,column}.powspctrm = mean(baselined_patients_timefreq{row,column}.powspctrm,1);
        baselined_patients_timefreq{row,column}.label={'General-ANT'};
    end
end
clear column row

clear nrem
% Divide the data into new structures based on the sleep-wake phases
for newrow=1:size(baselined_patients_timefreq,1)
    phasic{newrow}=baselined_patients_timefreq{newrow,1};
    tonic{newrow}=baselined_patients_timefreq{newrow,2};
    wake{newrow}=baselined_patients_timefreq{newrow,3};
   % nrem{newrow}=baselined_patients_timefreq{newrow,4};
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
cfg.parameter='powspctrm';
cfg.method = 'montecarlo'; 
cfg.statistic = 'ft_statfun_depsamplesT';  
cfg.avgoverfreq='no';
cfg.avgovertime='no';
cfg.correctm = 'cluster'; 
cfg.clusterstatistic = 'maxsum'; 
cfg.clusteralpha = 0.05; 
cfg.numrandomization = 'all'; 
cfg.latency=[0.05 0.65];
cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
cfg.uvar = 1; 
cfg.ivar = 2; 
cfg.numrandomization = 'all'; 

% Perform cluster-based statistics
%[stat_p_n_timefreq_cluster] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
%[stat_t_n_timefreq_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
%[stat_w_n_timefreq_cluster] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
 [stat_p_t_timefreq_cluster] = ft_freqstatistics(cfg, phasic{1,:}, tonic{1,:}); 
 [stat_p_w_timefreq_cluster] = ft_freqstatistics(cfg,  phasic{1,:}, wake{1,:}); 
[stat_t_w_timefreq_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, wake{1,:}); 
% % 
% % Save the results
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz'
save('stat_p_n_timefreq_cluster_with_baseline_7_45', 'stat_p_n_timefreq_cluster', '-v7.3');
save('stat_t_n_timefreq_cluster_with_baseline_7_45', 'stat_t_n_timefreq_cluster', '-v7.3');
save('stat_w_n_timefreq_cluster_with_baseline_7_45', 'stat_w_n_timefreq_cluster', '-v7.3');
% save('stat_p_t_timefreq_cluster_without_baseline_7_45', 'stat_p_t_timefreq_cluster', '-v7.3');
% save('stat_p_w_timefreq_cluster_without_baseline_7_45', 'stat_p_w_timefreq_cluster', '-v7.3');
% save('stat_t_w_timefreq_cluster_without_baseline_7_45', 'stat_t_w_timefreq_cluster', '-v7.3');
% % 

%%%%% Select clusters for the plot

% Select clusters for p_t: negative cluster 
p_t_neg_cluster_pvals = [stat_p_t_timefreq_cluster.negclusters(:).prob];
p_t_neg_clust         = find(p_t_neg_cluster_pvals < 0.05);
stat_p_t_timefreq_cluster.negclusterslabelmat=squeeze(stat_p_t_timefreq_cluster.negclusterslabelmat);
p_t_clustermatrix = stat_p_t_timefreq_cluster.negclusterslabelmat == p_t_neg_clust;

% Select clusters for p_w: negative cluster 
p_w_neg_cluster_pvals = [stat_p_w_timefreq_cluster.negclusters(:).prob];
p_w_neg_clust         = find(p_w_neg_cluster_pvals < 0.05);
stat_p_w_timefreq_cluster.negclusterslabelmat=squeeze(stat_p_w_timefreq_cluster.negclusterslabelmat);
p_w_clustermatrix_1 = stat_p_w_timefreq_cluster.negclusterslabelmat == p_w_neg_clust(1,1);
p_w_clustermatrix_2 = stat_p_w_timefreq_cluster.negclusterslabelmat == p_w_neg_clust(1,2);

% Select clusters for t_w: positive cluster
t_w_pos_cluster_pvals = [stat_t_w_timefreq_cluster.posclusters(:).prob];
t_w_pos_clust         = find(t_w_pos_cluster_pvals < 0.05);
stat_t_w_timefreq_cluster.posclusterslabelmat=squeeze(stat_t_w_timefreq_cluster.posclusterslabelmat);
t_w_clustermatrix = stat_t_w_timefreq_cluster.posclusterslabelmat == t_w_pos_clust;

% Select clusters for p_n: negative cluster
p_n_neg_cluster_pvals = [stat_p_n_timefreq_cluster.negclusters(:).prob];
p_n_neg_clust         = find(p_n_neg_cluster_pvals < 0.05);
stat_p_n_timefreq_cluster.negclusterslabelmat=squeeze(stat_p_n_timefreq_cluster.negclusterslabelmat);
p_n_clustermatrix = stat_p_n_timefreq_cluster.negclusterslabelmat == p_n_neg_clust;


%%%%% Perform operations on selected clusters
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz'
load 'grand_average_nrem_timefreq_baselined_7_45.mat'
load 'grand_average_phasic_timefreq_baselined_7_45.mat'
load 'grand_average_tonic_timefreq_baselined_7_45.mat'
load 'grand_average_wake_timefreq_baselined_7_45.mat'
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'powspctrm';
% Subtract phasic and tonic grand averaged data
grandaverage_phasic_vs_tonic = ft_math(cfg,phasic,tonic);
stat_p_t_timefreq_cluster.raw = (grandaverage_phasic_vs_tonic.powspctrm);
% Subtract phasic and wake grand averaged data
grandaverage_phasic_vs_wake = ft_math(cfg,phasic,wake);
stat_p_w_timefreq_cluster.raw = (grandaverage_phasic_vs_wake.powspctrm);
% Subtract tonic and wake grand averaged data
grandaverage_tonic_vs_wake = ft_math(cfg,tonic,wake);
stat_t_w_timefreq_cluster.raw = (grandaverage_tonic_vs_wake.powspctrm);
% Subtract phasic and nrem grand averaged data
grandaverage_phasic_vs_nrem = ft_math(cfg,phasic,nrem);
stat_p_n_timefreq_cluster.raw = (grandaverage_phasic_vs_nrem.powspctrm);
% Subtract tonic and nrem grand averaged data
grandaverage_tonic_vs_nrem = ft_math(cfg,tonic,nrem);
stat_t_n_timefreq_cluster.raw = (grandaverage_tonic_vs_nrem.powspctrm);
% Subtract wake and nrem grand averaged data
grandaverage_wake_vs_nrem = ft_math(cfg,wake,nrem);
stat_w_n_timefreq_cluster.raw = (grandaverage_wake_vs_nrem.powspctrm);

% Load significant differences from the previous analysis
cd 'd:\ANT_HEP\HEP_timefreq_512_70hz'
% load('stat_p_t_timefreq_with_baseline_7_45');
% load('stat_p_w_timefreq_with_baseline_7_45');
% load('stat_t_w_timefreq_with_baseline_7_45');
load('stat_p_n_timefreq_with_baseline_7_45');
load('stat_t_n_timefreq_with_baseline_7_45');
load('stat_w_n_timefreq_with_baseline_7_45');
% For phasic-nrem
pn_pvals = squeeze(stat_p_n_timefreq.prob);
pn_loc         = pn_pvals < 0.05 ;
simple_boundary_pn = bwboundaries(pn_loc);
% For phasic-tonic
pt_pvals = squeeze(stat_p_t_timefreq.prob);
pt_loc         = pt_pvals < 0.05 ;
simple_boundary_pt = bwboundaries(pt_loc);
% For phasic-wake
pw_pvals = squeeze(stat_p_w_timefreq.prob);
pw_loc         = pw_pvals < 0.05 ;
simple_boundary_pw = bwboundaries(pw_loc);
% For tonic-wake
tw_pvals = squeeze(stat_t_w_timefreq.prob);
tw_loc         = tw_pvals < 0.05 ;
simple_boundary_tw = bwboundaries(tw_loc);

%%%%% Plot the results

% Get outlines of clusters
boundary_pn = bwboundaries(p_n_clustermatrix);

boundary_pt = bwboundaries(p_t_clustermatrix);
boundary_pw = bwboundaries(p_w_clustermatrix_1);
boundary_pw_2 = bwboundaries(p_w_clustermatrix_2);
boundary_tw = bwboundaries(t_w_clustermatrix);

% Initialize some settings
font_size=10;
x_label='Time (sec)';
y_label='Frequency (Hz)';

averages=[grandaverage_phasic_vs_nrem grandaverage_tonic_vs_nrem grandaverage_wake_vs_nrem];
% averages=[grandaverage_phasic_vs_tonic grandaverage_phasic_vs_wake grandaverage_tonic_vs_wake];
titles=["Phasic-NREM", "Tonic-NREM", "Wake-NREM"];
%cluster_boundaries=[boundary_pn];
cluster_boundaries=[boundary_pt boundary_pw boundary_pw_2];
%simple_boundaries={simple_boundary_pn};
% simple_boundaries={simple_boundary_pt simple_boundary_pw simple_boundary_tw};

PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;
for i=1:3
    subplot(3,1,i)
    cfg = [];
    cfg.colormap ='jet';
    cfg.figure = 'gcf';
    fig1_comps.p1=ft_singleplotTFR(cfg, averages(i));
    xline(0, '--w', 'LineWidth',3);
    caxis([-1 1])
    xlim([-0.1 0.6])
    c = colorbar;
    c.Label.String = 'Power';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=font_size;
    fig1_comps.plotTitle=title(titles(i), 'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    
    hold on; % Allow the plot to overlay on the image
    
    % Plot the significant differences as a shaded area
    % Loop through boundaries and shade the areas
    for k=1:length(simple_boundaries{1,i})
        d=simple_boundaries{1,i}{k,1};
        fill(averages(i).time(d(:,2)), averages(i).freq(d(:,1)), 'black', 'FaceAlpha', 0.3);
    end

    % Plot the most significant clusters boundary as a line
    b = cluster_boundaries{1,i};
    plot(averages(i).time(b(:,2)), averages(i).freq(b(:,1)), 'w', 'LineWidth', 2);
    
    hold off
end


%% Segmented statistical analysis

% Open files WITHOUT baseline correction 
clear all
load('f:\ANT_HEP\HEP_timefreq_512_70hz\tonic_timefreq_without_baseline_7_45');
load('f:\ANT_HEP\HEP_timefreq_512_70hz\phasic_timefreq_without_baseline_7_45');
load('f:\ANT_HEP\HEP_timefreq_512_70hz\wake_timefreq_without_baseline_7_45');
load('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_without_baseline_7_45');


patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 
sleep_wake_phases = {'nrem'};
% sleep_wake_phases = {'phasic', 'tonic','wake'};

% Save into a big structure
% for row=1:length(sleep_wake_phases)
%     for column=1:length(patients)
%         if row==1
%             patients_timefreq{1,column} =phasic_timefreq{1,column};
%         elseif row==2
%             patients_timefreq{2,column}= tonic_timefreq{1, column};
%         else
%             patients_timefreq{3,column}= wake_timefreq{1, column};
%         end
%     end
% end
% Save into a big structure
for row=1:length(sleep_wake_phases)
    for column=1:length(patients)
        if row==1
            patients_timefreq{1,column} =nrem_timefreq{1,column};
        end
    end
end
clear column row

% Average the ANT channels
for row=1:size(patients_timefreq,1)
    for column=1:size(patients_timefreq,2)
        patients_timefreq{row,column}.powspctrm = mean(patients_timefreq{row,column}.powspctrm,1);
        patients_timefreq{row,column}.label={'General-ANT'};
    end
end
clear column row


% Segment into 75 ms segments
baseline_indices=1:8; % baseline: -125 ms - -50ms
start=19; % first segment starts at 50ms
% Creating the segments
for i=1:floor((length(patients_timefreq{1,1}.time)-start)/length(baseline_indices))
    timepoints=patients_timefreq{1,1}.time;
    segments{i}=timepoints(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
    segment_indices{i}=(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
end
baseline_window=timepoints(baseline_indices);


% Run the statistical analysis on each segment
for i=1:length(segments)

    % Create the structure for each segment(i)

    for j=1:length(sleep_wake_phases)
        % Cut every patient's data into segments
        % (new variable: phase x patient)
        for k=1:length(patients)
            act_segment{j,k}.label=patients_timefreq{j,k}.label;
            act_segment{j,k}.dimord=patients_timefreq{j,k}.dimord;
            act_segment{j,k}.freq=patients_timefreq{j,k}.freq;
            act_segment{j,k}.cfg=patients_timefreq{j,k}.cfg;
            act_segment{j,k}.time=patients_timefreq{j,k}.time(1,segment_indices{i});
            act_segment{j,k}.powspctrm=patients_timefreq{j,k}.powspctrm(:,:,segment_indices{i});
        
            bl_segment{j,k}.label=patients_timefreq{j,k}.label;
            bl_segment{j,k}.dimord=patients_timefreq{j,k}.dimord;
            bl_segment{j,k}.freq=patients_timefreq{j,k}.freq;
            bl_segment{j,k}.cfg=patients_timefreq{j,k}.cfg;
            bl_segment{j,k}.time=patients_timefreq{j,k}.time(1,segment_indices{i}); %using this time interval, because the statistical analysis runs only with the same time points
            bl_segment{j,k}.powspctrm=patients_timefreq{j,k}.powspctrm(:,:,baseline_indices);
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
    cfg.avgoverfreq='no';
    cfg.avgovertime='no';
    cfg.numrandomization = 'all'; 
%     cfg.correctm = 'cluster'; 
%     cfg.clusterstatistic = 'maxsum'; 
%     cfg.clusteralpha = 0.05; 
    cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
    cfg.alpha = 0.05; 
    cfg.uvar = 1; 
    cfg.ivar = 2; 

    % Run the analysis
     stat_n_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
%     stat_p_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
%     stat_t_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{2,:}, bl_segment{2,:}); 
%     stat_w_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{3,:}, bl_segment{3,:}); 

end

%Concatenate the segments 
% stats=[stat_p_timefreq_seg; stat_t_timefreq_seg; stat_w_timefreq_seg];
stats=[stat_n_timefreq_seg];


unified{1,1}=stats{1,1};
unified{2,1}=stats{2,1};
unified{3,1}=stats{3,1};

for i = 1:size(stats,1) 
   unified{i,1}.time=cat(2, stats{i,1}.time, stats{i,2}.time, stats{i,3}.time, stats{i,4}.time, stats{i,5}.time, stats{i,6}.time, stats{i,7}.time, stats{i,8}.time);
   unified{i,1}.prob=cat(3, stats{i,1}.prob, stats{i,2}.prob, stats{i,3}.prob, stats{i,4}.prob, stats{i,5}.prob, stats{i,6}.prob, stats{i,7}.prob, stats{i,8}.prob);
   unified{i,1}.stat=cat(3, stats{i,1}.stat, stats{i,2}.stat, stats{i,3}.stat, stats{i,4}.stat, stats{i,5}.stat, stats{i,6}.stat, stats{i,7}.stat, stats{i,8}.stat);
   unified{i,1}.ref=cat(3, stats{i,1}.ref, stats{i,2}.ref, stats{i,3}.ref, stats{i,4}.ref, stats{i,5}.ref, stats{i,6}.ref, stats{i,7}.ref, stats{i,8}.ref);
   unified{i,1}.cirange=cat(3, stats{i,1}.cirange, stats{i,2}.cirange, stats{i,3}.cirange, stats{i,4}.cirange, stats{i,5}.cirange, stats{i,6}.cirange, stats{i,7}.cirange, stats{i,8}.cirange);
   unified{i,1}.mask=cat(3, stats{i,1}.mask, stats{i,2}.mask, stats{i,3}.mask, stats{i,4}.mask, stats{i,5}.mask, stats{i,6}.mask, stats{i,7}.mask, stats{i,8}.mask);
end
%this structure can be now plotted

% Plot the results

% Initialize some settings
font_size = 10;
x_label = 'Time (sec)';
y_label = 'Frequency (Hz)';
stat=[unified{1,1}];
% stat=[unified{1,1}; unified{2,1}; unified{3,1}];
titles=["NREM - baseline"];
% titles=["Phasic - baseline", "Tonic - baseline", "Wake - baseline"];



% Plot prob, p values
PS = PLOT_STANDARDS();
figure;
for i=1:length(titles)
    for j=1:length(segments)
        fig1_comps.fig = gcf;
        subplot(3, 1, i)
        fig1_comps.p1=imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).prob));
        hold on
        set(gca,'YDir','normal')
        xline(0, '--k', 'LineWidth', 2)
        colormap(bone)
        caxis([0.001 0.05])
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
    end
end


% Plot stat, T values
PS = PLOT_STANDARDS();
figure;
for i=1:length(titles)
    fig1_comps.fig = gcf;
    subplot(3,1,i)
    fig1_comps.p1=imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).stat));
    set(gca,'YDir','normal')
    xline(0, '--w', 'LineWidth',2)
    caxis([-2 2])
    c = colorbar;
    c.Label.String = 'T value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title(titles(i),'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end



%% Plot that contains every result from time-frequency analysis 

clear all
%nonbaselined
    % open the grand averaged data
    cd 'f:\ANT_HEP\HEP_timefreq_512_70hz';
    load('grand_average_phasic_timefreq_nonbaselined_7_45');
    load('grand_average_tonic_timefreq_nonbaselined_7_45');
    load('grand_average_wake_timefreq_nonbaselined_7_45');
    
    % open the statistics
    load('stat_p_t_timefreq_without_baseline_7_45');
    load('stat_p_w_timefreq_without_baseline_7_45');
    load('stat_t_w_timefreq_without_baseline_7_45');
    
    % open the cluster based statistics
    load('stat_p_t_timefreq_cluster_without_baseline_7_45');
    load('stat_p_w_timefreq_cluster_without_baseline_7_45');
    load('stat_t_w_timefreq_cluster_without_baseline_7_45');

clear all
%baselined
    % open the grand averaged data
    cd 'd:\ANT_HEP\HEP_timefreq_512_70hz';
  %  load('grand_average_nrem_timefreq_baselined_7_45');
    load('grand_average_phasic_timefreq_baselined_7_45');
    load('grand_average_tonic_timefreq_baselined_7_45');
    load('grand_average_wake_timefreq_baselined_7_45');
    
    % open the statistics
   % load('stat_p_n_timefreq_with_baseline_7_45');
    %load('stat_t_n_timefreq_with_baseline_7_45');
    %load('stat_w_n_timefreq_with_baseline_7_45');
      
   load('stat_p_t_timefreq_with_baseline_7_45');
     load('stat_p_w_timefreq_with_baseline_7_45');
     load('stat_t_w_timefreq_with_baseline_7_45');
%     
    % open the cluster based statistics
  %  load('stat_p_n_timefreq_cluster_with_baseline_7_45');
   % load('stat_t_n_timefreq_cluster_with_baseline_7_45');
    %load('stat_w_n_timefreq_cluster_with_baseline_7_45');

 load('stat_p_t_timefreq_cluster_with_baseline_7_45');
     load('stat_p_w_timefreq_cluster_with_baseline_7_45');
     load('stat_t_w_timefreq_cluster_with_baseline_7_45');
% 

% Plot

tcl = tiledlayout(3,3);  
ax = gobjects(1,9);
index=1;
for i = 1:3
     ax(index) = nexttile;
    
    % Plot the grand averaged data in the first column
    % Initialize some settings
    phases=[phasic tonic wake];
%    phases=[nrem nrem nrem];
    titles=["Phasic", "Tonic", "Wake"];
    z_lim= [0.8 1.2];
    x_lim=[-0.13 0.6];
    font_size=10;
    x_label='Time (sec)';
    y_label='Frequency (Hz)';
    
    % Set up figure standards
    PS = PLOT_STANDARDS();
    fig1_comps.fig = gcf;
    colormap(ax(index), "jet")  % STEP 2
    cfg = [];
    % Pay attention to cfg.baseline, this is optional to use
    % cfg.baseline     = [-0.2 -0.05];
    %  cfg.baselinetype = 'relative';
    cfg.zlim=z_lim;
    cfg.xlim= x_lim;
    cfg.figure = 'gcf';
    fig1_comps.p1=ft_singleplotTFR(cfg, phases(i));
    xline(0, '--k', 'LineWidth',2);
    c = colorbar;
    c.Label.String = 'Relatív erő';
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
    stat=[stat_p_t_timefreq stat_p_w_timefreq stat_t_w_timefreq];
    titles=["Phasic-Tonic", "Phasic-Wake", "Tonic-Wake"];

    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(squeeze(stat(i).time), squeeze(stat(i).freq), squeeze(stat(i).stat));
    colormap(ax(index), "jet")  
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

    % Plot the p values in the third column

    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).prob));
    colormap(ax(index),"bone")
    set(gca,'YDir','normal')
    ax=gca;
    ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xticks([0 0.2 0.4 0.6])
    xline(0, '--k', 'LineWidth',2)
    caxis([0.001 0.05])
    c = colorbar;
    c.Label.String = 'p value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title({strcat(titles(i))},'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');


    % Step onto the next tile
%    index=index+1;
%      ax(index) = nexttile;
    
    %Plot the significant clusters in the fourth column

    %Initialize some settings    
    clusterstats=[stat_p_t_timefreq_cluster stat_p_w_timefreq_cluster stat_t_w_timefreq_cluster];
    titles=["Phasic-Tonic - cluster", "Phasic-Wake - cluster", "Tonic-Wake - cluster"];

%     % Set up figure standards
%     fig1_comps.fig = gcf;
%     fig1_comps.p1=imagesc(clusterstats(i).time, squeeze(clusterstats(i).freq), squeeze(clusterstats(i).prob));
%     hold on
%     colormap(ax(index),"bone")
%     set(gca,'YDir','normal')
%     ax=gca;
%     ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column
%     xline(0, '--k', 'LineWidth', 2)
%     caxis([0.001 0.05])
%     c = colorbar;
%     c.Label.String = 'p value';
%     c.Label.FontName= 'Times New Roman';
%     c.Label.FontSize=font_size;
%     fig1_comps.plotTitle=title( titles(i), 'fontsize', font_size);
%     fig1_comps.plotXLabel = xlabel(x_label) ;
%     fig1_comps.plotYLabel =ylabel(y_label);
%     set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
%     set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
%     clear fig1_comps

    % Plot the significant clusters inside
    hold on
%    clusterstats=[stat_p_n_timefreq_cluster stat_t_n_timefreq_cluster stat_w_n_timefreq_cluster];
    clusterstats=[stat_p_t_timefreq_cluster stat_p_w_timefreq_cluster stat_t_w_timefreq_cluster];

    simple_boundary = bwboundaries(squeeze(clusterstats(i).prob)<0.05);
    if size(simple_boundary,1) >0
        if size(simple_boundary,1) ==1
            for b=1:length(simple_boundary{1,1})
                        d=simple_boundary{1,1};
                        c = [0.8500 0.3250 0.0980]; %colour code
                        h1=plot(stat(i).time(d(:,2)), stat(i).freq(d(:,1)), 'Color', c, 'LineWidth', 2); %FaceAlpha', 0.3)
                    set(h1, 'DisplayName', 'significant cluster');
                            fig1_comps.plotLegend=legend(h1,'significant cluster', "Location","southoutside", "FontName", 'Times New Roman');


            end
        end
         if size(simple_boundary,1) ==2
            for b=1:length(simple_boundary{1,1})
                        d=simple_boundary{1,1};
                        c = [0.8500 0.3250 0.0980]; %colour code
                        h1=plot(stat(i).time(d(:,2)), stat(i).freq(d(:,1)), 'Color', c, 'LineWidth', 2); %FaceAlpha', 0.3)
                    set(h1, 'DisplayName', 'significant cluster');
                            fig1_comps.plotLegend=legend(h1,'significant clusters', "Location","southoutside", "FontName", 'Times New Roman');


            end
             for b=1:length(simple_boundary{2,1})
                        d=simple_boundary{2,1};
                        c = [0.8500 0.3250 0.0980]; %colour code
                        h2=plot(stat(i).time(d(:,2)), stat(i).freq(d(:,1)), 'Color', c, 'LineWidth', 2); %FaceAlpha', 0.3)
                    set(h2, 'DisplayName', 'significant cluster');
             fig1_comps.plotLegend=legend(h2,'significant clusters', "Location","southoutside", "FontName", 'Times New Roman');

             end
        end
        clear h1 h2 d
    end
%   
end


%% Rügers area
%whether at least half of these results were signiﬁcant at least at 1/2 of the conventional p = 0.05 signiﬁcance level
%at least one-third of them were signiﬁcant at least at 1/3 of the conventional p = 0.05 signiﬁcance level
%If both of these conditions were fulﬁlled, the area as a whole was considered signiﬁcant.

clear all
% open the statistics
cd 'f:\ANT_HEP\HEP_timefreq_512_70hz';
load('stat_p_t_timefreq_with_baseline_7_45');
load('stat_p_w_timefreq_with_baseline_7_45');
load('stat_t_w_timefreq_with_baseline_7_45');


stats=[stat_p_t_timefreq stat_p_w_timefreq stat_t_w_timefreq];

for s=1:size(stats,2)
    %choose the area where the p values are under 0.05
    areas = bwboundaries(squeeze(stats(s).prob)<0.05);
    
    %loop through the blobs
    for a=1:size(areas,1)
        d=areas{a,1};   
        %create the time-freq matrix, with the blob outlines being 1, everything
        %else 0
        M=size(stat_p_t_timefreq.time,2);
        N=size(stat_p_t_timefreq.freq,2);
        tfmatrix=zeros(N,M);
        %fill the blobs with 1s
        for i=1:N
            for j=1:M
                for k=1:size(d,1)
                     if i==d(k,1) && j==d(k,2)
                        tfmatrix(i,j)=1;
                     end
                end
           end
        end
    
        %search for the prob values in the places of ones
        filled_tfmatrix=imfill(tfmatrix); %basically this is the mask
        %create an inverse for this
        inverse_filled_tfmatrix=1-filled_tfmatrix;
        for i=1:N
            for j=1:M
                if inverse_filled_tfmatrix(i,j)==0
                    probs=squeeze(stats(s).prob);
                    prob_filled_tfmatrix(i,j)=probs(i,j);
                else
                    prob_filled_tfmatrix(i,j)=inverse_filled_tfmatrix(i,j);
                end
            end
        end
    
        %is this a significant rügers area?
        allnumbers(a,1)=0;
        allnumber=0;
        half05(a,1)=0;
        half=0;
        third05(a,1)=0;
        third=0;
        for i=1:N
            for j=1:M
                if prob_filled_tfmatrix(i,j)<0.05 
                    allnumber=allnumber+1;
                end
                if prob_filled_tfmatrix(i,j)<0.025
                    half=half+1;
                end
                if prob_filled_tfmatrix(i,j)<0.0167
                    third=third+1;
                end
            end
        end
        stat(s).allnumbers(a,1)=allnumber;
        stat(s).half05(a,1)=half;
        stat(s).third05(a,1)=third;
        if stat(s).half05(a,1)>=stat(s).allnumbers(a,1)/2 && stat(s).third05(a,1)>=stat(s).allnumbers(a,1)/3 && stat(s).allnumbers(a,1)~=0 && stat(s).half05(a,1)~=0 && stat(s).third05(a,1)~=0
            stat(s).sig(a)=1;
        else
           stat(s).sig(a)=0;
        end
        %save the masked matrix for each found area  
        stat(s).maskedmatrix{a}=filled_tfmatrix;
        clear prob_filled_tfmatrix inverse_filled_tfmatrix filled_tfmatrix
    end
end

%Create a matrix that only contains the significant Rüger area
%probabilities
%create the mask
for i=1:3
    stat(i).sigmaskedmatrix=zeros(N,M);
    for j=1:length(stat(i).sig)
        if  stat(i).sig(1,j)==1
            stat(i).sigmaskedmatrix=stat(i).sigmaskedmatrix+stat(i).maskedmatrix{1,j};
        end
    end
end

%fill it with p values
for i=1:3
    for j=1:size(stat(i).sigmaskedmatrix,1)
        for k=1:size(stat(i).sigmaskedmatrix,2)
             if stat(i).sigmaskedmatrix(j,k)==1
                probs=squeeze(stats(i).prob);
                stat(i).sigprobmatrix(j,k)=probs(j,k);
             else
                stat(i).sigprobmatrix(j,k)=1;             
             end
        end
    end
end


%Plot 
    figure

for i = 1:3
    subplot(3,1,i)
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(stats(i).time, squeeze(stats(i).freq), stat(i).sigprobmatrix);
    colormap("bone")
    set(gca,'YDir','normal')
    ax=gca;
    %ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xticks([0 0.2 0.4 0.6])
    xline(0, '--k', 'LineWidth',2)
    caxis([0.0001 0.05])
    c = colorbar;
    c.Label.String = 'p érték';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title('title', 'fontsize', 10);
    fig1_comps.plotXLabel = xlabel('time') ;
    fig1_comps.plotYLabel =ylabel('freq');
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 10, 'FontWeight' , 'bold');
end



%% Plot that contains every result from time-frequency analysis with Rügers area

%baselined
    % open the grand averaged data
    cd 'f:\ANT_HEP\HEP_timefreq_512_70hz';
    load('grand_average_phasic_timefreq_baselined_7_45');
    load('grand_average_tonic_timefreq_baselined_7_45');
    load('grand_average_wake_timefreq_baselined_7_45');

% Plot
figure
tcl = tiledlayout(3,3);  
ax = gobjects(1,9);
index=1;
for i = 1:3
     ax(index) = nexttile;
    
    % Plot the grand averaged data in the first column
    % Initialize some settings
    phases=[phasic tonic wake];
    titles=["Fázisos", "Tónusos", "Ébrenlét"];
    z_lim= [0.8 1.2];
    x_lim=[-0.13 0.6];
    font_size=10;
    x_label='Idő (mp)';
    y_label='Frekvencia (Hz)';
    
    % Set up figure standards
    PS = PLOT_STANDARDS();
    fig1_comps.fig = gcf;
    colormap(ax(index), "jet")  % STEP 2
    cfg = [];
    % Pay attention to cfg.baseline, this is optional to use
    % cfg.baseline     = [-0.2 -0.05];
    %  cfg.baselinetype = 'relative';
    cfg.zlim=z_lim;
    cfg.xlim= x_lim;
    cfg.figure = 'gcf';
    fig1_comps.p1=ft_singleplotTFR(cfg, phases(i));
    xline(0, '--k', 'LineWidth',2);
    c = colorbar;
    c.Label.String = 'Relatív erő';
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
    titles=["Fázisos-Tónusos", "Fázisos-Ébrenlét", "Tónusos-Ébrenlét"];

    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(squeeze(stats(i).time), squeeze(stats(i).freq), squeeze(stats(i).stat));
    colormap(ax(index), "jet")  
    set(gca,'YDir','normal')
    ax=gca;
    ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xline(0, '--k', 'LineWidth',2)
    caxis([-2 2])
    xticks([0 0.2 0.4 0.6])
    c = colorbar;
    c.Label.String = 'T érték';
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

    % Plot the p values in the third column

    % Set up figure standards
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(stats(i).time, squeeze(stats(i).freq), stat(i).sigprobmatrix);
    colormap(ax(index),"bone")
    set(gca,'YDir','normal')
    ax=gca;
    ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
    xticks([0 0.2 0.4 0.6])
    xline(0, '--k', 'LineWidth',2)
    caxis([0.001 0.05])
    c = colorbar;
    c.Label.String = 'p érték';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    fig1_comps.plotTitle=title({strcat(titles(i))},'fontsize', font_size);
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');  
end

%% PLOT FOR THE MEAN OF 20-40HZ, segmented


% Open files WITHOUT baseline correction 
clear all
load('d:\ANT_HEP\HEP_timefreq_512_70hz\tonic_timefreq_without_baseline_7_45');
load('d:\ANT_HEP\HEP_timefreq_512_70hz\phasic_timefreq_without_baseline_7_45');
load('d:\ANT_HEP\HEP_timefreq_512_70hz\wake_timefreq_without_baseline_7_45');
%load('d:\ANT_HEP\HEP_timefreq_512_70hz\nrem_timefreq_without_baseline_7_45');


patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 
%sleep_wake_phases = {'nrem'};
 sleep_wake_phases = {'phasic', 'tonic','wake'};

% Save into a big structure
for row=1:length(sleep_wake_phases)
    for column=1:length(patients)
        if row==1
            patients_timefreq{1,column} =phasic_timefreq{1,column};
        elseif row==2
            patients_timefreq{2,column}= tonic_timefreq{1, column};
        else
            patients_timefreq{3,column}= wake_timefreq{1, column};
        end
    end
end
% % Save into a big structure
% for row=1:length(sleep_wake_phases)
%     for column=1:length(patients)
%         if row==1
%             patients_timefreq{1,column} =nrem_timefreq{1,column};
%         end
%     end
% end
% clear column row

% Average the ANT channels, átlaolni 20-40hz között
for row=1:size(patients_timefreq,1)
    for column=1:size(patients_timefreq,2)
        patients_timefreq{row,column}.powspctrm = mean(patients_timefreq{row,column}.powspctrm,1);
        patients_timefreq{row,column}.label={'General-ANT'};
        patients_timefreq{row,column}.powspctrmnew = mean(patients_timefreq{row,column}.powspctrm(:,14:34,:));
        patients_timefreq{row,column}.powspctrm =  patients_timefreq{row,column}.powspctrmnew ;
        patients_timefreq{row,column}.freq = 30 ;
    end
end
clear column row
pre_atlag_ph=0;
pre_atlag_to=0;
pre_atlag_wa=0;
%14-34

for i=1:size(patients_timefreq,1)
    for j=1:size(patients_timefreq,2)
        if i==1
        atlag_ph_null= patients_timefreq{i,j}.powspctrm;
        pre_atlag_ph=atlag_ph_null+pre_atlag_ph;
        end
         if i==2
         atlag_to_null= patients_timefreq{i,j}.powspctrm;
        pre_atlag_to=atlag_to_null+pre_atlag_to;
         end 
         if i==3
          atlag_wa_null= patients_timefreq{i,j}.powspctrm;
        pre_atlag_wa=atlag_wa_null+pre_atlag_wa;
         end
    end
end

atlag_ph=squeeze(pre_atlag_ph/11);
atlag_to=squeeze(pre_atlag_to/11);
atlag_wa=squeeze(pre_atlag_wa/11);

plot(atlag_ph, 'b')
hold on
plot(atlag_to, 'r')
hold on
plot(atlag_wa, 'y')
legend ('phasic', 'tonic', 'wake')
% 
% % Segment into 75 ms segments
% baseline_indices=1:8; % baseline: -125 ms - -50ms
% start=19; % first segment starts at 50ms
% % Creating the segments, 
% for i=1:floor((length(patients_timefreq{1,1}.time)-start)/length(baseline_indices))
%     timepoints=patients_timefreq{1,1}.time;
%     segments{i}=timepoints(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
%     segment_indices{i}=(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
% end
% baseline_window=timepoints(baseline_indices);


% Segment into 150 ms segments
baseline_indices=1:16; % baseline: -125 ms - -50ms x2

start=19; % first segment starts at 50ms
% Creating the segments, 
for i=1:floor((length(patients_timefreq{1,1}.time)-start)/length(baseline_indices))
    timepoints=patients_timefreq{1,1}.time;
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
            act_segment{j,k}.label=patients_timefreq{j,k}.label;
            act_segment{j,k}.dimord=patients_timefreq{j,k}.dimord;
            act_segment{j,k}.freq=patients_timefreq{j,k}.freq;
            act_segment{j,k}.cfg=patients_timefreq{j,k}.cfg;
            act_segment{j,k}.time=patients_timefreq{j,k}.time(1,segment_indices{i});
            act_segment{j,k}.powspctrm=patients_timefreq{j,k}.powspctrm(:,:,segment_indices{i});
        
            bl_segment{j,k}.label=patients_timefreq{j,k}.label;
            bl_segment{j,k}.dimord=patients_timefreq{j,k}.dimord;
            bl_segment{j,k}.freq=patients_timefreq{j,k}.freq;
            bl_segment{j,k}.cfg=patients_timefreq{j,k}.cfg;
            bl_segment{j,k}.time=patients_timefreq{j,k}.time(1,segment_indices{i}); %using this time interval, because the statistical analysis runs only with the same time points
            c=patients_timefreq{j,k}.powspctrm(:,:,baseline_indices);
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
    cfg.avgoverfreq='no';
    cfg.avgovertime='no';
    cfg.numrandomization = 'all'; 
    cfg.correctm = 'cluster'; 
    cfg.clusterstatistic = 'maxsum'; 
    cfg.clusteralpha = 0.05; 
    cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
    cfg.alpha = 0.05; 
    cfg.uvar = 1; 
    cfg.ivar = 2; 

    % Run the analysis
    % stat_n_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
     stat_p_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
     stat_t_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{2,:}, bl_segment{2,:}); 
     stat_w_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{3,:}, bl_segment{3,:}); 
end

%Concatenate the segments 
 stats=[stat_p_timefreq_seg; stat_t_timefreq_seg; stat_w_timefreq_seg];
%stats=[stat_n_timefreq_seg];


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


plot(atlag_ph, 'b')
hold on
plot(atlag_to, 'r')
hold on
plot(atlag_wa, 'y')


grandaverage_phasic.avg=atlag_ph.';
grandaverage_phasic.time=patients_timefreq{1,1}.time;
grandaverage_phasic.label='ecg';
grandaverage_tonic.avg=atlag_to.';
grandaverage_tonic.time=patients_timefreq{1,1}.time;
grandaverage_tonic.label='ecg';
grandaverage_wake.avg=atlag_wa.';
grandaverage_wake.time=patients_timefreq{1,1}.time;
grandaverage_wake.label='ecg';


%baselinekorrigált
baseline_ph=mean(grandaverage_phasic.avg(1:8));
baseline_to=mean(grandaverage_tonic.avg(1:8));
baseline_wa=mean(grandaverage_wake.avg(1:8));

j=1;
for i=9:86
    grandaverage_phasic.corr(j)=grandaverage_phasic.avg(1,i)/baseline_ph;
    grandaverage_tonic.corr(j)=grandaverage_tonic.avg(1,i)/baseline_to;
    grandaverage_wake.corr(j)=grandaverage_wake.avg(1,i)/baseline_wa;
    j=j+1;
end
% Plot the results

%plotting
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;
 fig1_comps.p1=plot(grandaverage_phasic.time(19:79), grandaverage_phasic.corr(11:71), 'DisplayName', char(grandaverage_phasic.label) + ", phasic");
hold on
 fig1_comps.p2=plot(grandaverage_tonic.time(19:79), grandaverage_tonic.corr(11:71), 'DisplayName', char(grandaverage_tonic.label) + ", tonic");
hold on
 fig1_comps.p3=plot(grandaverage_wake.time(19:79), grandaverage_wake.corr(11:71), 'DisplayName', char(grandaverage_wake.label) + ", wake");
 hold on
% 
%  
% for i=1:size(unified{i,1}.time,2)
%     if unified{i,1}.prob(i) < 0.05
%        [row, col]=find(grandaverage_phasic.time==stat_p.time(i))       ;
%        hold on
%        fig1_comps.p5=plot(grandaverage_phasic.time(col), -0.5, 'k.', 'handlevisibility', 'off');
%        set(fig1_comps.p5, 'Color', PS.MyBlack, 'MarkerSize',10);
%     end
% end

areas = bwboundaries(squeeze(unified{1,1}.prob)<0.05);
for i=1:size(areas,1)
    bw_phasic(i,1)=min(areas{i,1}(:,1));
    bw_phasic(i,2)=max(areas{i,1}(:,1));
end


areas = bwboundaries(squeeze(unified{2,1}.prob)<0.05);
for i=1:size(areas,1)
    bw_tonic(i,1)=min(areas{i,1}(:,1));
    bw_tonic(i,2)=max(areas{i,1}(:,1));
end


areas = bwboundaries(squeeze(unified{3,1}.prob)<0.05);
for i=1:size(areas,1)
    bw_wake(i,1)=min(areas{i,1}(:,1));
    bw_wake(i,2)=max(areas{i,1}(:,1));
end


shaded = [grandaverage_phasic.time(1,bw_phasic(1,1)) grandaverage_phasic.time(1,bw_phasic(1,1)) grandaverage_phasic.time(1,bw_phasic(1,2)) grandaverage_phasic.time(1,bw_phasic(1,2))];



 xline(0.054, '--k', 'LineWidth',2, 'handlevisibility', 'off')
fig1_comps.plotTitle=title('Power modulations at post-baseline period - 20-40 Hz averaged');
fig1_comps.plotXLabel = xlabel('Time (sec)');
fig1_comps.plotYLabel =ylabel({'Amplitude (μV)'});


shaded = [grandaverage_phasic.time(1,bw_phasic(1,1)+18) grandaverage_phasic.time(1,bw_phasic(1,1)+18) grandaverage_phasic.time(1,bw_phasic(1,2)+18) grandaverage_phasic.time(1,bw_phasic(1,2)+18)];
fill(shaded, [0.91 0.955 0.955 0.91], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');


shaded = [grandaverage_phasic.time(1,bw_phasic(2,1)+18) grandaverage_phasic.time(1,bw_phasic(2,1)+18) grandaverage_phasic.time(1,bw_phasic(2,2)+18) grandaverage_phasic.time(1,bw_phasic(2,2)+18)];
fill(shaded, [0.94 0.97 0.97 0.94], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');


% shaded = [grandaverage_phasic.time(1,bw_phasic(3,1)+18) grandaverage_phasic.time(1,bw_phasic(3,1)+18) grandaverage_phasic.time(1,bw_phasic(3,2)+18) grandaverage_phasic.time(1,bw_phasic(3,2)+18)];
% fill(shaded, [0.945 0.96 0.96 0.945], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');

fig1_comps.plotLegend =legend('Phasic REM','Tonic REM','Wake', 'significant clusters');
 hold on
legendX0 = .78; legendY0 = .82; legendWidth = .05; legendHeight = .04;
set(fig1_comps.plotLegend, 'position', [legendX0, legendY0, legendWidth, ...
    legendHeight], 'Box', 'on');
set(fig1_comps.plotLegend, 'FontSize', 10, 'LineWidth', 0.5);

set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman','FontSize', 14);
set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 18, 'FontWeight' , 'bold');
set(fig1_comps.p1, 'Color', PS.Blue4, 'LineWidth', 1.3);
set(fig1_comps.p2, 'Color', PS.Red4, 'LineWidth',1.3);
set(fig1_comps.p3, 'Color', PS.Yellow4, 'LineWidth', 1.3);
xlim([0 0.65])
ylim([0.9 1.14])

fill(shaded, [-1 1 1 -1], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
