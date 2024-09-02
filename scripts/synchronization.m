%% ANT - SCALP synchronization
% phase locking value calculation
clear all
ft_defaults

outfolder = 'd:\ANT_HEP\HEP_synchronization_512_70hz'; 
patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 
sleep_wake_phases = {'nrem'};
% sleep_wake_phases = {'phasic', 'tonic', 'wake'};
frequency=7:1:45;
fs=512;
% Choose between the two options whether you'd like to use baseline
% correction or not.
tobaseline = "yes";
%tobaseline = "no";


for i=1:length(sleep_wake_phases)
    for j=1:length(patients)
        % Load the epoch data for the current patient and sleep-wake phase
        cd 'd:\ANT_HEP\HEP_epochdata_512_70hz';
        load(strcat(char(patients(j)), '_clean_', char(sleep_wake_phases(i)), '_hep_epochs'));
        
        % Load segment data based on the sleep-wake phase (only to create 
        % structure fields)
        if strcmpi(sleep_wake_phases(i),'tonic') == 1 
               cd 'f:\ANT_HEP\HEP_segments_512_70hz\tonic';
            elseif  strcmpi(sleep_wake_phases(i),'phasic') == 1
                cd 'f:\ANT_HEP\HEP_segments_512_70hz\phasic';
            elseif  strcmpi(sleep_wake_phases(i),'wake') == 1
                cd 'f:\ANT_HEP\HEP_segments_512_70hz\wake';
        else
            cd 'd:\ANT_HEP\HEP_segments_512_70hz\nrem';
        end
        segment=load(strcat(char(patients(j)), '_clean_', char(sleep_wake_phases(i))));
        fns=fieldnames(segment);

        % Create a structure for the analysis while using the epoch data
        file.label=segment.(fns{1}).label;
        file.fsample=fs;
        for k = 1:size(epochdata,1)
            file.trial{1,k}=squeeze(epochdata(:,:,k));
            file.time{1,k}=linspace(-0.2,0.8,fs);
        end
      
        for d = 0:size(epochdata,3)-1 
            file.sampleinfo(d+1,1)=d*fs+1;
            file.sampleinfo(d+1,2)=(d+1)*fs;
        end
       
        % Remove reference Tp9 & Tp10 and unnecessary channels         
        cfg=[];
        cfg.channel={'all', '-EEG Tp10-G2','-EEG Tp9-G2', '-EEG ECG-G2', '-EEG Zyg2-Zyg1', '-EEG ecg-G2'};
        if j==6 %Piri (patient no6) has an all 0 channel, th11
            cfg=[];
            cfg.channel={'all', '-EEG Tp10-G2','-EEG Tp9-G2', '-EEG Th11-G2', '-EEG ECG-G2', '-EEG Zyg2-Zyg1', '-EEG ecg-G2'};
        end
        if i==1
            nrem_clean=ft_selectdata(cfg,file);
        end   
%         if i==1
%             phasic_clean=ft_selectdata(cfg,file);
%             elseif i==2
%             tonic_clean=ft_selectdata(cfg,file);
%             else
%             wake_clean = ft_selectdata(cfg,file);
%         end      

        % Compute cross spectral density
        cfg              = [];
        cfg.output       = 'powandcsd';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.pad          = 2; %ONLY WORKS THIS WAY
        cfg.foi          = 7:1:45;                       % analysis 10 to 45 Hz in steps of 1 Hz    
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.15;  
        cfg.toi          = -0.125:0.01:0.725;                  % time window "slides" from -0.125 to 0.725 sec in steps of 0.01 sec (10 ms)       
        if i==1
            nrem_csd = ft_freqanalysis(cfg, nrem_clean);
        end        
% if i==1
%             phasic_csd = ft_freqanalysis(cfg, phasic_clean);
%             elseif i==2
%             tonic_csd  =  ft_freqanalysis(cfg, tonic_clean);
%             else
%             wake_csd  =  ft_freqanalysis(cfg, wake_clean);
%         end
        
        clear file
        % Exctract ANT-SCALP pairs
       if i==1
            cfg=[];
            cfg.channelcmb=nrem_csd.labelcmb([find(strncmpi(nrem_csd.labelcmb(:,1),'EEG BTH',7));... 
                    find(strncmpi(nrem_csd.labelcmb(:,1),'EEG JTH',7));...
                    find(strncmpi(nrem_csd.labelcmb(:,1),'EEG TH',6))],:); 
            nrem_ant_scalp=ft_selectdata(cfg,nrem_csd);
       end
%         if i==1
%             cfg=[];
%             cfg.channelcmb=phasic_csd.labelcmb([find(strncmpi(phasic_csd.labelcmb(:,1),'EEG BTH',7));... 
%                     find(strncmpi(phasic_csd.labelcmb(:,1),'EEG JTH',7));...
%                     find(strncmpi(phasic_csd.labelcmb(:,1),'EEG TH',6))],:); 
%             phasic_ant_scalp=ft_selectdata(cfg,phasic_csd);
%             elseif i==2
%             cfg=[];
%             cfg.channelcmb=tonic_csd.labelcmb([find(strncmpi(tonic_csd.labelcmb(:,1),'EEG BTH',7));... 
%                  find(strncmpi(tonic_csd.labelcmb(:,1),'EEG JTH',7));...
%                 find(strncmpi(tonic_csd.labelcmb(:,1),'EEG TH',6))],:); 
%             tonic_ant_scalp=ft_selectdata(cfg,tonic_csd);
%             else
%             cfg=[]; 
%             cfg.channelcmb=wake_csd.labelcmb([find(strncmpi(wake_csd.labelcmb(:,1),'EEG BTH',7));... 
%                  find(strncmpi(wake_csd.labelcmb(:,1),'EEG JTH',7));...
%                 find(strncmpi(wake_csd.labelcmb(:,1),'EEG TH',6))],:); 
%             wake_ant_scalp=ft_selectdata(cfg,wake_csd);
%         end

        % compute plv
        cfg=[];
        cfg.method     = 'plv';
        if i==1
            plv_nrem    = ft_connectivityanalysis(cfg, nrem_ant_scalp);
            plv_nrem.plvspctrm=abs(plv_nrem.plvspctrm);   
         end
%         if i==1
%             plv_phasic    = ft_connectivityanalysis(cfg, phasic_ant_scalp);
%             plv_phasic.plvspctrm=abs(plv_phasic.plvspctrm);   
%             elseif i==2
%             plv_tonic    = ft_connectivityanalysis(cfg, tonic_ant_scalp);
%             plv_tonic.plvspctrm=abs(plv_tonic.plvspctrm);   
%             else
%             plv_wake    = ft_connectivityanalysis(cfg, wake_ant_scalp);
%             plv_wake.plvspctrm=abs(plv_wake.plvspctrm);   
%         end
    
        % Average ANT - fr, cent, parietal pairs
        if i==1
            nrem_all_synch.fr(j,:,:)=squeeze(mean(plv_nrem.plvspctrm([find(strncmpi(plv_nrem.labelcmb(:,2),'EEG F',5))],:,:),1));
            nrem_all_synch.ce(j,:,:)=squeeze(mean(plv_nrem.plvspctrm([find(strncmpi(plv_nrem.labelcmb(:,2),'EEG C',5)); find(strncmpi(plv_nrem.labelcmb(:,2),'EEG T',5))],:,:),1)); 
            nrem_all_synch.pa(j,:,:)=squeeze(mean(plv_nrem.plvspctrm([find(strncmpi(plv_nrem.labelcmb(:,2),'EEG P',5)); find(strncmpi(plv_nrem.labelcmb(:,2),'EEG O',5))],:,:),1)); 
            nrem_all_synch.time=plv_nrem.time;
            nrem_all_synch.freq=plv_nrem.freq;
            nrem_all_synch.dimord=plv_nrem.dimord;
            nrem_all_synch.label{1,1}='ant-scalp synch';


%             phasic_all_synch.fr(j,:,:)=squeeze(mean(plv_phasic.plvspctrm([find(strncmpi(plv_phasic.labelcmb(:,2),'EEG F',5))],:,:),1));
%             phasic_all_synch.ce(j,:,:)=squeeze(mean(plv_phasic.plvspctrm([find(strncmpi(plv_phasic.labelcmb(:,2),'EEG C',5)); find(strncmpi(plv_phasic.labelcmb(:,2),'EEG T',5))],:,:),1)); 
%             phasic_all_synch.pa(j,:,:)=squeeze(mean(plv_phasic.plvspctrm([find(strncmpi(plv_phasic.labelcmb(:,2),'EEG P',5)); find(strncmpi(plv_phasic.labelcmb(:,2),'EEG O',5))],:,:),1)); 
%             phasic_all_synch.time=plv_phasic.time;
%             phasic_all_synch.freq=plv_phasic.freq;
%             phasic_all_synch.dimord=plv_phasic.dimord;
%             phasic_all_synch.label{1,1}='ant-scalp synch';



            elseif i==2
            tonic_all_synch.fr(j,:,:)=squeeze(mean(plv_tonic.plvspctrm([find(strncmpi(plv_tonic.labelcmb(:,2),'EEG F',5))],:,:),1));
            tonic_all_synch.ce(j,:,:)=squeeze(mean(plv_tonic.plvspctrm([find(strncmpi(plv_tonic.labelcmb(:,2),'EEG C',5)); find(strncmpi(plv_tonic.labelcmb(:,2),'EEG T',5))],:,:),1)); 
            tonic_all_synch.pa(j,:,:)=squeeze(mean(plv_tonic.plvspctrm([find(strncmpi(plv_tonic.labelcmb(:,2),'EEG P',5)); find(strncmpi(plv_tonic.labelcmb(:,2),'EEG O',5))],:,:),1)); 
            tonic_all_synch.time=plv_tonic.time;
            tonic_all_synch.freq=plv_tonic.freq;
            tonic_all_synch.dimord=plv_tonic.dimord;
            tonic_all_synch.label{1,1}='ant-scalp synch';


            else
            wake_all_synch.fr(j,:,:)=squeeze(mean(plv_wake.plvspctrm([find(strncmpi(plv_wake.labelcmb(:,2),'EEG F',5))],:,:),1));
            wake_all_synch.ce(j,:,:)=squeeze(mean(plv_wake.plvspctrm([find(strncmpi(plv_wake.labelcmb(:,2),'EEG C',5)); find(strncmpi(plv_wake.labelcmb(:,2),'EEG T',5))],:,:),1)); 
            wake_all_synch.pa(j,:,:)=squeeze(mean(plv_wake.plvspctrm([find(strncmpi(plv_wake.labelcmb(:,2),'EEG P',5)); find(strncmpi(plv_wake.labelcmb(:,2),'EEG O',5))],:,:),1)); 
            wake_all_synch.time=plv_wake.time;
            wake_all_synch.freq=plv_wake.freq;
            wake_all_synch.dimord=plv_wake.dimord;
            wake_all_synch.label{1,1}='ant-scalp synch';

        end   
    end
end   

cd (outfolder)
save ('ant_scalp_plv_without_baseline_nrem_7_45.mat', 'nrem_all_synch')
 

% Baseline correction

if tobaseline == "yes"
%     phases=[phasic_all_synch tonic_all_synch wake_all_synch];
    phases=[nrem_all_synch];
    baselined_nrem_all_synch=nrem_all_synch;
%     baselined_phasic_all_synch=phasic_all_synch;
%     baselined_tonic_all_synch=tonic_all_synch;
%     baselined_wake_all_synch=wake_all_synch;
    % Perform baseline correction on the data
    for i =1:length(phases)
       fn=fieldnames(phases(i));
        for j=1:3
            cfg=[];
            cfg.baseline     = [-0.2 -0.05];
            cfg.baselinetype = 'relative';
            cfg.parameter = fn(j);
            phases(i).dimord='chan_freq_time';
            [freq] = ft_freqbaseline(cfg, phases(i));
            if i==1
                 baselined_nrem_all_synch.(fn{j})=freq.(fn{j});
%                 baselined_phasic_all_synch.(fn{j})=freq.(fn{j});
%             elseif i==2
%                 baselined_tonic_all_synch.(fn{j})=freq.(fn{j});
%             else
%                 baselined_wake_all_synch.(fn{j})=freq.(fn{j});
            end
            clear freq
        end
    end
 
    %save the created variables
    save ('ant_scalp_plv_with_baseline_nrem_7_45.mat', 'baselined_nrem_all_synch')
    %this way the code will flow continously with the baselinecorrected data
     nrem_all_synch=baselined_nrem_all_synch;

    %     phasic_all_synch=baselined_phasic_all_synch;
%     tonic_all_synch=baselined_tonic_all_synch;
%     wake_all_synch=baselined_wake_all_synch;
end


% Average over patients
name=fieldnames(nrem_all_synch);
clear k
for k=1:3
      nrem_plv_av(k,:,:)=squeeze(mean(nrem_all_synch.(name{k}),1));
%     phasic_plv_av(k,:,:)=squeeze(mean(phasic_all_synch.(name{k}),1));
%     tonic_plv_av(k,:,:)=squeeze(mean(tonic_all_synch.(name{k}),1));
%     wake_plv_av(k,:,:)=squeeze(mean(wake_all_synch.(name{k}),1));
end
%now each row represents the grand average for the ANT - fr, cent, parietal pairs


% Plot  
% Some initial configurations
% phases=[phasic_plv_av; tonic_plv_av; wake_plv_av];
phases=[nrem_plv_av];
freq=nrem_all_synch.freq;
time=nrem_all_synch.time;
font_size=10;
x_label = 'Time (sec)';
y_label = 'Frequency (Hz)';
% phase_titles=["Phasic", "Tonic", "Wake"];
phase_titles=["NREM"];

place_titles=["Frontal EEG & ANT synchronization", "Central & ANT synchronization", "Parietal & ANT synchronization"];
clear i j
figure
for i=1:length(phases) %for phase
    for j=1:3 %for eeg area
        subplot (3,1, (i-1)*3+j)
        fig1_comps.fig = gcf;
        fig1_comps.p1=imagesc(squeeze(time), squeeze(freq), squeeze(phases((i-1)*3+j,:,:)));
        set(gca,'YDir','normal')
        xline(0, '--w', 'LineWidth',2)
        colormap="jet";
        c = colorbar;
        c.Label.String = 'Phase Locking Value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline=="yes"
           fig1_comps.plotTitle=title(strcat(phase_titles(i),' - ', place_titles(j), " - baseline corrected"),'fontsize', font_size);
           caxis([0.8 1.6])
        else
            fig1_comps.plotTitle=title(strcat(phase_titles(i),' - ', place_titles(j) ),'fontsize', font_size);
            caxis([0.13 0.28])
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    end
end


% Average over all scalp region-ANT pairs
% phasic_plv_av_all=mean(phasic_plv_av,1);
% tonic_plv_av_all=mean(tonic_plv_av,1);
% wake_plv_av_all=mean(wake_plv_av,1);
nrem_plv_av_all=mean(nrem_plv_av,1);
phases_avg=[nrem_plv_av_all];
% phases_avg=[phasic_plv_av_all; tonic_plv_av_all; wake_plv_av_all];


% Plot
clear i 
figure
for i=1:1 %for phase
    %subplot (3,1, i)
    fig1_comps.fig = gcf;
    fig1_comps.p1=imagesc(squeeze(time), squeeze(freq), squeeze(phases_avg(i,:,:)));
    set(gca,'YDir','normal')
    xline(0, '--w', 'LineWidth',2)
    colormap="jet";
    c = colorbar;
    c.Label.String = 'Phase Locking Value';
    c.Label.FontName= 'Times New Roman';
    c.Label.FontSize=10;
    if tobaseline=="yes"
        fig1_comps.plotTitle=title(strcat(phase_titles(i)," - Averaged EEG & ANT synchronization - baseline corrected"),'fontsize', font_size);
        caxis([0.8 1.4])
    else
        fig1_comps.plotTitle=title(strcat(phase_titles(i)," - Averaged EEG & ANT synchronization"),'fontsize', font_size);
        caxis([0.13 0.26])
    end
    fig1_comps.plotXLabel = xlabel(x_label) ;
    fig1_comps.plotYLabel =ylabel(y_label);
    set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
    set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
end



%% STATISTICS

clear all

% Choose one of the following
tobaseline = "yes"; %for baseline corrected data
tobaseline = "no"; %for non baseline corrected data

if tobaseline=="yes"
    load ('ant_scalp_plv_with_baseline_nrem_7_45.mat')
    %this way the code will flow continously with the baselinecorrected data
    nrem_all_synch=baselined_nrem_all_synch;
%     phasic_all_synch=baselined_phasic_all_synch;
%     tonic_all_synch=baselined_tonic_all_synch;
%     wake_all_synch=baselined_wake_all_synch;
else
    load 'F:\ANT_HEP\HEP_synchronization_512_70hz\ant_scalp_plv_without_baseline_phasic_tonic_wake_7_45.mat';
end

patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'}; 
sleep_wake_phases = {'nrem'};

% sleep_wake_phases = {'phasic', 'tonic', 'wake'};


% Create a structure so that the statistical analysis can run (divide the
% structure into a cell array that contains as much 1x1 structures as many
% patients we have)
nrem_all_synch.dimord='chan_freq_time';
phasic_all_synch.dimord='chan_freq_time';
tonic_all_synch.dimord='chan_freq_time';
wake_all_synch.dimord='chan_freq_time';

f_names=fieldnames(nrem_all_synch);
for i=1:length(patients)
    for j=4:length(f_names)
        nrem{1,i}.(f_names{j,1})=nrem_all_synch.(f_names{j,1});
%         phasic{1,i}.(f_names{j,1})=phasic_all_synch.(f_names{j,1});
%         tonic{1,i}.(f_names{j,1})=tonic_all_synch.(f_names{j,1});
%         wake{1,i}.(f_names{j,1})=wake_all_synch.(f_names{j,1});
    end
   for k=1:3
       % fr, ce, pa - regions of the scalp
       nrem{1,i}.(f_names{k,1})=nrem_all_synch.(f_names{k,1})(i,:,:);
%         phasic{1,i}.(f_names{k,1})=phasic_all_synch.(f_names{k,1})(i,:,:);
%         tonic{1,i}.(f_names{k,1})=tonic_all_synch.(f_names{k,1})(i,:,:);
%         wake{1,i}.(f_names{k,1})=wake_all_synch.(f_names{k,1})(i,:,:);
   end
    %create a field with the averaged eeg values (over the fr ce pa regions)
%     phasic{1,i}.scalpavg=mean([phasic_all_synch.(f_names{1,1})(i,:,:); phasic_all_synch.(f_names{2,1})(i,:,:); phasic_all_synch.(f_names{3,1})(i,:,:)],1);
%     tonic{1,i}.scalpavg=mean([tonic_all_synch.(f_names{1,1})(i,:,:); tonic_all_synch.(f_names{2,1})(i,:,:); tonic_all_synch.(f_names{3,1})(i,:,:)],1);
%     wake{1,i}.scalpavg=mean([wake_all_synch.(f_names{1,1})(i,:,:); wake_all_synch.(f_names{2,1})(i,:,:); wake_all_synch.(f_names{3,1})(i,:,:)],1);
    nrem{1,i}.scalpavg=mean([nrem_all_synch.(f_names{1,1})(i,:,:); nrem_all_synch.(f_names{2,1})(i,:,:); nrem_all_synch.(f_names{3,1})(i,:,:)],1);

end
%quickly save the newly created structure
if tobaseline=="yes"
    cd 'd:\ANT_HEP\HEP_synchronization_512_70hz'
    save("patientdata_nrem_synch_with_baseline_7_45", 'nrem', '-v7.3');
    save("patientdata_phasic_synch_with_baseline_7_45", 'phasic', '-v7.3');
    save("patientdata_tonic_synch_with_baseline_7_45", 'tonic', '-v7.3');
    save("patientdata_wake_synch_with_baseline_7_45", 'wake', '-v7.3');
else
    cd 'f:\ANT_HEP\HEP_synchronization_512_70hz'
    save("patientdata_phasic_synch_without_baseline_7_45", 'phasic', '-v7.3');
    save("patientdata_tonic_synch_without_baseline_7_45", 'tonic', '-v7.3');
    save("patientdata_wake_synch_without_baseline_7_45", 'wake', '-v7.3');
end
load("patientdata_phasic_synch_with_baseline_7_45.mat")
load("patientdata_tonic_synch_with_baseline_7_45.mat")
load("patientdata_wake_synch_with_baseline_7_45.mat")
load("patientdata_nrem_synch_with_baseline_7_45.mat")

% Create the design matrix
subj = length(patients);
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;


parameters=["fr" "ce" "pa" "scalpavg"];

for i=1:length(parameters)
    % Set analysis parameters
    cfg = []; 
    cfg.design = design;
    cfg.parameter= convertStringsToChars(parameters(i)); 
    cfg.method = 'montecarlo'; 
    cfg.statistic = 'ft_statfun_depsamplesT';  
    cfg.avgoverfreq='no';
    cfg.avgovertime='no';
    cfg.numrandomization = 'all'; 
    cfg.alpha = 0.05; 
    cfg.latency=[0.05 0.65];
    cfg.tail = 0;  
    cfg.uvar = 1; 
    cfg.ivar = 2; 
    
    
    % Perform cluster-based statistics
    [stat_p_n_synch] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
    [stat_t_n_synch] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
    [stat_w_n_synch] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
%     [stat_p_t_synch] = ft_freqstatistics(cfg, phasic{1,:}, tonic{1,:}); 
%     [stat_p_w_synch] = ft_freqstatistics(cfg,  phasic{1,:}, wake{1,:}); 
%     [stat_t_w_synch] = ft_freqstatistics(cfg,  tonic{1,:}, wake{1,:}); 
    
    % Save the results
    cd 'd:\ANT_HEP\HEP_synchronization_512_70hz'
    if tobaseline == "yes"
        save(strcat("stat_p_n_snych_with_baseline_7_45_", parameters(i)), 'stat_p_n_synch', '-v7.3');
        save(strcat("stat_t_n_synch_with_baseline_7_45_", parameters(i)), 'stat_t_n_synch', '-v7.3');
        save(strcat("stat_w_n_synch_with_baseline_7_45_", parameters(i)), 'stat_w_n_synch', '-v7.3');
    else
        save(strcat("stat_p_t_snych_without_baseline_7_45_", parameters(i)), 'stat_p_t_synch', '-v7.3');
        save(strcat("stat_p_w_synch_without_baseline_7_45_", parameters(i)), 'stat_p_w_synch', '-v7.3');
        save(strcat("stat_t_w_synch_without_baseline_7_45_", parameters(i)), 'stat_t_w_synch', '-v7.3');
    end

    labels={{'ANT & Frontal EEG synchronization'},{'ANT & Central EEG synchronization'},{'ANT & Parietal EEG synchronization'}, {'ANT & Averaged EEG synchronization'}};

    %%%% Plot
    
    % Initialize some settings
    font_size = 10;
    x_label = 'Time (sec)';
    y_label = 'Frequency (Hz)';
    stat=[stat_p_n_synch stat_t_n_synch stat_w_n_synch];
    titles=[" in Phasic-NREM", " in Tonic-NREM", " in Wake-NREM"];
    
    % Plot prob, p values
    PS = PLOT_STANDARDS();
    figure;
    for j=1:3
        fig1_comps.fig = gcf;
        subplot(3, 1, j)
        fig1_comps.p1=imagesc(squeeze(stat(j).time), squeeze(stat(j).freq), squeeze(stat(j).prob));
        set(gca,'YDir','normal')
        xline(0, '--k', 'LineWidth', 2)
        colormap(bone)
        caxis([0.001 0.05])
        c = colorbar;
        c.Label.String = 'p value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=font_size;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j), " - baseline corrected"), 'fontsize', font_size);
        else
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j)), 'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
        clear fig1_comps
    end
    
    
    % Plot stat, T values
    PS = PLOT_STANDARDS();
    figure;
    for k=1:3
        fig1_comps.fig = gcf;
        subplot(3,1,k)
        fig1_comps.p1=imagesc(squeeze(stat(k).time), squeeze(stat(k).freq), squeeze(stat(k).stat));
        set(gca,'YDir','normal')
        xline(0, '--w', 'LineWidth',2)
        caxis([-1 2])
        c = colorbar;
        c.Label.String = 'T value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j), " - baseline corrected"),'fontsize', font_size);
        else
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j)),'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    end

end


%% Cluster Analysis


for i=1:length(parameters)
    % Set analysis parameters
    cfg = []; 
    cfg.design = design;
    cfg.parameter= convertStringsToChars(parameters(i)); 
    cfg.method = 'montecarlo'; 
    cfg.statistic = 'ft_statfun_depsamplesT';  
    cfg.avgoverfreq='no';
    cfg.avgovertime='no';
    cfg.numrandomization = 'all'; 
    cfg.correctm = 'cluster'; 
    cfg.clusterstatistic = 'maxsum'; 
    cfg.clusteralpha = 0.05; 
    cfg.latency=[0.05 0.65];
    cfg.tail = 0;  
    cfg.uvar = 1; 
    cfg.ivar = 2; 
    
    
    % Perform cluster-based statistics
     [stat_p_n_synch_cluster] = ft_freqstatistics(cfg, phasic{1,:}, nrem{1,:}); 
    [stat_t_n_synch_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, nrem{1,:}); 
    [stat_w_n_synch_cluster] = ft_freqstatistics(cfg,  wake{1,:}, nrem{1,:}); 
%     [stat_p_t_synch_cluster] = ft_freqstatistics(cfg, phasic{1,:}, tonic{1,:}); 
%     [stat_p_w_synch_cluster] = ft_freqstatistics(cfg,  phasic{1,:}, wake{1,:}); 
%     [stat_t_w_synch_cluster] = ft_freqstatistics(cfg,  tonic{1,:}, wake{1,:}); 
    
    % Save the results
    cd 'd:\ANT_HEP\HEP_synchronization_512_70hz'
    if tobaseline == "yes"
        save(strcat("stat_p_n_snych_with_baseline_7_45_cluster_", parameters(i)), 'stat_p_n_synch_cluster', '-v7.3');
        save(strcat("stat_t_n_synch_with_baseline_7_45_cluster_", parameters(i)), 'stat_t_n_synch_cluster', '-v7.3');
        save(strcat("stat_w_n_synch_with_baseline_7_45_cluster_", parameters(i)), 'stat_w_n_synch_cluster', '-v7.3');
%         save(strcat("stat_p_t_snych_with_baseline_7_45_cluster_", parameters(i)), 'stat_p_t_synch_cluster', '-v7.3');
%         save(strcat("stat_p_w_synch_with_baseline_7_45_cluster_", parameters(i)), 'stat_p_w_synch_cluster', '-v7.3');
%         save(strcat("stat_t_w_synch_with_baseline_7_45_cluster_", parameters(i)), 'stat_t_w_synch_cluster', '-v7.3');
    else
         save(strcat("stat_p_t_snych_without_baseline_7_45_cluster_", parameters(i)), 'stat_p_t_synch_cluster', '-v7.3');
        save(strcat("stat_p_w_synch_without_baseline_7_45_cluster_", parameters(i)), 'stat_p_w_synch_cluster', '-v7.3');
        save(strcat("stat_t_w_synch_without_baseline_7_45_cluster_", parameters(i)), 'stat_t_w_synch_cluster', '-v7.3');
    end


    %%%% Plot
    
    % Initialize some settings
    font_size = 10;
    x_label = 'Time (sec)';
    y_label = 'Frequency (Hz)';
    stat=[stat_p_n_synch_cluster stat_t_n_synch_cluster stat_w_n_synch_cluster];
    titles=[" in Phasic-NREM", " in Tonic-NREM", " in Wake-NREM"];
    labels={{'ANT & Frontal EEG synchronization'},{'ANT & Central EEG synchronization'},{'ANT & Parietal EEG synchronization'}, {'ANT & Averaged EEG synchronization'}};

    
    % Plot prob, p values
    PS = PLOT_STANDARDS();
    figure;
    for j=1:3
        fig1_comps.fig = gcf;
        subplot(3, 1, j)
        fig1_comps.p1=imagesc(squeeze(stat(j).time), squeeze(stat(j).freq), squeeze(stat(j).prob));
        set(gca,'YDir','normal')
        xline(0, '--k', 'LineWidth', 2)
        colormap(bone)
        caxis([0.001 0.05])
        c = colorbar;
        c.Label.String = 'p value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=font_size;
        if tobaseline== "yes"
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j), " - baseline corrected"), 'fontsize', font_size);
        else
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j)), 'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
        clear fig1_comps
    end
    
    
    % Plot stat, T values
    PS = PLOT_STANDARDS();
    figure;
    for k=1:3
        fig1_comps.fig = gcf;
        subplot(3,1,k)
        fig1_comps.p1=imagesc(squeeze(stat(k).time), squeeze(stat(k).freq), squeeze(stat(k).stat));
        set(gca,'YDir','normal')
        xline(0, '--w', 'LineWidth',2)
        caxis([-1 2])
        c = colorbar;
        c.Label.String = 'T value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j), " - baseline corrected"),'fontsize', font_size);
        else
            fig1_comps.plotTitle=title(strcat(labels{i}, titles(j)),'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    end

end


%% Plot that contains every result from previous sections 

clear all

% Choose one of the followings
tobaseline = "yes"; %for baseline corrected data
tobaseline = "no"; %for non baseline corrected data

cd 'd:\ANT_HEP\HEP_synchronization_512_70hz'
parameters=["fr" "ce" "pa" "scalpavg"];
area_titles=["Frontal EEG & ANT synch - ", "Central EEG & ANT synch - ", "Parietal EEG & ANT synch - ", "Averaged EEG & ANT synch - "];
phase_titles=["NREM"];
% phase_titles=["Phasic", "Tonic", "Wake"];
stat_phase_titles=["Phasic-NREM", "Tonic-NREM", "Wake-NREM"]; 


% Plot

%loop through the different scalp areas
for i=1:length(parameters) 
figure
    % open the data to be grand averaged
    if tobaseline == "yes"
        load("ant_scalp_plv_with_baseline_phasic_tonic_wake_7_45.mat");
        %this way the code will flow continously with the baselinecorrected data
        phasic_all_synch=baselined_phasic_all_synch;
        tonic_all_synch=baselined_tonic_all_synch;
        wake_all_synch=baselined_wake_all_synch;
    else
        load("ant_scalp_plv_without_baseline_phasic_tonic_wake_7_45.mat");
    end
    % average    
    name=fieldnames(phasic_all_synch);
    for k=1:3
        phasic_all_synch.(name{k})=mean(phasic_all_synch.(name{k}),1);
        tonic_all_synch.(name{k})=mean(tonic_all_synch.(name{k}),1);
        wake_all_synch.(name{k})=mean(wake_all_synch.(name{k}),1);
    end
    clear k

    % Average over all scalp region-ANT pairs
    phasic_all_synch.scalpavg=mean([phasic_all_synch.fr; phasic_all_synch.ce; phasic_all_synch.pa],1);
    tonic_all_synch.scalpavg=mean([tonic_all_synch.fr; tonic_all_synch.ce; tonic_all_synch.pa],1);
    wake_all_synch.scalpavg=mean([wake_all_synch.fr; wake_all_synch.ce; wake_all_synch.pa],1);
    
    %create a powspctrm field because plotting functions need that field
    phasic_all_synch.powspctrm=phasic_all_synch.(parameters(i));
    tonic_all_synch.powspctrm=tonic_all_synch.(parameters(i));
    wake_all_synch.powspctrm=wake_all_synch.(parameters(i));

    % open the statistics
    if tobaseline == "yes"
        load(strcat("stat_p_t_snych_with_baseline_7_45_", parameters(i)));
        load(strcat("stat_p_w_synch_with_baseline_7_45_", parameters(i)));
        load(strcat("stat_t_w_synch_with_baseline_7_45_", parameters(i)));
    else
        load(strcat("stat_p_t_snych_without_baseline_7_45_", parameters(i)));
        load(strcat("stat_p_w_synch_without_baseline_7_45_", parameters(i)));
        load(strcat("stat_t_w_synch_without_baseline_7_45_", parameters(i)));
    end
    
    % open the cluster based statistics
    if tobaseline == "yes"
        load(strcat("stat_p_t_snych_with_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_p_w_synch_with_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_t_w_synch_with_baseline_7_45_cluster_", parameters(i)));
    else
        load(strcat("stat_p_t_snych_without_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_p_w_synch_without_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_t_w_synch_without_baseline_7_45_cluster_", parameters(i)));
    end
           
    tcl = tiledlayout(3,4);  
    ax = gobjects(1,12);
    index=1;
    %loop through the phases
    for j = 1:3
        ax(index) = nexttile;
        
        % Plot the grand averaged data in the first column     
        % Some initial configurations
        phases=[phasic_all_synch; tonic_all_synch; wake_all_synch];
        freq=phasic_all_synch.freq;
        time=phasic_all_synch.time;
        font_size=10;
        x_label = 'Time (sec)';
        y_label = 'Frequency (Hz)';
         if tobaseline =="yes"
            z_lim=[0.8 1.6];
        else
           z_lim=[0.13 0.28];
        end
        x_lim=[-0.2 0.6];
        % Set up figure standards
        PS = PLOT_STANDARDS();
        fig1_comps.fig = gcf;
        colormap(ax(index), "jet") 
        cfg = [];
        % Pay attention to cfg.baseline, this is optional to use
        % cfg.baseline     = [-0.2 -0.05];
        %  cfg.baselinetype = 'relative';
        cfg.zlim=z_lim;
        cfg.xlim= x_lim;
        cfg.figure = 'gcf';
        fig1_comps.p1=ft_singleplotTFR(cfg, phases(j));  

        xline(0, '--w', 'LineWidth',2);
        c = colorbar;
        c.Label.String = 'PLV';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=font_size;
        if tobaseline =="yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i), phase_titles(j)), " Grand averaged - baseline corrected"}, 'fontsize', font_size);
        else
           fig1_comps.plotTitle=title({strcat(area_titles(i), phase_titles(j)), " Grand averaged"}, 'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
        
        % Step onto the next tile
        index=index+1;
        ax(index) = nexttile;
    
        % Plot the T values in the second column
        % Initialize some settings
        stat=[stat_p_t_synch stat_p_w_synch stat_t_w_synch];

        % Set up figure standards
        fig1_comps.fig = gcf;
        fig1_comps.p1=imagesc(squeeze(stat(j).time), squeeze(stat(j).freq), squeeze(stat(j).stat));
        xline(0, '--w', 'LineWidth',2);
        colormap(ax(index), "jet")  
        set(gca,'YDir','normal')
         ax=gca;
        ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
        xline(0, '--k', 'LineWidth',2)
        caxis([-2 2])
        c = colorbar;
        c.Label.String = 'T value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "T values - baseline corrected"},'fontsize', font_size);
            caxis([-3 4])
        else
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "T values"},'fontsize', font_size);
            caxis([-2 2])
        end      
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
        fig1_comps.p1=imagesc(stat(j).time, squeeze(stat(j).freq), squeeze(stat(j).prob));
        colormap(ax(index),"bone")
        set(gca,'YDir','normal')
         ax=gca;
        ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
        xline(0, '--k', 'LineWidth', 2)
        caxis([0.001 0.05])
        c = colorbar;
        c.Label.String = 'p value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "p values - baseline corrected"},'fontsize', font_size);
        else
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "p values"},'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');


        % Step onto the next tile
        index=index+1;
        ax(index) = nexttile;
    
        % Plot the significant clusters in the fourth column
    
        %Initialize some settings    
        clusterstats=[stat_p_t_synch_cluster stat_p_w_synch_cluster stat_t_w_synch_cluster];
        % Set up figure standards
        fig1_comps.fig = gcf;
        fig1_comps.p1=imagesc(clusterstats(j).time, squeeze(clusterstats(j).freq), squeeze(clusterstats(j).prob));
        hold on
        colormap(ax(index),"bone")
        set(gca,'YDir','normal')
         ax=gca;
        ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
        xline(0, '--k', 'LineWidth', 2)
        caxis([0.001 0.05])
        c = colorbar;
        c.Label.String = "cluster's p value";
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=font_size;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "clustered p values - baseline corrected"}, 'fontsize', font_size);
        else
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "clustered p values"}, 'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
        clear fig1_comps
       
    end
end



%% Second version with 3x3 layout and clustered p-s overlapping normal p-s

clear all

% Choose one of the followings
%tobaseline = "yes"; %for baseline corrected data
tobaseline = "no"; %for non baseline corrected data

cd 'f:\ANT_HEP\HEP_synchronization_512_70hz'
parameters=["fr" "ce" "pa" "scalpavg"];
area_titles=["Frontal EEG & ANT synch - ", "Central EEG & ANT synch - ", "Parietal EEG & ANT synch - ", "Averaged EEG & ANT synch - "];
phase_titles=["Phasic", "Tonic", "Wake"];
stat_phase_titles=["Phasic-Tonic", "Phasic-Wake", "Tonic-Wake"]; 


% Plot

%loop through the different scalp areas
for i=1:length(parameters) 
figure
    % open the data to be grand averaged
    if tobaseline == "yes"
        load("ant_scalp_plv_with_baseline_phasic_tonic_wake_7_45.mat");
        %this way the code will flow continously with the baselinecorrected data
        phasic_all_synch=baselined_phasic_all_synch;
        tonic_all_synch=baselined_tonic_all_synch;
        wake_all_synch=baselined_wake_all_synch;
    else
        load("ant_scalp_plv_without_baseline_phasic_tonic_wake_7_45.mat");
    end
    % average    
    name=fieldnames(phasic_all_synch);
    for k=1:3
        phasic_all_synch.(name{k})=mean(phasic_all_synch.(name{k}),1);
        tonic_all_synch.(name{k})=mean(tonic_all_synch.(name{k}),1);
        wake_all_synch.(name{k})=mean(wake_all_synch.(name{k}),1);
    end
    clear k

    % Average over all scalp region-ANT pairs
    phasic_all_synch.scalpavg=mean([phasic_all_synch.fr; phasic_all_synch.ce; phasic_all_synch.pa],1);
    tonic_all_synch.scalpavg=mean([tonic_all_synch.fr; tonic_all_synch.ce; tonic_all_synch.pa],1);
    wake_all_synch.scalpavg=mean([wake_all_synch.fr; wake_all_synch.ce; wake_all_synch.pa],1);
    
    %create a powspctrm field because plotting functions need that field
    phasic_all_synch.powspctrm=phasic_all_synch.(parameters(i));
    tonic_all_synch.powspctrm=tonic_all_synch.(parameters(i));
    wake_all_synch.powspctrm=wake_all_synch.(parameters(i));

    % open the statistics
    if tobaseline == "yes"
        load(strcat("stat_p_t_snych_with_baseline_7_45_", parameters(i)));
        load(strcat("stat_p_w_synch_with_baseline_7_45_", parameters(i)));
        load(strcat("stat_t_w_synch_with_baseline_7_45_", parameters(i)));
    else
        load(strcat("stat_p_t_snych_without_baseline_7_45_", parameters(i)));
        load(strcat("stat_p_w_synch_without_baseline_7_45_", parameters(i)));
        load(strcat("stat_t_w_synch_without_baseline_7_45_", parameters(i)));
    end
    
    % open the cluster based statistics
    if tobaseline == "yes"
        load(strcat("stat_p_t_snych_with_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_p_w_synch_with_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_t_w_synch_with_baseline_7_45_cluster_", parameters(i)));
    else
        load(strcat("stat_p_t_snych_without_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_p_w_synch_without_baseline_7_45_cluster_", parameters(i)));
        load(strcat("stat_t_w_synch_without_baseline_7_45_cluster_", parameters(i)));
    end
           
    tcl = tiledlayout(3,3);  
    ax = gobjects(1,9);
    index=1;
    %loop through the phases
    for j = 1:3
        ax(index) = nexttile;
        
        % Plot the grand averaged data in the first column     
        % Some initial configurations
        phases=[phasic_all_synch; tonic_all_synch; wake_all_synch];
        freq=phasic_all_synch.freq;
        time=phasic_all_synch.time;
        font_size=10;
        x_label = 'Time (sec)';
        y_label = 'Frequency (Hz)';
         if tobaseline =="yes"
            z_lim=[0.8 1.6];
        else
           z_lim=[0.13 0.28];
        end
        x_lim=[-0.2 0.6];
        % Set up figure standards
        PS = PLOT_STANDARDS();
        fig1_comps.fig = gcf;
        colormap(ax(index), "jet") 
        cfg = [];
        % Pay attention to cfg.baseline, this is optional to use
        % cfg.baseline     = [-0.2 -0.05];
        %  cfg.baselinetype = 'relative';
        cfg.zlim=z_lim;
        cfg.xlim= x_lim;
        cfg.figure = 'gcf';
                fig1_comps.p1=ft_singleplotTFR(cfg, phases(j));  

        xline(0, '--w', 'LineWidth',2);
        c = colorbar;
        c.Label.String = 'PLV';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=font_size;
        if tobaseline =="yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i), phase_titles(j)), " Grand averaged - baseline corrected"}, 'fontsize', font_size);
        else
           fig1_comps.plotTitle=title({strcat(area_titles(i), phase_titles(j)), " Grand averaged"}, 'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
            
        % Step onto the next tile
        index=index+1;
        ax(index) = nexttile;
    
        % Plot the T values in the second column
        % Initialize some settings
        stat=[stat_p_t_synch stat_p_w_synch stat_t_w_synch];

        % Set up figure standards
        fig1_comps.fig = gcf;
        fig1_comps.p1=imagesc(squeeze(stat(j).time), squeeze(stat(j).freq), squeeze(stat(j).stat));
        colormap(ax(index), "jet")  
        set(gca,'YDir','normal')
        ax=gca;
        ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
        xline(0, '--k', 'LineWidth', 2)
        c = colorbar;
        c.Label.String = 'T value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        caxis([-3 4])
        if tobaseline == "yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "T values - baseline corrected"},'fontsize', font_size);
        else
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "T values"},'fontsize', font_size);
        end
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
        fig1_comps.p1=imagesc(stat(j).time, squeeze(stat(j).freq), squeeze(stat(j).prob));
        colormap(ax(index),"bone")
        set(gca,'YDir','normal')
        ax=gca;
        ax.XAxis.Limits = x_lim; % if you need the plot to have the same axis values as in the first column 
        xline(0, '--k', 'LineWidth', 2)
        caxis([0.001 0.05])
        c = colorbar;
        c.Label.String = 'p value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        if tobaseline == "yes"
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "p values - baseline corrected"},'fontsize', font_size);
        else
            fig1_comps.plotTitle=title({strcat(area_titles(i),stat_phase_titles(j)), "p values"},'fontsize', font_size);
        end
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    
        % Plot the significant clusters inside
        hold on
        clusterstats=[stat_p_t_synch_cluster stat_p_w_synch_cluster stat_t_w_synch_cluster];

        simple_boundary = bwboundaries(squeeze(clusterstats(j).prob)<0.05);
        if size(simple_boundary,1) >0
            for b=1:length(simple_boundary{1,1})
                d=simple_boundary{1,1};
                c = [0.8500 0.3250 0.0980]; %colour code
                h1=plot(stat(j).time(d(:,2)), stat(j).freq(d(:,1)), 'Color', c, 'LineWidth', 2); %FaceAlpha', 0.3)
            end
            set(h1, 'DisplayName', 'significant cluster');
            fig1_comps.plotLegend=legend(h1,'significant clusters', "Location","southoutside", "FontName", 'Times New Roman');
            clear h1 d
        end
        
    end
end


%% Segmented statistical analysis

clear all
%working specifically with non-baselinecorrected data
% Open the correct files
cd 'f:\ANT_HEP\HEP_synchronization_512_70hz'
load("patientdata_phasic_synch_without_baseline_7_45.mat");
load("patientdata_tonic_synch_without_baseline_7_45.mat");
load("patientdata_wake_synch_without_baseline_7_45.mat");


patients = {'BA', 'HaJu','KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
sleep_wake_phases = {'phasic', 'tonic','wake'};

% Save into a big structure
for row=1:length(sleep_wake_phases)
    for column=1:length(patients)
        if row==1
            patients_synch{1,column} =phasic{1,column};
        elseif row==2
            patients_synch{2,column}= tonic{1, column};
        else
            patients_synch{3,column}= wake{1, column};
        end
    end
end
clear column row


% Segment into 75 ms segments
baseline_indices=1:8; % baseline: -125 ms - -50ms
start=19; % first segment starts at 50ms
% Creating the segments
for i=1:floor((length(patients_synch{1,1}.time)-start)/length(baseline_indices))
    timepoints=patients_synch{1,1}.time;
    segments{i}=timepoints(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
    segment_indices{i}=(start+((i-1)*length(baseline_indices)):(start+(i*length(baseline_indices))-1));
end
baseline_window=timepoints(baseline_indices);


%Create an analysis for each scalp area: fr,ce,pa,scalpavg
parameters=["fr" "ce" "pa" "scalpavg"];

%loop through each area consecutively

for p=1:length(parameters) 
    % Run the statistical analysis on each segment
    for i=1:length(segments)
    
        % Create the structure for each segment(i)
    
        for j=1:length(sleep_wake_phases)
            % Cut every patient's data into segments
            % (new variable: phase x patient)
            for k=1:length(patients)
                act_segment{j,k}.label=patients_synch{j,k}.label;
                act_segment{j,k}.dimord=patients_synch{j,k}.dimord;
                act_segment{j,k}.freq=patients_synch{j,k}.freq;
                act_segment{j,k}.time=patients_synch{j,k}.time(1,segment_indices{i});
                act_segment{j,k}.powspctrm=patients_synch{j,k}.(parameters{1,p})(:,:,segment_indices{i});
            
                bl_segment{j,k}.label=patients_synch{j,k}.label;
                bl_segment{j,k}.dimord=patients_synch{j,k}.dimord;
                bl_segment{j,k}.freq=patients_synch{j,k}.freq;
                bl_segment{j,k}.time=patients_synch{j,k}.time(1,segment_indices{i}); %using this time interval, because the statistical analysis runs only with the same time points
                bl_segment{j,k}.powspctrm=patients_synch{j,k}.(parameters{1,p})(:,:,baseline_indices);
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
        cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
        cfg.clusteralpha = 0.05; 
        cfg.tail = 0; % -1, 1 or 0 (default = 0); one-sided or two-sided test 
        cfg.alpha = 0.05; 
        cfg.uvar = 1; 
        cfg.ivar = 2; 
    
        % Run the analysis
        stat_p_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{1,:}, bl_segment{1,:}); 
        stat_t_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{2,:}, bl_segment{2,:}); 
        stat_w_timefreq_seg{i} = ft_freqstatistics(cfg, act_segment{3,:}, bl_segment{3,:}); 
    
    end
    
    %Concatenate the segments 
    stats=[stat_p_timefreq_seg; stat_t_timefreq_seg; stat_w_timefreq_seg];
    
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
    stat=[unified{1,1}; unified{2,1}; unified{3,1}];
    titles=["Phasic time segments-baseline window", "Tonic time segments-baseline window", "Wake time segments-baseline window"];
    parameter_titles=["Frontal EEG & ANT synch. - ","Central EEG & ANT synch. - ", "Parietal EEG & ANT synch. - ", "Averaged EEG & ANT synch. - "];
    
    
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
            fig1_comps.plotTitle=title(strcat(parameter_titles(p), titles(i)), 'fontsize', font_size);
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
    for i=1:3
        fig1_comps.fig = gcf;
        subplot(3,1,i)
        fig1_comps.p1=imagesc(stat(i).time, squeeze(stat(i).freq), squeeze(stat(i).stat));
        set(gca,'YDir','normal')
        xline(0, '--w', 'LineWidth',2)
        caxis([-4 4])
        c = colorbar;
        c.Label.String = 'T value';
        c.Label.FontName= 'Times New Roman';
        c.Label.FontSize=10;
        fig1_comps.plotTitle=title(strcat(parameter_titles(p), titles(i)),'fontsize', font_size);
        fig1_comps.plotXLabel = xlabel(x_label) ;
        fig1_comps.plotYLabel =ylabel(y_label);
        set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
        set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', font_size, 'FontWeight' , 'bold');
    end
end
