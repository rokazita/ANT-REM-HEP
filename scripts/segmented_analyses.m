%% TFA segmented statistical analysis

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

%% ITC Segmented statistical analysis

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


%% SYNCHRONIZATION Segmented statistical analysis

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
