%% HEP analysis
%% creating the structures for hep analysis
ft_defaults
clc;
clear all;
close all;
%EÉ excluded
patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};

% loop_version = {'phasic', 'tonic', 'wake'};
loop_version = {'nrem'};


%how many trials we have
alltrials=zeros(size(loop_version,2), size(patients,2));

for i=1:length(loop_version)
    for j=1:length(patients)
        cd('d:\ANT_HEP\HEP_epochdata_512_70hz')   % PATH!
        %load epoch files
        load(strcat(char(patients(j)),'_clean_', char(loop_version(i)),'_hep_epochs.mat'));
        %load the clean mat file - for the channel labels !!!! the corresponding one
        cd (strcat('d:\ANT_HEP\HEP_segments_512_70hz\',char(loop_version(i))));   % PATH!
        clean=load(strcat(char(patients(j)),'_clean_', char(loop_version(i)),'.mat'));
        %create for hep analyses in fieldtrip = timelock analyses
        hep=[];
        hep.epochs=epochdata;
        %hep.avg=mean(epochdata,3);
        for k=1:size(hep.epochs,1)
            hep.avg(k,:)=mean(epochdata(k,:,:),3);
        end
        hep.var = var(epochdata,0,3); 
        hep.fsample=512;  %!!!!!512 for downsampled or 1024 for FSI or 2048 for others
        fns=fieldnames(clean);
        hep.label=clean.(fns{1}).label;
        hep.dimord='chan_time';
        hep.time=linspace(-0.2,0.8,512); %!!!!!ADJUST SAMPLING RATE
        hep_all{1,j}=hep;
        alltrials(i,j)=size(epochdata,3);
    end
    cd('d:\ANT_HEP\HEP_epochdata_512_70hz')   % PATH!
    save(strcat('hep_all_', char(loop_version(i))), 'hep_all', '-v7.3');
end
%save(strcat('hep_all_trial_numbers'), 'alltrials', '-v7.3');

%latex_table = latex(sym(alltrials));


%% Correction of inverse signals:
cd('d:\ANT_HEP\HEP_epochdata_512_70hz')   % PATH!

patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
%loop_version = {'phasic', 'tonic','wake'};
loop_version = {'nrem'};

clear i j
for i = 1:length(loop_version)
    load(strcat('hep_all_', char(loop_version(i))), 'hep_all');
    every_hep_all(i,:) = hep_all;
end

%KB(25-{'EEG BTh0-BTh1'}), 
%MAFE (25-{'EEG BTH1-BTH2'}),
%PJ(27-{'EEG JTh10-JTh11'}), 
%TöTa(25-{'EEG JTh10-JTh11'}), 
%ToZa (26-{'EEG Bth2-Bth3'}, 28-{'EEG Jth9-Jth10'})
for k=1:size(loop_version,2)
    every_hep_all{k,3}.avg(25,:)=-every_hep_all{k,3}.avg(25,:);
    every_hep_all{k,3}.var(25,:)=-every_hep_all{k,3}.var(25,:);
end

for k=1:size(loop_version,2)
    every_hep_all{k,5}.avg(25,:)=-every_hep_all{k,5}.avg(25,:);
    every_hep_all{k,5}.var(25,:)=-every_hep_all{k,5}.var(25,:);
end

for k=1:size(loop_version,2)
    every_hep_all{k,7}.avg(27,:)=-every_hep_all{k,7}.avg(27,:);
    every_hep_all{k,7}.var(27,:)=-every_hep_all{k,7}.var(27,:);
end

for k=1:size(loop_version,2)
    every_hep_all{k,9}.avg(25,:)=-every_hep_all{k,9}.avg(25,:);
    every_hep_all{k,9}.var(25,:)=-every_hep_all{k,9}.var(25,:);
end

for k=1:size(loop_version,2)
    every_hep_all{k,10}.avg(26,:)=-every_hep_all{k,10}.avg(26,:);
    every_hep_all{k,10}.var(26,:)=-every_hep_all{k,10}.var(26,:);
end

for k=1:size(loop_version,2)
    every_hep_all{k,10}.avg(28,:)=-every_hep_all{k,10}.avg(28,:);
    every_hep_all{k,10}.var(28,:)=-every_hep_all{k,10}.var(28,:);
end
%rename Piri's 0 line
for k=1:size(loop_version,2)
    every_hep_all{k,6}.label(20,:)={'ref_0'};
end

save('nrem_every_hep_all_without_baselinecorr', 'every_hep_all', '-v7.3');

%baseline correction
for j = 1:length(patients)
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        cfg=[];
        cfg.baseline     = [-0.2 -0.05];
        cfg.channel      = index;
        every_hep_all{k,j} = ft_timelockbaseline(cfg, every_hep_all{k,j});
    end
end


save('nrem_every_hep_all_with_baselinecorr', 'every_hep_all', '-v7.3');


%% HEP ANALYSIS VISUALIZATION individually ANT CHANNEL

clear all
%with baseline
% load('f:\ANT_HEP\HEP_epochdata_512_70hz\every_hep_all_with_baselinecorr.mat')
load('d:\ANT_HEP\HEP_epochdata_512_70hz\nrem_every_hep_all_with_baselinecorr.mat')

%without baselinecorr
load('f:\ANT_HEP\HEP_epochdata_512_70hz\every_hep_all_without_baselinecorr.mat')

patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
%loop_version = {'phasic', 'tonic','wake'};
loop_version = {'nrem'};

%microsates individually

%first
 PS = PLOT_STANDARDS();
 figure;
 fig1_comps.fig = gcf;
 ind = reshape(1:length(patients)*length(loop_version),length(patients), length(loop_version)).';
 count=1;
clear k j
 for j= 1:length(patients)
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        time_hep = every_hep_all{k,j}.time;
        if k==1
            mean_phasic_hep = every_hep_all{1,j}.avg(index,:);
            subplot(length(loop_version),length(patients),ind(count));
            count=count+1;
            fig1_comps.p1=plot(time_hep,mean_phasic_hep);
            %fig1_comps.plotLegend=legend('phasic');
            hold on
        elseif k==2
            mean_tonic_hep = every_hep_all{2,j}.avg(index,:);
            subplot(length(loop_version),length(patients),ind(count));
            count=count+1;
            fig1_comps.p2=plot(time_hep,mean_tonic_hep);
            %fig1_comps.plotLegend=legend('tonic');
            hold on
        else
            mean_wake_hep = every_hep_all{3,j}.avg(index,:);
            subplot(length(loop_version),length(patients),ind(count));
            count=count+1;
            fig1_comps.p3=plot(time_hep,mean_wake_hep);
            %fig1_comps.plotLegend=legend('wake');
            hold on
        end
        fig1_comps.plotTitle=title(strcat(patients(j), '- ANT'), 'fontsize', 10);
       fig1_comps.plotXLabel = xlabel('Time (sec)') ;
       fig1_comps.plotYLabel =ylabel('Amplitude (μV)');
       set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
       set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 10, 'FontWeight' , 'bold');
   % fig1_comps.plotLegend =legend(strcat(loop_version(k)), 'location','northeast','orientation','horizontal');
    end
      set(fig1_comps.p1, 'Color', PS.Blue4, 'LineWidth', 0.6);
      %set(fig1_comps.p2, 'Color', PS.Green4, 'LineWidth', 0.6);
      %set(fig1_comps.p3, 'Color', PS.Red4, 'LineWidth', 0.6);
 end
% fig1_comps.plotLegend =legend('phasic', 'tonic','wake', 'location','northeast','orientation','horizontal');
 %set(fig1_comps.plotLegend, 'FontName', 'Times New Roman');
%saveas(gcf,'microstates_individually','epsc')

clearvars -except patients loop_version every_hep_all
 
%microstates in one figure
 PS = PLOT_STANDARDS();
 figure;
 fig1_comps.fig = gcf;
 count=1;
 for j= 1:length(patients)
    subplot(5,3,j)  %numbers to be overwritten
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        time_hep = every_hep_all{k,j}.time;
        if k==1
            mean_phasic_hep = every_hep_all{1,j}.avg(index,:);
            fig1_comps.p1=plot(time_hep,mean_phasic_hep);
            hold on
        elseif k==2
            mean_tonic_hep = every_hep_all{2,j}.avg(index,:);
            fig1_comps.p2=plot(time_hep,mean_tonic_hep);
            hold on
        else
            mean_wake_hep = every_hep_all{3,j}.avg(index,:);
            fig1_comps.p3=plot(time_hep,mean_wake_hep);
            hold on
        end
        %used only for correct legend labels
%         L=plot(nan,nan, 'b');
%         M=plot(nan,nan, 'r');
%         G=plot(nan,nan, 'g');
       fig1_comps.plotTitle=title(strcat(patients(j), '- ANT'));
       fig1_comps.plotXLabel = xlabel('Time (sec)') ;
       fig1_comps.plotYLabel =ylabel('Amplitude (μV)');
       set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman');
       set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 10, 'FontWeight' , 'bold');
%         legend([L(1), M(1), G(1)], 'phasic', 'tonic', 'wake');
    end
     set(fig1_comps.p1, 'Color', PS.Blue4, 'LineWidth', 0.6);
     %set(fig1_comps.p2, 'Color', PS.Green4, 'LineWidth', 0.6);
     %set(fig1_comps.p3, 'Color', PS.Red4, 'LineWidth', 0.6);
 end
 
 
 %% GRAND AVERAGE - ANT CHANNEL
ft_defaults
clc;
clear all;
close all;

%with baseline
load('d:\ANT_HEP\HEP_epochdata_512_70hz\every_hep_all_with_baselinecorr.mat')
load('d:\ANT_HEP\HEP_epochdata_512_70hz\nrem_every_hep_all_with_baselinecorr.mat')
%without baselinecorr
load('d:\ANT_HEP\HEP_epochdata_512_70hz\every_hep_all_without_baselinecorr.mat')
load('d:\ANT_HEP\HEP_epochdata_512_70hz\nrem_every_hep_all_without_baselinecorr.mat')


patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
loop_version = {'phasic', 'tonic','wake'};
loop_version = {'nrem'};


%----renaming ANT-channel
for j = 1:length(patients)
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        if length(index) == 1
            every_hep_all{k,j}.label{index,1} = 'Averagedant';
        elseif length(index) == 2
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:))/2;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant'; 
        else 
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            every_hep_all{k,j}.label{index(3,1),1} = 'General-ANT3';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:)+every_hep_all{k,j}.avg(index(3,1),:))/3;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant';
        end
    end
end

cfg=[];
cfg.channel={'Averagedant'};
grandaverage_phasic = ft_timelockgrandaverage(cfg, every_hep_all{1,:});
grandaverage_tonic = ft_timelockgrandaverage(cfg, every_hep_all{2,:});
grandaverage_wake = ft_timelockgrandaverage(cfg, every_hep_all{3,:});
grandaverage_nrem = ft_timelockgrandaverage(cfg, every_hep_all{1,:});




%save
save('grandaverage_phasic_without_baselinecorr', 'grandaverage_phasic', '-v7.3');
save('grandaverage_tonic_without_baselinecorr', 'grandaverage_tonic', '-v7.3');
save('grandaverage_wake_without_baselinecorr', 'grandaverage_wake', '-v7.3');
save('grandaverage_nrem_without_baselinecorr', 'grandaverage_nrem', '-v7.3');


load('grandaverage_nrem_with_baselinecorr.mat');
load('grandaverage_phasic_with_baselinecorr.mat');
load('grandaverage_tonic_with_baselinecorr.mat');
load('grandaverage_wake_with_baselinecorr.mat');


%plotting
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;


 fig1_comps.p1=plot(grandaverage_phasic.time, grandaverage_phasic.avg, 'DisplayName', char(grandaverage_phasic.label) + ", phasic");
hold on
 fig1_comps.p2=plot(grandaverage_tonic.time, grandaverage_tonic.avg, 'DisplayName', char(grandaverage_tonic.label) + ", tonic");
hold on
 fig1_comps.p3=plot(grandaverage_wake.time, grandaverage_wake.avg, 'DisplayName', char(grandaverage_wake.label) + ", wake");
 hold on
 fig1_comps.p4=plot(grandaverage_nrem.time, grandaverage_nrem.avg, 'DisplayName', char(grandaverage_nrem.label) + ", nrem");
hold on
fig1_comps.plotLegend =legend('Fázisos REM','Tónusos REM','Ébrenlét','NREM');
 hold on

% Shade area between mean and standard deviation
fill([grandaverage_phasic.time, fliplr(grandaverage_phasic.time)], [grandaverage_phasic.avg + sqrt(grandaverage_phasic.var), fliplr(grandaverage_phasic.avg - sqrt(grandaverage_phasic.var))], PS.Blue4, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
fill([grandaverage_tonic.time, fliplr(grandaverage_tonic.time)], [grandaverage_tonic.avg + sqrt(grandaverage_tonic.var), fliplr(grandaverage_tonic.avg - sqrt(grandaverage_tonic.var))], PS.Red4, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
fill([grandaverage_wake.time, fliplr(grandaverage_wake.time)], [grandaverage_wake.avg + sqrt(grandaverage_wake.var), fliplr(grandaverage_wake.avg - sqrt(grandaverage_wake.var))], PS.Yellow4, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');
fill([grandaverage_nrem.time, fliplr(grandaverage_nrem.time)], [grandaverage_nrem.avg + sqrt(grandaverage_nrem.var), fliplr(grandaverage_nrem.avg - sqrt(grandaverage_nrem.var))], PS.Purple4, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

xline(0, '--k', 'LineWidth',2, 'handlevisibility', 'off')
fig1_comps.plotTitle=title('HEP analízis - ANT csatorna');
fig1_comps.plotXLabel = xlabel('Idő (sec)');
fig1_comps.plotYLabel =ylabel({'Amplitúdó (μV)'});

legendX0 = .78; legendY0 = .82; legendWidth = .05; legendHeight = .04;
set(fig1_comps.plotLegend, 'position', [legendX0, legendY0, legendWidth, ...
    legendHeight], 'Box', 'on');
set(fig1_comps.plotLegend, 'FontSize', 10, 'LineWidth', 0.5);

set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman','FontSize', 14);
set(fig1_comps.plotTitle, 'FontName', 'Times New Roman','FontSize', 18, 'FontWeight' , 'bold');
set(fig1_comps.p1, 'Color', PS.Blue4, 'LineWidth', 1.3);
set(fig1_comps.p2, 'Color', PS.Red4, 'LineWidth',1.3);
set(fig1_comps.p3, 'Color', PS.Yellow4, 'LineWidth', 1.3);
set(fig1_comps.p4, 'Color', PS.Purple4, 'LineWidth', 1.3);

%ylim([-2.5 4]);

%plot the clusters too (both for baselined and nonbaselined
%with baseline
 load('d:\ANT_HEP\HEP_epochdata_512_70hz\ant_hep_stat_sigcluster_times_baselinecorrected.mat')
load('d:\ANT_HEP\HEP_epochdata_512_70hz\nrem_ant_hep_stat_sigcluster_times_baselinecorrected.mat')
%without baselinecorr
%load('f:\ANT_HEP\HEP_epochdata_512_70hz\ant_hep_stat_sigcluster_times_nonbaselinecorrected.mat')


for i=1:size(stat_p_n.time,2)
    if stat_p_n.prob(i) < 0.05
       [row, col]=find(grandaverage_phasic.time==stat_p_n.time(i))       ;
       hold on
       fig1_comps.p5=plot(grandaverage_phasic.time(col), -0.5, 'k.', 'handlevisibility', 'off');
       set(fig1_comps.p5, 'Color', PS.MyBlack, 'MarkerSize',10);
    end
    if stat_p_t.prob(i) < 0.05
      [row, col]=find(grandaverage_phasic.time==stat_p_t.time(i))       ;
      hold on
      fig1_comps.p6=plot(grandaverage_phasic.time(col), -2, 'r.' );    
     set(fig1_comps.p6, 'Color', PS.Purple4, 'MarkerSize',10);
    end
    if stat_t_n.prob(i) < 0.05
       [row, col]=find(grandaverage_phasic.time==stat_t_n.time(i))       ;
       hold on
       fig1_comps.p7=plot(grandaverage_phasic.time(col), -1, '.', 'handlevisibility', 'off');
       set(fig1_comps.p7, 'Color', PS.Red4, 'MarkerSize',10);
    end
    if stat_p_w.prob(i) < 0.05 
      [row, col]=find(grandaverage_phasic.time==stat_p_t.time(i))       ;
      hold on
      fig1_comps.p8=plot(grandaverage_phasic.time(col), -2.2, '.' );
      set(fig1_comps.p8, 'Color', PS.Yellow1, 'MarkerSize',10);
   end
end
legend boxoff 
a=axes('position',get(gca,'position'),'visible','off');
legend(a, [fig1_comps.p5 fig1_comps.p6 fig1_comps.p7 fig1_comps.p8], 'Fázisos vs. NREM', 'Fázisos vs Tónusos', 'Tónusos vs. NREM', 'Fázisos vs. Ébrenlét','Location','southeast', 'FontSize', 10)






%% próba

load('grandaverage_nrem_with_baselinecorr.mat');
load('grandaverage_phasic_with_baselinecorr.mat');
load('grandaverage_tonic_with_baselinecorr.mat');
load('grandaverage_wake_with_baselinecorr.mat');


% Plotting
PS = PLOT_STANDARDS();
figure;
fig1_comps.fig = gcf;

% Main plots
fig1_comps.p1 = plot(grandaverage_phasic.time, grandaverage_phasic.avg, 'DisplayName', char(grandaverage_phasic.label) + ", phasic");
hold on;
fig1_comps.p2 = plot(grandaverage_tonic.time, grandaverage_tonic.avg, 'DisplayName', char(grandaverage_tonic.label) + ", tonic");
hold on;
% fig1_comps.p3 = plot(grandaverage_wake.time, grandaverage_wake.avg, 'DisplayName', char(grandaverage_wake.label) + ", wake");
% hold on;
fig1_comps.p4 = plot(grandaverage_nrem.time, grandaverage_nrem.avg, 'DisplayName', char(grandaverage_nrem.label) + ", nrem");
hold on;

fig1_comps.plotLegend = legend('Phasic REM','Tonic REM','NREM');
hold on;

% Shade area between mean and standard deviation
fill([grandaverage_phasic.time, fliplr(grandaverage_phasic.time)], [grandaverage_phasic.avg + sqrt(grandaverage_phasic.var), fliplr(grandaverage_phasic.avg - sqrt(grandaverage_phasic.var))], PS.Blue4, 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'HandleVisibility', 'off');
fill([grandaverage_tonic.time, fliplr(grandaverage_tonic.time)], [grandaverage_tonic.avg + sqrt(grandaverage_tonic.var), fliplr(grandaverage_tonic.avg - sqrt(grandaverage_tonic.var))], PS.Red4, 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'HandleVisibility', 'off');
% fill([grandaverage_wake.time, fliplr(grandaverage_wake.time)], [grandaverage_wake.avg + sqrt(grandaverage_wake.var), fliplr(grandaverage_wake.avg - sqrt(grandaverage_wake.var))], PS.Yellow4, 'EdgeColor', 'none', 'FaceAlpha', 0.1, 'HandleVisibility', 'off');
fill([grandaverage_nrem.time, fliplr(grandaverage_nrem.time)], [grandaverage_nrem.avg + sqrt(grandaverage_nrem.var), fliplr(grandaverage_nrem.avg - sqrt(grandaverage_nrem.var))], PS.Purple4, 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'HandleVisibility', 'off');

xline(0, '--k', 'LineWidth', 2, 'HandleVisibility', 'off');
fig1_comps.plotTitle = title('HEP analysis - ANT channel');
fig1_comps.plotXLabel = xlabel('Time (sec)');
fig1_comps.plotYLabel = ylabel({'Amplitude (μV)'});

legendX0 = .78; legendY0 = .82; legendWidth = .05; legendHeight = .04;
set(fig1_comps.plotLegend, 'position', [legendX0, legendY0, legendWidth, legendHeight], 'Box', 'on');
set(fig1_comps.plotLegend, 'FontSize', 10, 'LineWidth', 0.5);

set([fig1_comps.plotXLabel, fig1_comps.plotYLabel], 'FontName', 'Times New Roman', 'FontSize', 14);
set(fig1_comps.plotTitle, 'FontName', 'Times New Roman', 'FontSize', 18, 'FontWeight', 'bold');
set(fig1_comps.p1, 'Color', PS.Blue4, 'LineWidth', 1.3);
set(fig1_comps.p2, 'Color', PS.Red4, 'LineWidth', 1.3);
% set(fig1_comps.p3, 'Color', PS.Yellow4, 'LineWidth', 1.3);
set(fig1_comps.p4, 'Color', PS.Purple4, 'LineWidth', 1.3);

%ylim([-2.5 4]);

% Plot the clusters too (both for baselined and nonbaselined)
% with baseline
load('d:\ANT_HEP\HEP_epochdata_512_70hz\ant_hep_stat_sigcluster_times_baselinecorrected.mat')
load('d:\ANT_HEP\HEP_epochdata_512_70hz\nrem_ant_hep_stat_sigcluster_times_baselinecorrected.mat')

% without baselinecorr
% load('f:\ANT_HEP\HEP_epochdata_512_70hz\ant_hep_stat_sigcluster_times_nonbaselinecorrected.mat')

for i = 1:size(stat_p_n.time, 2)
    if stat_p_n.prob(i) < 0.05
        [row, col] = find(grandaverage_phasic.time == stat_p_n.time(i));
        hold on;
        fig1_comps.p5 = plot(grandaverage_phasic.time(col), -1.5, 'k.', 'HandleVisibility', 'off');
        set(fig1_comps.p5, 'Color', PS.MyBlack, 'MarkerSize', 10);
    end
    if stat_p_t.prob(i) < 0.05
        [row, col] = find(grandaverage_phasic.time == stat_p_t.time(i));
        hold on;
        fig1_comps.p6 = plot(grandaverage_phasic.time(col), -1.7, 'r.', 'HandleVisibility', 'off');
        set(fig1_comps.p6, 'Color', PS.Blue4, 'MarkerSize', 10);
    end
    if stat_t_n.prob(i) < 0.05
        [row, col] = find(grandaverage_phasic.time == stat_t_n.time(i));
        hold on;
        fig1_comps.p7 = plot(grandaverage_phasic.time(col), -1.9, '.', 'HandleVisibility', 'off');
        set(fig1_comps.p7, 'Color', PS.Red4, 'MarkerSize', 10);
    end
%     if stat_p_w.prob(i) < 0.05
%         [row, col] = find(grandaverage_phasic.time == stat_p_t.time(i));
%         hold on;
%         fig1_comps.p8 = plot(grandaverage_phasic.time(col), -2.1, '.', 'HandleVisibility', 'off');
%         set(fig1_comps.p8, 'Color', PS.Green4, 'MarkerSize', 10);
%     end
end

legend boxoff;
a = axes('position', get(gca, 'position'), 'visible', 'off');
% legend(a, [fig1_comps.p5 fig1_comps.p6 fig1_comps.p7 fig1_comps.p8], 'Phasic vs. NREM', 'Phasic vs Tonic', 'Tonic vs. NREM', 'Phasic vs. Awake', 'Location', 'southeast', 'FontSize', 10);
% legend(a, [fig1_comps.p6 fig1_comps.p8], 'Phasic vs Tonic',  'Phasic vs. Awake', 'Location', 'southeast', 'FontSize', 10);
legend(a, [fig1_comps.p5 fig1_comps.p6 fig1_comps.p7], 'Phasic vs. NREM', 'Phasic vs Tonic', 'Tonic vs. NREM', 'Location', 'southeast', 'FontSize', 10);


%% STATISTICAL ANALYSIS
ft_defaults
clc;
clear all;
close all;

%without baselinecorr
cd 'd:\ANT_HEP\HEP_epochdata_512_70hz'

patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
% loop_version = {'phasic', 'tonic','wake'};

nrem=load('nrem_every_hep_all_with_baselinecorr.mat');
rem=load('every_hep_all_with_baselinecorr.mat');
%or
nrem=load('nrem_every_hep_all_without_baselinecorr.mat');
load('every_hep_all_without_baselinecorr.mat')

for i=1:length(patients)
    rem.every_hep_all(4,i)=nrem.every_hep_all(1,i);
end
loop_version = {'phasic', 'tonic','wake','nrem'};
for j=1:length(loop_version)
    for i=1:length(patients)
        every_hep_all(j,i)=rem.every_hep_all(j,i);
    end
end


%----renaming ANT-channel
for j = 1:length(patients)
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        if length(index) == 1
            every_hep_all{k,j}.label{index,1} = 'Averagedant';
        elseif length(index) == 2
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:))/2;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant'; 
        else 
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            every_hep_all{k,j}.label{index(3,1),1} = 'General-ANT3';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:)+every_hep_all{k,j}.avg(index(3,1),:))/3;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant';
        end
    end
end



% %-----baseline correction
% 
% for j = 1:length(patients)
%     for k = 1:length(loop_version)
%         cfg=[];
%         cfg.baseline     = [-0.2 -0.05];
%         cfg.channel      = 'General-ANT*';
%         every_hep_all{k,j} = ft_timelockbaseline(cfg, every_hep_all{k,j});
%     end
% end

%design matrix
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


latency = [0.35 0.65];
cfg = []; 
cfg.design = design;
cfg.channel = 'Averagedant'; 
cfg.latency = latency;
%cfg.parameter='averaged';
%cfg.latency = 'all'; 
cfg.method = 'montecarlo'; % use the Monte Carlo Method to calculate the significance probability 
cfg.statistic = 'ft_statfun_depsamplesT';  %paired samples t test
cfg.clusteralpha = 0.05; 
cfg.correctm = 'cluster'; 
cfg.clusterstatistic = 'maxsum'; 
% cfg.minnbchan = 1; 
cfg.numrandomization = 'all'; 
cfg.tail = 0;
cfg.clustertail = 0; 
cfg.alpha = 0.05; 
cfg.uvar = 1; 
cfg.ivar = 2; 
 cfg.neighbours = []; 

%cfg.clusterthreshold = 'nonparametric_individual';
% 
[stat_p_n] = ft_timelockstatistics(cfg, every_hep_all{1,:}, every_hep_all{4,:}); %phasic vs nrem
[stat_t_n] = ft_timelockstatistics(cfg, every_hep_all{2,:}, every_hep_all{4,:}); %tonic vs nrem
[stat_w_n] = ft_timelockstatistics(cfg, every_hep_all{3,:}, every_hep_all{4,:}); %wake vs nrem

% [stat_p_t] = ft_timelockstatistics(cfg, every_hep_all{1,:}, every_hep_all{2,:}); %phasic vs tonic
% [stat_p_w] = ft_timelockstatistics(cfg, every_hep_all{1,:}, every_hep_all{3,:}); 
% [stat_t_w] = ft_timelockstatistics(cfg, every_hep_all{2,:}, every_hep_all{3,:}); 

%there are significant clusters both with and without baselinecorrection

% 
% [stat_p_t] = ft_timelockstatistics(cfg, grandaverage_phasic,grandaverage_tonic); %phasic vs tonic
% [stat_p_w] = ft_timelockstatistics(cfg, grandaverage_phasic, grandaverage_wake); 
% [stat_t_w] = ft_timelockstatistics(cfg, grandaverage_tonic, grandaverage_wake); 
% 
% grandaverage_phasic.label={'ANT'}
% grandaverage_tonic.label={'ANT'}

%there are clusters when baselinecorrection is applied
j=1;
k=1;
for i=1:size(stat_t_n.prob,2)
    if stat_p_n.prob(i) < 0.05
          p_n_ido(j)=stat_p_n.time(i);
          j=j+1;
    end
     if stat_t_n.prob(i) < 0.05
          t_n_ido(k)=stat_t_n.time(i);
          k=k+1;
    end
end
save("nrem_ant_hep_stat_sigcluster_times_baselinecorrected", "stat_p_n", "stat_t_n", "p_n_ido", "t_n_ido")


%% average those parts that are significantly different
cd D:\ANT_HEP\HEP_epochdata_512_70hz
load('d:\ANT_HEP\HEP_epochdata_512_70hz\every_hep_all_with_baselinecorr.mat')


patients = {'BA', 'HaJu', 'KB', 'KEA', 'MaFe', 'PiRi', 'PJ', 'TI', 'TöTa', 'ToZa', 'FSI'};
loop_version = {'phasic', 'tonic','wake'};


%----renaming ANT-channel
for j = 1:length(patients)
    for k = 1:length(loop_version)
        ANT_channel = strfind(lower(every_hep_all{k,j}.label),'th');
        index= find(~cellfun(@isempty, ANT_channel));
        if length(index) == 1
            every_hep_all{k,j}.label{index,1} = 'Averagedant';
        elseif length(index) == 2
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:))/2;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant'; 
        else 
            every_hep_all{k,j}.label{index(1,1),1} = 'General-ANT1';
            every_hep_all{k,j}.label{index(2,1),1} = 'General-ANT2';
            every_hep_all{k,j}.label{index(3,1),1} = 'General-ANT3';
            channum=size(every_hep_all{k,j}.avg,1)+1;
            every_hep_all{k,j}.avg(channum,:)=(every_hep_all{k,j}.avg(index(1,1),:)+every_hep_all{k,j}.avg(index(2,1),:)+every_hep_all{k,j}.avg(index(3,1),:))/3;
            every_hep_all{k,j}.label{channum,1} = 'Averagedant';
        end
    end
end


for i=1:11
    for j=1:size(every_hep_all{1, i}.label,1)
        if(every_hep_all{1, i}.label(j,1)=="Averagedant")
            % phasic-tonic cluster (424.3–529.9 ms) time:319-374 index
            phasic_pt(i)=mean(every_hep_all{1, i}.avg(j,319:374));
            tonic_pt(i)=mean(every_hep_all{2, i}.avg(j,319:374));
            %phasic-wake cluster (373.4-488.8 ms) time: 294-353 index
            phasic_pw(i)=mean(every_hep_all{1, i}.avg(j,294:353));
            wake_pw(i)=mean(every_hep_all{3, i}.avg(j,294:353));
        end
    end
end




