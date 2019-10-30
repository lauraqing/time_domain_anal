%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qing Wu Phase kurtosis calculation on labdata_pf     %
% v1: 2013-11-11                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load impulse noise waveform %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% constant definition
    p0=2e-5;             % referrence pressure p0
    dis1=2e-5;           % Microphone 20usec
    dis2=5e-1;           % Microphone 20usec
    step=6;              % To find the A-duration length by intersection
    Fs=44100;            % Sampling Frequency of energy caculation
    t_tot=33.11;         % ave_A_dur=sum(a_dur)/10=9.46ms; t+/t-=1/2.5; t_tot=33.11ms
    v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];

    % use execl to draw durations 
    filename = 'kurtosis_freq_domain_test.xlsx';
    sheet=1; % classic defined kurtosis
%     sheet2=2; % time cycle for classic defined kurtosis in time domain
    % acoustic defined kurtosis
    xlRange='B2';

for i=1:20
    filepath=(['C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\',num2str(i),'\']);
    %% Ininitial the temp variances
    peak=0;n1=0;n2=0;x_0=0;abs_x_0=0;n0=0;n0_0=0;a=0;k=0;x_01=0;temp_neg=0;a_dur1=zeros(10,1);
    n00=0;a_dur2=zeros(10,1);n0_1=0;n_half=0;M=0;temp_com_dB=0;temp_peak=zeros(1,10);
    N=zeros(10,1);N1=0;N2=0;ave_x_kur=0;min_x_kur=0;max_x_kur=0;tem_max_x_kur=0;
    kurtosis_x_imp=zeros(10,1);dif_x_kur=0;num_norm=0;cons_ave_norm=0;
    for j=1:10       
    %% Load the data of impulse noise
       if(j~=10)
          x_imp_temp(:,:,j)=load([filepath,'testdata_00',num2str(j),'.lvm']); %test_00i.lvm
       else
          x_imp_temp(:,:,j)=load([filepath,'testdata_0',num2str(j),'.lvm']); %test_00i.lvm 
       end 
    x_imp(:,j)=x_imp_temp(:,2,j);
    temp_peak(j)=max(x_imp(:,j)); % Pick out all the peaks  
    L=length(x_imp);              %L is max length of impulse noise x 
    dB(j)=20*log10(max(x_imp(:,j)/(p0)));
    x(:,j)=x_imp(:,j);
    %% FFT
    n=4000;
    N_max=65536;
    % Consider the total secend if is 1s
%     t_tot=1;
%     t=[0:t_tot/Fs:n*t_tot/Fs];    % 0.0610 s
%     F=(Fs/n)*([1:n]-1); % calculate the real frequency used for FFT of 4000 points
    % Total time if is 6.683s
%     T_tot=6.683;                    % unit is 's'
    T_tot=1;                    % unit is 's'
    t_tot=T_tot*n*1000/N_max;            % window, unit is 'ms'
    Fs=N_max/T_tot;
    t=[0:1/Fs:n*1000/Fs];  
    F=(Fs/n)*([1:n]-1)/1000;             % calculate the real frequency used for FFT of 4000 points, unit in 1K
     
    n_step=16;          
    % defined a kurtosis wave in phase domain
        for h=1:n/n_step
           x_kur=x((h-1)*n_step+1:h*n_step,j);     
           y=x_kur;
           y_fft = fft(y, n_step);                      % small window FFT
           Y =abs(y_fft);
           Y_temp(:,h)=Y/(n_step/2);   
        end   
        for l=1:n_step
            Y_fdk_ave_temp=sum(Y_temp(l,:))/(n/n_step);
            dif_Y_fdk_temp=(Y_temp(l,:)-Y_fdk_ave_temp);
            Y_kdf_temp(l)=((sum((dif_Y_fdk_temp).^4))*((n/n_step)-1))/((sum(dif_Y_fdk_temp.^2)).^2); 
%         % acoustic kurtosis
%         kurtosis_x_imp1(j)=((sum((Y_temp).^4))*(L-1))/((sum(Y_temp.^2)).^2);   
        end
% % Debug plot the phase kurotsis: correct WQ. 2013-11-11    
%     figure;plot(F(1:n_step/2),Y_kdf_temp(1:n_step/2));
%     xlabel('KHz');title('Phase kurtosis in 61ms window');ylabel('Frequency Domian Kurtosis K(w)');
    %% Kurtosis calculation in frequenct domain use the same segamentation as classic kurtosis  
    %% phase classic kurtosis
    Y_kdf(:,j)=Y_kdf_temp;    
    end
    % Phase kurtosis calculation in 100 data by probability idea
        for p=1:n_step
            Y_fdk_ave(p,i)=sum(Y_kdf(p,:))/10;
            std_Y_fdk_ave(p,i)=std(Y_kdf(p,:));
        end
%     ave_kur_x_fdk_acoustic(i)=sum(kurtosis_x_imp1)/10;
%     std_ave_kur_x_fdk_acoustic(i)=std(ave_kur_x_fdk_acoustic);
% % Debug plot the phase kurotsis in statistics method: correct!~ QW 2013-11-11 
%     figure;plot(F(1:n_step/2),Y_fdk_ave(1:n_step/2,i));
%     xlabel('KHz');title('Phase kurtosis in 61ms window');ylabel('Frequency Domian Kurtosis K(w)');
%     figure;errorbar(F(1:n_step/2),Y_fdk_ave(1:n_step/2,i),std_Y_fdk_ave(1:n_step/2,i));
%     xlabel('KHz');title('Phase kurtosis in 61ms window');ylabel('Frequency Domian Kurtosis K(w)');

% write data into excel table only draw average peak vs voltage
% Phase defined classic kurtosis
    x_wrt_k=[Y_fdk_ave(:,i),std_Y_fdk_ave(:,i)]; 
    xlswrite(filename,x_wrt_k,sheet,xlRange);     % classic defined kurtosis
    sheet=sheet+1;
end
% % phase acoustic kurtosis
% x_wrt_k=[v',ave_kur_x_fdk_acoustic',std_ave_kur_x_fdk_acoustic']; 
% xlswrite(filename,x_wrt_k,sheet2,xlRange);     % classic defined kurtosis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures  Plot     pf data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
FF=F*4*1000/n_step; % F scaled

for i=1:20;
figure;plot(FF(1:n_step/2+1),Y_fdk_ave(1:n_step/2+1,i));
% figure;plot(Y_fdk_ave(1:n_step/2,i));
% set(gca,'xlim', [1, 10]);
title(['Frequency Domain Kurtosis in ',num2str(t_tot),' ms window @ ',num2str(v(i)),'v']);xlabel('Frequency (KHz)');ylabel('Frequency Domian Kurtosis K(w)');
%     figure;errorbar(F(1:n_step/2),Y_fdk_ave(1:n_step/2,i),std_Y_fdk_ave(1:n_step/2,i));
%     xlabel('KHz');title('Phase kurtosis in 61ms window');ylabel('Frequency Domian Kurtosis K(w)');
saveas(gcf,['C:\Documents and Settings\qing\Desktop\kurtosis\phase kurtosis plot\v_',num2str(v(i)),'.emf']);
end
