%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do acoustic wave data analysis in time domain: PF data               %
% Qing Wu 2013-10-27                                                      %
% Key parameters analysis: E/Laeq                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
%% Load impulse noise waveform %
% constant definition
    p0=2e-5;             % referrence pressure p0
    dis1=2e-5;           % Microphone 20usec
    dis2=5e-1;           % Microphone 20usec
    step=6;              % To find the A-duration length by intersection
    Fs=44100;            % Sampling Frequency of energy caculation
    t_tot=33.11;         % ave_A_dur=sum(a_dur)/10=9.46ms; t+/t-=1/2.5; t_tot=33.11ms
    T=(1/Fs)*t_tot;      % For energy caculation, test time
    t1=0;                % For LAeq calculation
    % t2=10;             % Already tested the time won't bother LAeq too much! Q.W.
    t2=60;               % consider the 8 hour for LAeq_8h, use 10min
    con=(1/(t2-t1))*(1/(p0^2));  % constant 4.1667e7
    con1=10*log10(600/(8*3600)); % consider the 8 hour for LAeq_8h, use 10min
    v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];

for i=1:20
    filepath=(['C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\',num2str(i),'\']);
    %% Ininitial the temp variances
    peak=0;n1=0;n2=0;x_0=0;abs_x_0=0;n0=0;n0_0=0;a=0;k=0;x_01=0;temp_neg=0;a_dur1=zeros(10,1);
    n00=0;a_dur2=zeros(10,1);n0_1=0;n_half=0;M=0;temp_com_dB=0;temp_peak=zeros(1,10);
    %     x_kur=zeros(:,10); tem_x_kur=zeros(:,10);% Define kurtosis calculation in small regions 
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
%     %% Normalization for impulse noise
    M=length(x_imp); 
    dB(j)=20*log10(max(x_imp(:,j)/(p0)));
%     num_norm=i*25;   % Select window for normalization as number from 25 to 500
%     cons_ave_norm(j)=sum(x_imp((max(x_imp(:,j))-num_norm:max(x_imp(:,j))+num_norm),j))/num_norm;
%     x_imp_norm(:,j)=x_imp(:,j)-cons_ave_norm(j);
% %     figure;subplot(211); plot(x_imp(1:2000));title(['Orignal Impulse noise (peak =',num2str(dB),'dB)']);xlabel('Time (ms)');ylabel('Amplitude');
% %     subplot(212); plot(x_imp_norm(1:2000));title(['Normalized Impulse noise (peak =',num2str(dB),'dB)']);xlabel('Time (ms)');ylabel('Amplitude');
%     x(:,j)=x_imp_norm(:,j);   
    x(:,j)=x_imp(:,j);  
    %% Find the true peak pressure and A-duration lengths: use extropolation method
    n1=find(x(:,j)==max(x(:,j))); % find the peak postion
    n1=max(n1);
    n2=find(x(:,j)==min(x(:,j))); 
    x_0=x((n1:n2),j); 
    abs_x_0=abs(x_0);
    n0=find(abs_x_0==min(abs_x_0));% find the A-duration time for impulse noise x
    n0_0=n1+n0;                    % find the first zero point position after peak
    a_dur1(j)=n0_0-n1;                % length of A-duration afer peak: A_dur1 
    % spline interpolation to find x_c
    a=n1+dis1;
    x_c=interp1((n1:n2),x((n1:n2),j),a,'spline');
    k=(max(x(:,j))-x_c)/dis1;              % find the slope in scale 2us
    % A-duration length
    x_01=x(n1-step:n1,j);                 % assume the peak form like sinewave
    abs_x_01=abs(x_01);
    temp_neg=find(abs_x_01==min(abs_x_01));   % first negative point
    x_1=x(n1-temp_neg:n1,j);            % find the linear increasing sequence
    abs_x_1=abs(x_1);

    n00=find(abs_x_1==min(abs_x_1));% find the A-duration time for impulse noise x
    n00=min(n00);  
    a_dur2(j)=length(x_1)-n00;           % find the first zero point position before peak
    n0_1=n1-a_dur2(j);
    n_half=round(n1-a_dur2(j)/2);
%     temp_tru_peak(j)=max(x(:,j))*1.3; % Increase 2dB:1.585
    if max(x(:,j))/2>abs(x(n_half,j)-max(x(:,j)))>k*a_dur2(j)/2
        b=n1-dis2;
        x_half=interp1((n1-temp_neg:n1),x_1,b,'spline');
        temp_tru_peak(j)=x_half+k*a_dur2(j)/2;   
    end
    if 0< abs(x(n_half,j)-max(x(:,j)))< k*a_dur2(j)/2
        temp_tru_peak(j)=x(n_half,j)+k*a_dur2(j)/2;
    end 
    if abs(x(n_half,j)-max(x(:,j)))>max(x(:,j))/2
        x_half=max(x(:,j));
        temp_tru_peak(j)=x_half+k*a_dur2(j)/2;
    end 
    if abs(x(n_half,j)-max(x(:,j)))==0
        temp_tru_peak(j)=max(x(:,j))*1.6; % Increase 2dB:1.585
    end

    % Find out the true peak
    dB_tru_peak(j)=20*log10(temp_tru_peak(j)/p0);
    temp_com_dB(j)=dB_tru_peak(j)-dB(j);
    % A-duration length
    temp_a_dur(j)=a_dur1(j)+a_dur2(j);
    % negative-duration length
    temp_neg_dur(j)=2*max(n2-n0_0);
    
    %% Calculation of Energy and LAeq (in 8 hour)
    s(j)=sum(x(:,j).^2)*(1/M); %M=65536
    temp_E(j)=s(j)*T;              
    % compute the LAeq
    % s=sum(gy3(1:512000).^2)*(1/Fs);
    temp_LAeq(j)=10*log10(con*s(j));
    temp_LAeq_8h(j)=temp_LAeq(j)+con1;       
    %% Kurtosis calculation in time phase   
    % classic kurtosis
    min_x_kur=n0_1-1;
    tem_max_x_kur=n0_1+length(x_0)+1;
    tem_x_kur=zeros(tem_max_x_kur-min_x_kur+1,10);
    tem_x_kur(:,j)=x((min_x_kur:tem_max_x_kur),j);
    N1=length(tem_x_kur(:,j));
    max_x_kur=min_x_kur+2*N1;%2*length(x_0); % Test cycle for kurtosis, unit in 'ms'
    x_kur=x((min_x_kur:max_x_kur),j);
    N2=length(x_kur);
    ave_x_kur(j)=sum(x_kur)/N2;
    dif_x_kur=(x_kur-ave_x_kur(j));
    kurtosis_x_imp(j)=((sum((dif_x_kur).^4))*(N2-1))/((sum(dif_x_kur.^2)).^2);    
    % acoustic kurtosis

    % N_sec=100; 
% dif_x_imp_sec=zeros(1,N_sec);
% for i=1:M_imp/N_sec
%    
%     min_x_imp_sec=1+(i-1)*N_sec;
%     max_x_imp_sec=N_sec+(i-1)*N_sec;
%    
%     ave_x_imp_sec=sum(x_imp(min_x_imp_sec:max_x_imp_sec))/N_sec;
%     dif_x_imp_sec(:)=(x_imp(min_x_imp_sec:max_x_imp_sec)-ave_x_imp_sec); 
%     sd_x_imp_sec=std(x_imp(min_x_imp_sec:max_x_imp_sec));
%    
%     kurtosis_x_imp_sec(i)=((sum((dif_x_imp_sec).^4))/(N_sec-1))/(sd_x_imp_sec)^4;    
%     max_kurtosis_x_imp_sec=max(kurtosis_x_imp_sec);
% end  
%     
% % subplot(2,1,2);
% figure;stem(kurtosis_x_imp_sec);title(['Kurtosis of Impulse Noise (N=',num2str(N_sec),') max kurtosis=', num2str(max_kurtosis_x_imp_sec)]);
    end
    % Peak & compensation True Peak
    ave_peak(i)=sum(temp_peak)/10;
    com_dB(i)=sum(temp_com_dB)/10;
    tru_peak(i)=sum(temp_tru_peak)/10;
    dB_peak=20*log10(peak/(2e-5));
    dB_ave_peak(i)=20*log10(ave_peak(i)/(2e-5));
    std_peak_dB(i)=((sum((dB_peak-dB_ave_peak(i)).^2))/10)^.5;  % use std func2: n
    std_peak(i)=((sum((peak-ave_peak(i)).^2))/10)^.5;  % use std func2: n
    % Durations
    a_dur(i)=sum(temp_a_dur)/10;
    neg_dur(i)=sum(temp_neg_dur)/10;
    % Energy
    E(i)=sum(temp_E)/10;
    LAeq(i)=sum(temp_LAeq)/10;
    LAeq_8h(i)=sum(temp_LAeq_8h)/10;
    % Classic kurtosis
    ave_kur_x_imp(i)=sum(kurtosis_x_imp)/10;
    kur_N(i)=sum(N2)/10;
    % Normalization
    cons_norm(:,i)=cons_ave_norm;
    max_cons_norm(i)=find(cons_ave_norm==max(cons_ave_norm));
%     % Plot comparation to normalization
%     figure;plot(cons_norm(:,i));title(['constant voltage =',num2str(v(i)),'(v)'])
%     figure;subplot(211); plot(x_imp((1:2000),max_cons_norm(i)));title(['Orignal Impulse noise(voltage =',num2str(v(i)),'v)']);xlabel('Time (ms)');ylabel('Amplitude');
%     subplot(212); plot(x_imp_norm((1:2000),max_cons_norm(i)));title(['Normalized Impulse noise (voltage =',num2str(v(i)),'v)']);xlabel('Time (ms)');ylabel('Amplitude');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures  Plot     Ad data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% figure; stem(v,ave_peak);title('Impulse noise Amplitude vs. Voltage');xlabel('voltage (v)');ylabel('Amplitude (Kpa)');
% figure; stem(v,tru_peak);title('Impulse noise true peak presure vs. Voltage');xlabel('voltage (v)');ylabel('True peak (Kpa)');
% figure;stem(v,a_dur);title('A-duration vs. Voltage');xlabel('voltage (v)');ylabel('A-duraton (s)');
% figure;stem(v,neg_dur);title('Negative-duration vs. Voltage');xlabel('voltage (v)');ylabel('Negative-duraton (s)');
% figure;stem(v,com_dB);title('True peak compensation in dB vs. Voltage');xlabel('voltage (v)');ylabel('True peak compensation (dB)');
% figure;stem(v,E);title('Energy vs. Voltage');xlabel('voltage (v)');ylabel('Energy (Pa*Pa*s)');
figure;stairs(v,LAeq);title('LAeq vs. Voltage');xlabel('voltage (v)');ylabel('LAeq (Pa*Pa*s)');
figure;stairs(v,LAeq_8h);title('LAeq8h in dB vs. Voltage');xlabel('voltage (v)');ylabel('LAeq_8h (Pa*Pa*s)');
% figure;stem(v,ave_kur_x_imp);title('Classic kurtosis vs. Voltage');xlabel('voltage (v)');ylabel('Average Classic Kurtosis');
% figure;stem(ave_kur_x_imp,kur_N);title('Classic kurtosis vs. Time cycle');xlabel('Average Classic Kurtosis');ylabel('Time cycle (ms)');
% figure;stem(v,kur_N);title('Time cycle vs. Voltage');xlabel('voltage (v)');ylabel('Time cycle (ms)');
