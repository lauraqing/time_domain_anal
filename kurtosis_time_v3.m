%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do acoustic wave data analysis in time domain: PF data               %
% Qing Wu 2013-11-12 version 001                                          %
% Key parameters:Kurtosis in time domain to calculate the same M          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;
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
    T=(1/Fs)*t_tot;      % For energy caculation, test time
    t1=0;                % For LAeq calculation
    % t2=10;             % Already tested the time won't bother LAeq too much! Q.W.
    t2=60;               % consider the 8 hour for LAeq_8h, use 10min
    con=(1/(t2-t1))*(1/(p0^2));  % constant 4.1667e7
    con1=10*log10(600/(8*3600)); % consider the 8 hour for LAeq_8h, use 10min
    v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];
    M = [0.4 0.5 0.6 0.8 1.0 3.0]; %ms
    % use execl to draw durations 
    filename = 'kurtosis_time_domain_M_1.2.xlsx';
    sheet1=1; % classic defined kurtosis
    sheet2=2; % time cycle for classic defined kurtosis in time domain
    % acoustic defined kurtosis
    xlRange='B2';
    
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
    L=length(x_imp); 
    dB(j)=20*log10(max(x_imp(:,j)/(p0)));
 
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
    temp_tru_peak(j)=max(x(:,j))*1.3; % Increase 2dB:1.585     
    %% Kurtosis calculation in time phase   
    % classic kurtosis  %defined a kurtosis wave
    min_x_kur=n0_1-1;
    max_x_kur=n0_1+round(1.2*L/1000)+1;
    x_kur(:,j)=x((min_x_kur:max_x_kur),j);
    N2(j)=length(x_kur);
    ave_x_kur(j)=sum(x_kur(:,j))/N2(j);    
    dif_x_kur=zeros(N2(j));
    dif_x_kur=(x_kur(:,j)-ave_x_kur(j));
    % time kurtosis for a unique M
     kurtosis_x_imp(j)=((sum((dif_x_kur).^4))*(N2(j)-1))/((sum(dif_x_kur.^2)).^2);    
    end
    % Classic kurtosis
    ave_kur_x_imp(i)=sum(kurtosis_x_imp)/10;
    std_ave_kur_x_imp(i)=std(kurtosis_x_imp);   
end

% write data into excel table only draw average peak vs voltage
% classic defined kurtosis
x_wrt_k=[v',ave_kur_x_imp',std_ave_kur_x_imp']; 
xlswrite(filename,x_wrt_k,sheet1,xlRange);     % classic defined kurtosis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures  Plot     pf data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% figure;stem(v,ave_kur_x_imp);title('Classic kurtosis vs. Voltage');xlabel('voltage (v)');ylabel('Average Classic Kurtosis');
figure;errorbar(v,ave_kur_x_imp,std_ave_kur_x_imp);title('Classic kurtosis vs. Voltage');xlabel('voltage (v)');ylabel('Average Classic Kurtosis');
% figure;stem(ave_kur_x_imp,kur_N);title('Classic kurtosis vs. Time cycle');xlabel('Average Classic Kurtosis');ylabel('Time cycle (ms)');
% figure;stem(v,kur_N);title('Time cycle vs. Voltage');xlabel('voltage (v)');ylabel('Time cycle (ms)');

