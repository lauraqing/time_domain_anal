%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do acoustic wave data analysis in time domain: PF data               %
% Qing Wu 2013-11-3    version 2                                          %
% Key parameters:Duration figures                                         %
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
%   v = [1    2   3   4   5   6   7  8   9   10   11  12  13  14  15  16  17  18  19  20];  

    % use execl to draw durations 
    filename = 'durations_v3.xlsx';
    sheet1=1;sheet2=2;xlRange='B2';    

for i=1:20
    filepath=(['C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\',num2str(i),'\']);
%     filepath=('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\11\');
%     %% Ininitial the temp variances
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
    
    % A-duration length
    temp_a_dur(j)=a_dur1(j)+a_dur2(j);
  if (i<=8)  
      % negative-duration length
      temp_neg_dur(j)=2*max(n2-n0_0);
  else
    %find negative duration:half_neg_dur defined from negative peak to 0
    neg_dur_con=[0;0;0;0;0;0;0;0;5;5;5;6;6;6;7;7;7;8;8;8];
    neg_x_left=x((n0_0:n0_0+neg_dur_con(i)*ceil(temp_a_dur(j))),j);
    abs_neg_x_left=abs(neg_x_left);
    n0_2=find(abs_neg_x_left==min(abs_neg_x_left));
    n0_2=n0_0+n0_2;
    n0_2_temp=[abs(x(n0_2-1)); abs(x(n0_2)); abs(x(n0_2+1)) ];
    n0_2_idx=find(n0_2_temp==min(n0_2_temp));
    n0_2_cons=[-1; 0 ; 1];
    n0_2=max(n0_2+n0_2_cons(n0_2_idx));
    %Debuge okay: QW 2013-11-13
    value_n0_2_1=x(n0_2-1);
    value_n0_2=x(n0_2);
    value_n0_2_2=x(n0_2+1);
%     figure;plot(x);title(['impulse noise @',num2str(v(i)),' v']);
    half_neg_dur=max(n0_2-n0_0);
    
    % negative-duration length
    temp_neg_dur(j)=max(n2-n0_0)+half_neg_dur+1;
    end
    %duration transfer into time scale
    temp_a_dur(j)=(a_dur1(j)+a_dur2(j))/65.536;
    temp_neg_dur(j)=temp_neg_dur(j)/65.536;
  
    end
    % Durations
    a_dur(i)=sum(temp_a_dur)/10;
    std_a_dur(i)=std(temp_a_dur);
    neg_dur(i)=sum(temp_neg_dur)/10;
    std_neg_dur(i)=std(neg_dur);

end
% write data into excel table only draw average peak vs voltage
x_wrt_dB=[v',a_dur',std_a_dur'];
x_wrt=[v',neg_dur',std_neg_dur'];
header_dB=['','V','A','s'];
header=[' ','V','N','s'];
% xlswrite(filename,header,x_wrt,sheet,xlRange);
xlswrite(filename,x_wrt_dB,sheet1,xlRange);
xlswrite(filename,x_wrt,sheet2,xlRange);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures  Plot     pf data   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% figure;stem(v,a_dur);title('A-duration vs. Voltage');xlabel('voltage (v)');ylabel('A-duraton (ms)');
% figure;stem(v,neg_dur);title('Negative-duration vs. Voltage');xlabel('voltage (v)');ylabel('Negative-duraton (ms)');
figure;errorbar(v,a_dur,std_a_dur);title('A-duration vs. Voltage');xlabel('voltage (v)');ylabel('A-duraton (ms)');
figure;errorbar(v,neg_dur,std_neg_dur);title('Negative-duration vs. Voltage');xlabel('voltage (v)');ylabel('Negative-duraton (ms)');
% figure;stem(v,com_dB);title('True peak compensation in dB vs. Voltage');xlabel('voltage (v)');ylabel('True peak compensation (dB)');