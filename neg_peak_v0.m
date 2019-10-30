% Read negative peak
% Qing Wu 2013-11-12 
% Adding xls. wrting function, success

clc;
clear all;
close all;
% constant definition
    p0=2e-5;             % referrence pressure p0
    dis1=2e-5;           % Microphone 20usec
    dis2=5e-1;           % Microphone 20usec
    step=6;              % To find the A-duration length by intersection
    Fs=44100;            % Sampling Frequency of energy caculation
    t_tot=33.11;         % ave_A_dur=sum(a_dur)/10=9.46ms; t+/t-=1/2.5; t_tot=33.11ms
  v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];
% Output data into .xls
filename = 'negative_peak_dB.xlsx';
sheet1=1;sheet2=2;xlRange='B2';

for i=1:20
    %%data from 20131022
    filepath=(['C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\',num2str(i),'\']);
    peak=0;
    
   for j=1:10     
       if(j~=10)
          x(:,:,j)=load([filepath,'testdata_00',num2str(j),'.lvm']); %testdata_00i.lvm
       else
          x(:,:,j)=load([filepath,'testdata_0',num2str(j),'.lvm']); %testdata_00i.lvm 
       end
          x_imp(:,j)=x(:,2,j);    
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
    
    
    min_x_kur=n0_0+1;
    max_x_kur=2000;
    x_tem=zeros(length(x((min_x_kur:max_x_kur),j)),10);
    x_tem(:,j)=x((min_x_kur:max_x_kur),j);
    x_tem(:,j)=abs(x_tem(:,j));
   
    neg_peak(j)=max(x_tem(:,j));    
   end
     ave_peak(i)=sum(neg_peak)/10;
     dB_peak=20*log10(neg_peak/(2e-5));
     dB_ave_peak(i)=20*log10(ave_peak(i)/(2e-5));
     std_peak_dB(i)=((sum((dB_peak-dB_ave_peak(i)).^2))/10)^.5;  % use std func2: n
     std_peak(i)=((sum((peak-ave_peak(i)).^2))/10)^.5;  % use std func2: n

     max_peak(i)=max(peak);
end

% write data into excel table only draw average peak vs voltage
x_wrt_dB=[v',dB_ave_peak',std_peak_dB'];
% header_dB=['','V','P','s'];
% header=[' ','V','P','s'];
xlswrite(filename,x_wrt_dB,sheet1,xlRange);

x_wrt=[v',ave_peak',std_peak'];
xlswrite(filename,x_wrt,sheet2,xlRange);

figure; stem(v,dB_ave_peak);title('negative peak vs. voltage');ylabel('peak (dB)');xlabel('voltage (v)');
% figure; stem(v,max_peak);title('max peak vs. voltage');ylabel('peak (Pa)');xlabel('voltage (v)');
