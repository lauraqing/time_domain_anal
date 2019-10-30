% 10*20 data pf to normalize
% Qing Wu 
% Plot normalization pic

clc;
close all;
clear all;
    
p0=2e-5;             % referrence pressure p0
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
    %% Normalization for impulse noise
    M=length(x_imp); 
    dB(j)=20*log10(max(x_imp(:,j)/(p0)));
    num_norm=30+i;   % Select window for normalization as number from 15 to 300
    cons_ave_norm(j)=sum(x_imp((1:num_norm),j))/num_norm;
    x_imp_norm(:,j)=x_imp(:,j)-cons_ave_norm(j);

    end
    cons_norm(:,i)=cons_ave_norm;
    max_cons_norm(i)=find(cons_ave_norm==max(cons_ave_norm));
    
    % Plot comparation to normalization
    figure;stem(cons_norm(:,i));title(['constant voltage =',num2str(v(i)),'(v)'])
%     figure;subplot(211); plot(x_imp((1:2000),max_cons_norm(i)));title(['Orignal Impulse noise(voltage =',num2str(v(i)),'v)']);xlabel('Time (ms)');ylabel('Amplitude');
%     subplot(212); plot(x_imp_norm((1:2000),max_cons_norm(i)));title(['Normalized Impulse noise (voltage =',num2str(v(i)),'v)']);xlabel('Time (ms)');ylabel('Amplitude');
%     
end  