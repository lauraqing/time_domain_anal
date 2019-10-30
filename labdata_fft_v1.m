% Qing Wu FFT test on labdata_pf
% select the process to do FFT for acoustic wave
% v1: 2013-11-5 last revised
close all;
clear all;
clc;
% v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];
% v(1)=[1  2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17 18  19  20]
% v=[0.3 0.8 1.2 2.0 3.0 5.0];
v=[0.3 0.8 1.2 2.0 3.0 6.0];
A = 2000;
p0=2e-5;
    %  130 dB Peak SPL @ 0.3 v
    x1=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\1\testdata_003.lvm');  
    % 139.8 dB Peak SPL @ 0.8v
    x2=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\3\testdata_007.lvm');  
    % 143.45 dB Peak SPL @ 1.2v
    x3=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\5\testdata_007.lvm');  
    % 149 dB Peak SPL @ 2v
    x4=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\8\testdata_009.lvm');  
%     % % 151dB @ 3v
%     x5=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\10\testdata_005.lvm'); 
    % % 155 dB Peak SPL @ 3v
    x5=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\12\testdata_008.lvm'); 
%     % 155dB @ 5v
%     x6=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\14\testdata_005.lvm');  
    % 158.8 dB Peak SPL @ 6v
    x6=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\16\testdata_010.lvm');  

    data=[x1(:,2) x2(:,2) x3(:,2) x4(:,2) x5(:,2) x6(:,2)];
    
for i=1:6
    y=data(:,i);   
%     y_dB=0;
    
    y_dB(i)=20*log10(max(y)/p0)
    L = length(y);

    y_fft = fft(y, L);               %2
    Y =abs(y_fft(1:L));
    Y_temp=Y/(L/2);
    Y_dB =10*log10(Y_temp(1:L/2).^2); %3

    y_dB_SPL=20*log10(y/p0);
    y_dB_SPL_fft = fft(y_dB_SPL, L);
    Y_dB_SPL_fft = abs(y_dB_SPL_fft(1:L))/(L/2);
    Y_dB_SPL_fft_dB = 10*log10(Y_dB_SPL_fft(1:L/2).^2/406);   %4   
  
    % display a 10ms original impulse wave
    step=600; % 65536->1s; ??->10ms??: answer is 600! Q.W.
    y_max=find(y==max(y));
    shift_cons=120;
    y_start=y_max-shift_cons;
    y_end=y_max+step-shift_cons;
    y_len=length(y(y_start:y_end));
    t=0:1/y_len:1;
    t=t*10;
%     figure;plot(t(1:601),y(y_start:y_end));   
%     set(gca,'ylim', [-400, 1600]);
%     set(gca,'FontSize',20); 
%     saveas(gcf,['C:\Documents and Settings\qing\Desktop\FFT_pf\temp\Imp',num2str(i),'.emf']);

    % refer energy calculation
    % select 10 ms
%     y_ref=y(y_start:y_end); % select 10 ms
%     L_ref=y_len;
%     y_fft_ref= fft(y_ref, L_ref);   
%     Y_ref=abs(y_fft_ref(1:L_ref));
%     Y_ref_temp=Y_ref/(L_ref/2);
%     ref(i)=sum(abs(Y_ref_temp.^2))/ 406  %ms->s should take 1e-3 to multiply
    % select 33 ms
    y_ref=data((200:2362),i); % select 33 ms
    L_ref=length(data((200:2363),i));
    y_fft_ref= fft(y_ref, L_ref);   
    Y_ref=abs(y_fft_ref(1:L_ref));
    Y_ref_temp=Y_ref/(L_ref/2);
    ref(i)=sum(abs(Y_ref_temp.^2.*(33*10e-3)))/ 406  %ms->s should take 1e-3 to multiply
%     ref_nor(i)=ref(i)
    
% %     figure; semilogx(Y_dB);title(['FFT Spectrum (',num2str(v(i)),' v)']);
%     figure; semilogx(Y_dB);title('FFT Spectrum ');
%     set(gca,'xlim', [10, 30000],'ylim', [-150, 10]);
%     ylabel('Power spectrum (dB)');xlabel('Frequency (Hz)');
%     saveas(gcf,['C:\Documents and Settings\qing\Desktop\FFT_pf\temp\FFT',num2str(i),'.emf']);
end
% figure; stem(y_dB,ref);xlabel('dB SPL');ylabel('ref (*10e-3)');title('ref vs peak dB SPL');

