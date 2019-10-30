% Plot one cycle waveform by Qing Wu 2013-11-6

clc;
clear all;
close all;
v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];   
%     v = [1    2   3   4   5   6   7  8   9   10   11  12  13  14  15  16  17  18  19  20];  

for i=1:20
file=(['C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\',num2str(i),'\testdata_004.lvm']);

Fs=65536;
t=0:1/Fs:1;
x=load(file);
x_imp=x(:,2);

dB=20*log10(max(x_imp)/(2e-5))

figure;
plot(t(1:501),x_imp(600:1100),'linewidth',2); set(gca,'FontSize',20);  % Dr Qin wants font 4
xlabel('TIME IN SECONDS');ylabel('PRESSURE IN PASCALS ');title(['wavefrom in voltage ',num2str(v(i)),'v'])
grid on;
saveas(gcf,['C:\Documents and Settings\qing\Desktop\waveform_',num2str(v(i)),'.emf']);

end
