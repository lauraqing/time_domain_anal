% test frequency
% Qing Wu 2013-11-11

    n1=4000;
    N_max1=65536;
    % Consider the total secend if is 1s
    T_tot1=1;                    % unit is 's'
    t_tot1=T_tot1*n1*1000/N_max1;            % window, unit is 'ms'
    Fs1=N_max1/T_tot1;
    t1=[0:1/Fs:n1*1000/Fs];  
    F1=(Fs1/n1)*([1:n1]-1)/1000;             % calculate the real frequency used for FFT of 4000 points, unit in 1K
     
    step1=10;
    figure;stem(F1(1:n1/2));title(['OUR:frequency value step',num2str(step1),'@selected time ',num2str(t_tot1),'ms']);    
    ylabel('Frequency (KHz)');
    
    n=1024;
    Fs=20000;
    t_tot=51.2;  %ms
    
    % Consider the total secend if is 1s
%     T_tot=1;                    % unit is 's'
    t_tot=1*n*1000/Fs;            % window, unit is 'ms'
    t=[0:1/Fs:n*1000/Fs];  
    F=(Fs/n)*([1:n]-1)/1000;             % calculate the real frequency used for FFT of 4000 points, unit in 1K
     
    step=n;
    figure;stem(F(1:n));title(['OTHER:frequency value step',num2str(step),'@selected time ',num2str(t_tot),'ms']);
    ylabel('Frequency (KHz)');