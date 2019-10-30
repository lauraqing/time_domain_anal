% demo program to normalized one program
% Qing Wu 2013-10-24 

clc;
clear all;
close all;
v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];

j=2;dis1=2e-5;
    step=6;              % To find the A-duration length by intersection
    
x=load('C:\Documents and Settings\qing\My Documents\MATLAB\time_domain_anal\test result_pf_20131022\6\testdata_002.lvm');
figure;subplot(211);plot(x(1:2500,2));title('Orignal Signal');xlabel('Time (ms)');ylabel('Amplitude (Pa)');

%% Zero point recognization
    %Find the true peak pressure and A-duration lengths: use extropolation method
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
    n0_1=n1-a_dur2(j);                   % The first zero point
    n_half=round(n1-a_dur2(j)/2);

%% Normalize is shift dc to 0
    cons_loc=n0_1/10;
    cons_norm=sum(x(1:cons_loc,2))/(cons_loc+1)
    x=x-cons_norm;
subplot(212);plot(x(1:2500,2));title('Normalized Signal');xlabel('Time (ms)');ylabel('Amplitude (Pa)');

