% Read peak and true peak_pf_v0.20131022
% Qing Wu 2013-10-24 
% Adding xls. wrting function, success

clc;
clear all;
close all;

v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0];
% v = [0.3 0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.5 3.0 3.5 4.0 4.5 5.0]; % decide use 14 data
% Output data into .xls
filename = 'peak_dB_v_pf_20data.xlsx';
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
          x_temp(:,j)=x(:,2,j);
          peak(j)=max(x_temp(:,j));      

   end
     ave_peak(i)=sum(peak)/10;
     dB_peak=20*log10(peak/(2e-5));
     dB_ave_peak(i)=20*log10(ave_peak(i)/(2e-5));
     std_peak_dB(i)=((sum((dB_peak-dB_ave_peak(i)).^2))/10)^.5;  % use std func2: n
     std_peak(i)=((sum((peak-ave_peak(i)).^2))/10)^.5;  % use std func2: n

%      std_peak(i)=std(peak);    % use std func1: n-1
     max_peak(i)=max(peak);

     %%draw each peak
%      figure; stem(peak);title(['(',num2str(i),')v',num2str(v(i))]);
     %%Wave form draw
%     figure; plot(x(:,10));title(['(',num2str(i),')v',num2str(v(i))]);
%        xlRange='A1:J1';
%      %%use to put all data into an excel table
%      xlRange='A1';
%      xlswrite(filename,x_temp,sheet,xlRange);
%      sheet=sheet+1; 
end

% write data into excel table only draw average peak vs voltage
x_wrt_dB=[v',dB_ave_peak',std_peak_dB'];
x_wrt=[v',ave_peak',std_peak'];
header_dB=['','V','P','s'];
header=[' ','V','P','s'];
% xlswrite(filename,header,x_wrt,sheet,xlRange);
xlswrite(filename,x_wrt_dB,sheet1,xlRange);
xlswrite(filename,x_wrt,sheet2,xlRange);

figure; stem(v,dB_ave_peak);title('peak vs. voltage');ylabel('peak (dB)');xlabel('voltage (v)');
% figure; stem(v,max_peak);title('max peak vs. voltage');ylabel('peak (Pa)');xlabel('voltage (v)');

% e = 0.1 * std(ave_peak);
% errorbar(v,ave_peak,e);title('peak vs. voltage');ylabel('peak (Pa)');xlabel('voltage (v)');

% % test bench
%  x=load('C:\Documents and Settings\qing\My Documents\LAB_QING\lab data backup\test result_pengfei\14\test_010.lvm');
% figure;plot(x(:,2));
% 
% 
%     
%     filepath=('D:\test result_pengfei\1.2\');
%     peak=0;
%    for j=1:10 
%        
%        if(j~=10)
%           x(:,:,j)=load([filepath,'test_00',num2str(j),'.lvm']); %test_00i.lvm
%        else
%           x(:,:,j)=load([filepath,'test_0',num2str(j),'.lvm']); %test_00i.lvm 
%        end
%           x_temp(:,:,j)=x(:,:,j);
%           peak(j)=max(x_temp(:,2,j));              
%    end
%    
% figure; stem(peak);title('(3)v1.2');title(['(',num2str(i),')v',num2str(v(i))]);
