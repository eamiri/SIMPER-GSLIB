% close all
clc
clear all
nTimestep=15*2*24*365+1;
fc1=0.002; %low-pass frequency
%% Inputs
subplot(3,1,1)
% BCtop=readtable('PlateauBC_Intact.csv');
tempData=readtable('PLT2_SHLW.csv','header',1);
tempData=removevars(tempData,'Var3');
varNames = {'timestamp','time'};
for i=3:size(tempData,2)
    varNames{i} = strcat('temp', num2str(i-2));
end
tempData.Properties.VariableNames = varNames;
BCtop = table('Size', [size(tempData,1),2], 'VariableTypes', {'double','double'}, 'VariableNames', {'time','temperature'});
smoothTemp=zeros(nTimestep,size(tempData,2)-1);
smoothTemp(:,1)=(0:0.5:(nTimestep-1)*0.5)';
for n=3:size(tempData,2) 
    BCtop.time=tempData.time;
    BCtop.temperature=tempData{:,n};
    BCcor=[BCtop.time(1),BCtop.temperature(1)];
    for i=2:(length(BCtop.temperature))
        dt=BCtop.time(i) - BCtop.time(i-1);
        if dt==1
            BCcor=[BCcor;[0.5*(BCtop.time(i) + BCtop.time(i-1)),0.5*(BCtop.temperature(i) + BCtop.temperature(i-1));BCtop.time(i),BCtop.temperature(i)]];
        else
            break;
        end
    end
    BCcor = [BCcor;[BCtop.time(i:end),BCtop.temperature(i:end)]];
    for i=1:length(BCcor(:,1))
        if (BCcor(i,2)<-12)
            try
                BCcor(i,2)=BCcor(i+2*24*365,2);
            catch
                BCcor(i,2)=BCcor(i-2*24*365,2);
            end
        end
    end
    % BCtop2=zeros(2*24*365*15,2);
    % for i=1:length(BCtop2(:,2))
    % %     ibc=ceil(i/(24));% FOR DAILY TEMPERATURE DATA
    %     BCtop2(i,:)=[i,BCcor(i,2)];
    % end
    % BCtop=BCtop2;
%     BCtop=BCcor(1:2*24*365*15,:);
    
    fs =1;  % Specify sampling rate - how many datapoints per second
    y = BCcor(1:nTimestep,2);   % Corresponding vector that you want to calculate fft for and filter
    
    
    % If you only want to filter a subsection of
    % data,uncomment these rows below
    % t_start = 25*fs;
    % t_end = 35*fs;
    % y = y(t_start:t_end)
    
    % figure
    subplot(3,1,1)
    plot(y)
    hold on;
    title('Raw Data')
    L = length(y);
    w = hann(L);
    y2 = y.*w;
    Y = fft(y2);
    P2 = abs(Y/L);
    P1 = P2(1:(floor(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(L/2))/L;
    subplot(3,1,2)
    plot(f,P1)
    hold on;
%     fc1 = input('What frequency would you like to low-pass your data at?');
    fnorm=fc1/(fs/2);
    [B,A] = butter(6,fnorm,'low');
    y= filtfilt(B,A,y);
    subplot(3,1,3)
    plot(y)
    hold on;
    smoothTemp(:,n-1) = y;
end



