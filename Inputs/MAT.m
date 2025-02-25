clc
clear
% MEAN ANNUAL TEMPERATURE CALCULATOR
fileid=fopen('PlateauBC_Corrected.dat');
DATA=textscan(fileid,'%f %f', 'headerLines', 3);
BCtop2=[DATA{1}, DATA{2}];
BCtop=zeros(length(BCtop2(:,1)),2);
BCtop(:,1)=(1:length(BCtop2(:,1)))';
BCtop(:,2)=[BCtop2(1:2:end,2);BCtop2(1:2:end,2)];
fclose('all');
% xverts = [BCtop(1:end-1,1)'; BCtop(1:end-1,1)'; BCtop(2:end,1)'; BCtop(2:end,1)'];
% yverts = [zeros(1,length(BCtop(:,1))-1); BCtop(1:end-1,2)'; BCtop(2:end,2)'; zeros(1,length(BCtop(:,1))-1)];
% p = patch(xverts,yverts,'b','LineWidth',1.5);

corBCtop=BCtop;
maxTemp = max(corBCtop(:,2));
minTemp = min(corBCtop(:,2));
coeff=zeros(length(corBCtop(:,1)),2);
for i=1:length(corBCtop(:,1))
    if corBCtop(i,2) > 0
        Temp = corBCtop(i,2);
        coeff(i,1)=i;
        coeff(i,2) = 3.2*exp(-((Temp-maxTemp)/maxTemp)^10);
        corBCtop(i,2) = corBCtop(i,2) / coeff(i,2);
    else
        Temp = corBCtop(i,2);
        coeff(i,1)=i;
        coeff(i,2) = 0.3*exp(-((Temp-minTemp)/minTemp)^10);
        corBCtop(i,2) = corBCtop(i,2) / coeff(i,2);
    end
end

plot(corBCtop(:,1),corBCtop(:,2));
hold on;
plot(BCtop(:,1),BCtop(:,2));
xlim([0 3*365*24]);


% corBCtop(:,2)=BCtop(:,2)-3.95786; %ZERO
% corBCtop(:,2)=BCtop(:,2)-2.95786; %Positive 1
% corBCtop(:,2)=BCtop(:,2)-1.95786; %Positive 2
% corBCtop(:,2)=BCtop(:,2)-0.95786; %Positive 3
% corBCtop(:,2)=BCtop(:,2)-6.95786; %Negative 3
% corBCtop(:,2)=BCtop(:,2) * 2.5;

meanAnnTemp = trapz(corBCtop(:,1),corBCtop(:,2))/(max(corBCtop(:,1)))



