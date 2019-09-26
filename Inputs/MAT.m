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
xverts = [BCtop(1:end-1,1)'; BCtop(1:end-1,1)'; BCtop(2:end,1)'; BCtop(2:end,1)'];
yverts = [zeros(1,length(BCtop(:,1))-1); BCtop(1:end-1,2)'; BCtop(2:end,2)'; zeros(1,length(BCtop(:,1))-1)];
p = patch(xverts,yverts,'b','LineWidth',1.5);

corBCtop=BCtop;
% corBCtop(:,2)=BCtop(:,2)-2.0943; %ZERO
%corBCtop(:,2)=BCtop(:,2)+0.9058; %Positive 3
% corBCtop(:,2)=BCtop(:,2)-5.0943; %Negative 3
% corBCtop(:,2)=BCtop(:,2) * 2.5;

meanAnnTemp = trapz(corBCtop(:,1),corBCtop(:,2))/(max(corBCtop(:,1)))



