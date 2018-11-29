clc
clear

syms Tgau

syms Sres Tliq Tsol

% Sres = 0.1;
% Rsat = 2;
% Wpar = 1;
% Mpar = 1;

npor = 0.8;
Dwat = 1000;
Dice = 920;
Cice = 2108;
Cwat = 4182;
Kwat = 0.6;
Kice = 2.14;
Lhea = 334000;
Csol = 630;
Ksol = 0.05;
Dsol = 250;

Tsol = -1;
Tliq = 0;

SFC = 0.5 * (1 + sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol)));
SFC_Sres = 0.5*(1+Sres)+0.5*(1-Sres)*sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol));

HS(Tgau) = (heaviside(Tgau + 1) - heaviside(Tgau));

Swat = simplify(HS(Tgau)*(1-SFC)*(SFC) + Sres)
dSwat = diff(Swat, Tgau)
ISwat = int(Swat, [Tsol Tgau])
IdSwatxT = int(diff(Swat, Tgau) * Tgau, [Tsol Tgau])
ISwatxT = int(Swat * Tgau, [Tsol Tgau])

Sice = (HS(Tgau)) *(1 - SFC_Sres) * (1 - SFC) + (1 - Sres) * heaviside(-Tgau-1)
dSice = diff(Sice, Tgau)
ISice = int(Sice, [Tsol Tgau])
IdSicexT = int(diff(Sice, Tgau) * Tgau, [Tsol Tgau])
ISicexT = int(Sice * Tgau, [Tsol Tgau])


Cpar = npor * (Swat * Dwat * Cwat + Sice * Dice * Cice) + (1.0 - npor) * Dsol * Csol + npor * Dice * Lhea * (-dSice);
ICparxT = npor * (ISwatxT * Dwat * Cwat + ISicexT * Dice * Cice) + (1 - npor) * Dsol * Csol * Tgau * Tgau / 2.0 + npor * Dice * Lhea * (-IdSicexT);
ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice) + (1.0 - npor) * Dsol * Csol * Tgau + npor * Dice * Lhea*(1 - Sice - Sres);
Kpar = (Kwat ^ (Swat * npor)) * (Kice ^ (Sice * npor))* (Ksol ^ (1.0 - npor));

subplot(2,2,1); 
% fplot(Sice, [-1.5, 0.5]);
% 
% subplot(2,2,2); 
% fplot(ISice, [-1.5, 0.5]);
% 
% subplot(2,2,3); 
% fplot(ISicexT, [-1.5, 0.5]);
% 
% subplot(2,2,4); 
% fplot(IdSicexT , [-1.5, 0.5]);