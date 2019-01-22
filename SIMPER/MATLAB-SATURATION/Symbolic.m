clc
clear

syms Tgau

% syms Sres Tliq Tsol npor Dwat Dice Cice Cwat Kwat Kice Lhea Csol Ksol Dsol Cair Kair Dair
% 
TempG = -0.5;
Sres = 0.0;
Rsat = 2;
Wpar = 1;
Mpar = 1;
npor = 0.8;
Dair = 1.3;
Cair = 1.0;
Kair = 0.024;
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
Tsol = -1.0;
Tliq = 0.0;

HS(Tgau) = (heaviside(Tgau - Tsol) - heaviside(Tgau - Tliq));
% fplot(HS,[-1.5 1])

% UNSATURATED CONDITION
SFC = 0.5 * (1 + sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol)));
SFC_Sres = 0.5*(1+Sres)+0.5*(1-Sres)*sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol));
Sice = (HS(Tgau)) *(1 - SFC_Sres) * (1 - SFC) + (1 - Sres) * heaviside(-Tgau+Tsol)
Swat = (HS(Tgau)*(1-SFC)*(SFC) + Sres)
Sair = 1 - Swat - Sice

% SATURATED CONDITION
% SFC_Sres = 0.5*(1+Sres)+0.5*(1-Sres)*sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol));
% Sice = simplify(HS(Tgau)*(1.0 - SFC_Sres) + (1-Sres) * heaviside(-Tgau + Tsol))
% Swat = simplify(HS(Tgau)*(SFC_Sres) + 1.0 * heaviside(Tgau - Tliq) + Sres * (heaviside(- Tgau + Tsol)))
% Sair = 1 - Swat - Sice

% WATER
dSwat = simplify(diff(Swat, Tgau))
ISwat = int(Swat, Tgau);
ISwat = simplify(ISwat - subs(ISwat, Tgau, Tsol))

IdSwatxT = int(diff(Swat, Tgau) * Tgau, Tgau);
IdSwatxT = IdSwatxT - subs(IdSwatxT, Tgau, Tsol);

ISwatxT = int(Swat * Tgau, Tgau);
ISwatxT = simplify(ISwatxT - subs(ISwatxT,Tgau,Tsol))

% ICE
dSice = simplify(diff(Sice, Tgau))
ISice = int(Sice, Tgau);
ISice = simplify(ISice - subs(ISice,Tgau,Tsol))

IdSice = int(dSice, Tgau);
IdSice = simplify(IdSice - subs(IdSice, Tgau, Tsol))

IdSicexT = int(diff(Sice, Tgau) * Tgau, Tgau);
IdSicexT = simplify(simplify(IdSicexT - subs(IdSicexT,Tgau,Tsol)))

ISicexT = int(Sice * Tgau, Tgau);
ISicexT = simplify(simplify(ISicexT - subs(ISicexT, Tgau, Tsol)))

% AIR
dSair = simplify(diff(Sair, Tgau))
ISair = int(Sair, Tgau);
ISair = simplify(ISair - subs(ISair,Tgau,Tsol))

IdSairxT = int(diff(Sice, Tgau) * Tgau, Tgau);
IdSairxT = simplify(simplify(IdSairxT - subs(IdSairxT,Tgau,Tsol)))

ISairxT = int(Sice * Tgau, Tgau);
ISairxT = simplify(simplify(ISairxT - subs(ISairxT, Tgau, Tsol)))

% % % % % % % % % % % % % % % % 

Cpar = npor * (Swat * Dwat * Cwat + Sice * Dice * Cice + Sair * Dair * Cair) + (1.0 - npor) * Dsol * Csol + npor * Dice * Lhea * (-dSice);
% ICparxT = int(Cpar * Tgau, Tgau ,[Tsol Tgau]);
% ICpar = int(Cpar, Tgau ,[Tsol Tgau]);
ICparxT = npor * (ISwatxT * Dwat * Cwat + ISicexT * Dice * Cice + ISairxT * Dair * Cair) + (1 - npor) * Dsol * Csol * (Tgau * Tgau - Tsol*Tsol)/ 2.0 + npor * Dice * Lhea * (-IdSicexT);
ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice + ISair * Dair * Cair) + (1.0 - npor) * Dsol * Csol * (Tgau-Tsol) + npor * Dice * Lhea*(1 - Sice - Sres);
Kpar = (Kwat ^ (Swat * npor)) * (Kice ^ (Sice * npor))*(Kair ^ (Sair * npor))* (Ksol ^ (1.0 - npor));

subplot(2,2,1); 
fplot(Sice, [Tsol-0.5, Tliq+0.5]);

subplot(2,2,2); 
fplot(Swat, [Tsol-0.5, Tliq+0.5]);

subplot(2,2,3); 
fplot(Cpar, [Tsol-0.5, Tliq+0.5]);

subplot(2,2,4); 
fplot(Kpar , [Tsol-0.5, Tliq+0.5]);
% 
% 
Swat = eval(subs(Swat, Tgau, TempG))
dSwat = eval(subs(dSwat, Tgau, TempG))
ISwat = eval(subs(ISwat, Tgau, TempG))
IdSwatxT = eval(subs(IdSwatxT, Tgau, TempG))
ISwatxT = eval(subs(ISwatxT, Tgau, TempG))

Sice = eval(subs(Sice, Tgau, TempG))
dSice = eval(subs(dSice, Tgau, TempG))
ISice = eval(subs(ISice, Tgau, TempG))
IdSice = eval(subs(IdSice, Tgau, TempG))
IdSicexT = eval(subs(IdSicexT, Tgau, TempG))
ISicexT = eval(subs(ISicexT, Tgau, TempG))


Cpar = eval(subs(Cpar, Tgau, TempG))
ICparxT = eval(subs(ICparxT, Tgau, TempG))
ICpar = eval(subs(ICpar, Tgau, TempG))
Kpar = eval(subs(Kpar, Tgau, TempG))

Sair = eval(subs(Sair, Tgau, TempG))
ISair = eval(subs(ISair, Tgau, TempG))
ISairxT = eval(subs(ISairxT, Tgau, TempG))


















