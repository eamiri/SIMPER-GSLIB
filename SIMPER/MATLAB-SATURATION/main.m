clc
clear

Sres = 0.1;
Rsat = 2;
Wpar = 1;
Mpar = 1;

Tsol = -1;
Tliq = 0;
Tgau = (-1.5:0.001:0.5)';
Sw = zeros(length(Tgau), 1);
ISw = zeros(length(Tgau), 1);
ISwxT = zeros(length(Tgau), 1);
IdSwxT = zeros(length(Tgau), 1);

Si = zeros(length(Tgau), 1);
ISi = zeros(length(Tgau), 1);
ISixT = zeros(length(Tgau), 1);
IdSixT = zeros(length(Tgau), 1);

for i=1:length(Tgau)    
    Sw(i,1) = (1 - SATUR(Tgau(i,1), Tsol, Tliq, 0, Rsat, Wpar, Mpar)) * SATUR(Tgau(i,1), Tsol, Tliq, 0, Rsat, Wpar, Mpar) + Sres;
    ISw(i,1) = ISATUR(Tgau(i,1), Tsol, Tliq, Sres, Rsat, Wpar, Mpar);
    ISwxT(i,1) = ISATURxT(Tgau(i,1), Tsol, Tliq, Sres, Rsat, Wpar, Mpar);
    IdSwxT(i,1) = IdSATURxT(Tgau(i,1), Tsol, Tliq, Sres, Rsat, Wpar, Mpar);
    
    Si(i,1) = (1 - SATUR(Tgau(i,1), Tsol, Tliq, Sres, Rsat, Wpar, Mpar)) * (1 - SATUR(Tgau(i,1), Tsol, Tliq, 0, Rsat, Wpar, Mpar));
%     ISi(i,1) = Tgau(i,1) - ISw(i,1);
    ISi(i,1) =1 - Si(i,1) - Sw(i,1);
end



subplot(2,2,1); 
plot(Tgau,Sw);

subplot(2,2,2); 
plot(Tgau,Si);

subplot(2,2,3); 
plot(Tgau,ISw);

subplot(2,2,4); 
plot(Tgau,IdSwxT);
