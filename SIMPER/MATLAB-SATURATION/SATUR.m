function Sw=SATUR(Tgau,Tsol,Tliq,Sres,Rsat,Wpar,Mpar)
switch Rsat
    case 1
        if Tgau<=Tsol
            Sw=Sres;
        elseif Tgau<Tliq && Tgau>Tsol
            Sw=(1-Sres)*exp(-((Tgau-Tliq)/Wpar)^2)+Sres;
        else
            Sw=1.0;
        end
    case 0
        Btem=(Sres-1)/Mpar;
        if Tgau>=Btem
            Sw=Mpar*Tgau+1;
            if Tgau>Tliq
                Sw=1.0;
            end
        else
            Sw=Sres;
        end
    case 2
        if Tgau<=Tsol
            Sw=Sres;
        elseif Tgau<Tliq && Tgau>Tsol
            Sw=0.5*(1+Sres)+0.5*(1-Sres)*sin(pi*(Tgau-0.5*(Tliq+Tsol))/(Tliq-Tsol));
        else
            Sw=1.0;
        end
end
end