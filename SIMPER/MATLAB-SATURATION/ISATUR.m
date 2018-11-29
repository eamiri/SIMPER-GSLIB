function ISw=ISATUR(Tgau,Tsol,Tliq,Sres,Rsat,Wpar,Mpar)
switch Rsat
    case 1 % Exponential FTC
        if Tgau<=Tsol
            ISw=Sres*Tgau;
        elseif Tgau<Tliq && Tgau>Tsol
            CSw1=Sres*Tsol - (Sres*Tsol - (pi^(1/2)*erfi((Tsol - Tliq)*(-1/Wpar^2)^(1/2))*(Sres - 1))/(2*(-1/Wpar^2)^(1/2)));
            ISw=CSw1 + Sres*Tgau - (pi^(1/2)*erfi((Tgau - Tliq)*(-1/Wpar^2)^(1/2))*(Sres - 1))/(2*(-1/Wpar^2)^(1/2));
        else
            CSw1p=Sres*Tsol-(Sres*Tsol - (pi^(1/2)*erfi((Tsol - Tliq)*(-1/Wpar^2)^(1/2))*(Sres - 1))/(2*(-1/Wpar^2)^(1/2)));
            CSw1=CSw1p + Sres*Tliq - (pi^(1/2)*erfi((Tliq - Tliq)*(-1/Wpar^2)^(1/2))*(Sres - 1))/(2*(-1/Wpar^2)^(1/2));
            CSw2=CSw1-(Tliq);
            ISw=CSw2+Tgau;
        end
    case 0 % Linear FTC
        Btem=(Sres-1)/Mpar;
        if Tgau>=Btem
            CSw1=Sres*Btem-(Mpar*Btem^2/2+Btem);
            ISw=CSw1+Mpar*Tgau^2/2+Tgau;
            if Tgau>Tliq
                CSw2=(CSw1+Mpar*Tliq^2/2+Tliq)-Tliq;
                ISw=CSw2+Tgau;
            end
        else
            ISw=Sres*Tgau;
        end
    case 2 % Sinosodial FTC
        if Tgau<=Tsol
            ISw=Sres*Tgau;
        elseif Tgau<Tliq && Tgau>Tsol
            CSw1=Sres*Tsol - ((Tsol*(Sres + 1))/2 + ((Tliq - Tsol)*(Sres - 1))/(2*pi));
            ISw=CSw1 + (Tgau*(Sres + 1))/2 + (cos((pi*(Tliq - 2*Tgau + Tsol))/(4*(Tliq - Tsol)))^2*(Tliq - Tsol)*(Sres - 1))/pi;
        else
            CSw1p=Sres*Tsol - ((Tsol*(Sres + 1))/2 + ((Tliq - Tsol)*(Sres - 1))/(2*pi));
            CSw1=CSw1p + (Tliq*(Sres + 1))/2 + ((Tliq - Tsol)*(Sres - 1))/(2*pi);
            CSw2=CSw1-(Tliq);
            ISw=CSw2+Tgau;
        end
        
end
end

%TODO: check the constants coming out of the integral of g(T)