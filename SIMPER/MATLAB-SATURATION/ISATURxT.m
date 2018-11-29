function ISw=ISATURxT(Tgau,Tsol,Tliq,Sres,Rsat,Wpar,Mpar)
switch Rsat
    case 1
        if Tgau<=Tsol
            ISw=Sres*Tgau^2/2;
        elseif Tgau<Tliq && Tgau>Tsol
            CSw1=Sres*Tsol^2/2-((Sres*Tsol^2)/2 + (Wpar^2*exp(-(Tsol - Tliq)^2/Wpar^2)*(Sres - 1))/2 + ...
                (Tliq*pi^(1/2)*erf((Tsol - Tliq)*(-1/Wpar^2)^(1/2)*1i)*(Sres - 1)*1i)/(2*(-1/Wpar^2)^(1/2)));
            ISw=CSw1+(Sres*Tgau^2)/2 + (Wpar^2*exp(-(Tgau - Tliq)^2/Wpar^2)*(Sres - 1))/2 + ...
                (Tliq*pi^(1/2)*erf((Tgau - Tliq)*(-1/Wpar^2)^(1/2)*1i)*(Sres - 1)*1i)/(2*(-1/Wpar^2)^(1/2));
        else
            CSw1p=Sres*Tsol^2/2-((Sres*Tsol^2)/2 + (Wpar^2*exp(-(Tsol - Tliq)^2/Wpar^2)*(Sres - 1))/2 + ...
                (Tliq*pi^(1/2)*erf((Tsol - Tliq)*(-1/Wpar^2)^(1/2)*1i)*(Sres - 1)*1i)/(2*(-1/Wpar^2)^(1/2)));
            CSw1=CSw1p+(Sres*Tliq^2)/2 + (Wpar^2*exp(-(Tliq - Tliq)^2/Wpar^2)*(Sres - 1))/2 + ...
                (Tliq*pi^(1/2)*erf((Tliq - Tliq)*(-1/Wpar^2)^(1/2)*1i)*(Sres - 1)*1i)/(2*(-1/Wpar^2)^(1/2));
            CSw2=CSw1-(Tliq^2/2);
            ISw=CSw2+Tgau^2/2;
        end
    case 0
        Btem=(Sres-1)/Mpar;
        if Tgau>=Btem
            CSw1=Sres*Btem^2/2-(Mpar*Btem^3/3+Btem^2/2);
            ISw=CSw1+Mpar*Tgau^3/3+Tgau^2/2;
            if Tgau>Tliq
                CSw2=(CSw1+Mpar*Tliq^3/3+Tliq^2/2)-Tliq^2/2;
                ISw=CSw2+Tgau^2/2;
            end
        else
            ISw=Sres*Tgau^2/2;
        end
    case 2
       if Tgau<=Tsol
            ISw=Sres*Tgau^2/2;
        elseif Tgau<Tliq && Tgau>Tsol
            CSw1=Sres*Tsol^2/2 - (Tsol^2*(Sres/4 + 1/4) + ((Tliq - Tsol)^2*(Sres - 1))/(2*pi^2));
            ISw=CSw1 + Tgau^2*(Sres/4 + 1/4) + (sin((pi*(Tliq - 2*Tgau + Tsol))/(2*(Tliq ...
                - Tsol)))*(Tliq - Tsol)^2*(Sres - 1))/(2*pi^2) + (Tgau*cos((pi*(Tliq ...
                - 2*Tgau + Tsol))/(2*(Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1))/(2*pi);
        else
            CSw1p=Sres*Tsol^2/2 - (Tsol^2*(Sres/4 + 1/4) + ((Tliq - Tsol)^2*(Sres - 1))/(2*pi^2));
            CSw1=CSw1p + Tliq^2*(Sres/4 + 1/4) - ((Tliq - Tsol)^2*(Sres - 1))/(2*pi^2);
            CSw2=CSw1-(Tliq^2/2);
            ISw=CSw2+Tgau^2/2;
        end 
end
end

%TODO: check the constants coming out of the integral of g(T)