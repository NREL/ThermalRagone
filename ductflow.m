function [Nu_T, Nu_H, f] = ductflow(Re_bar, Re, Pr, LoverD, Asp, relRough)
if Asp<0 || Asp >1
    disp('Error - Aspect Ratio must be Between 0 and 1')
end
if relRough<0 || relRough>0.05 
    disp('Error - Relative Roughness must be between 0 and 0.05')
end

if Re_bar < 2300                                                           	%Laminar Flow
    %Friction Factor
    f_fd = 24*(1-1.3553*Asp+1.9467*Asp^2-1.7012*Asp^3+0.9564*Asp^4-...      %Fully Developed Friction Factor (EES)
        0.2537*Asp^5);
    x_plus = LoverD./Re;                                                    %Dimensionless Length (-)
    fR=(3.44./sqrt(x_plus)+(1.25./(4.*x_plus)+f_fd-3.44./sqrt(x_plus))./... %Rectangular Duct Friction Factor
        (1+0.00021.*x_plus.^(-2)));
    f = 4.*fR./Re;                                                          %Equations for Friction Factor from Eqn 3.158 in Kakac, Shah and Aung (EES)
    
    %Nusselt Number
    %Fully Developed Flow
    Nu_T_fd = 7.541*(1-2.610*Asp+4.970*Asp^2-5.119*Asp^3+2.702*Asp^4-...    %Eq. 3.159 in Kakac, Shah and Aung (EES)
        0.548*Asp^5);
    Nu_H_fd = 8.235*(1-2.0421*Asp+3.0853*Asp^2-2.4765*Asp^3+...             %Eq. 3.161 in Kakac, Shah and Aung (EES)
        1.0578*Asp^4-0.1861*Asp^5); 
    
    %Developing Flow (Wibulswas (1966) - EES)
    x_star = LoverD./(Re.*Pr);                                              %Dimensionless Length (-)
	lnx_star = log(x_star);                                                 %Natural Log of Dimensionless Length (-)
    
        %Constant T, Pr = 0.72
        a_T = 0.0357122+0.460756236*Asp-0.314865737*Asp^2;                  %Constant Temperature Constant 1
        %b_T = 0.602877+0.0337485/(0.1+Asp)+0.0377031*Asp - why in EES?
        if (Asp<0.167)
            b_T = 0.940362+Asp*(0.739606-0.940362)/0.167;                   %Constant Temperature Constant 2
        else
            b_T = 0.801105912 - 0.419264242*Asp + 0.293641181*Asp^2;        %Constant Temperature Constant 2
        end
        DNusselt_T = a_T.*exp(-b_T.*lnx_star);                              %Constant Temperature Nusselt Number for Pr = 0.72
    
        %Constant H, Pr=0.72
        a_H = 0.113636994 + 0.712134212*Asp - 0.392104717*Asp^2;            %Constant Heat Flux Constant 1
        if (Asp<0.25)
            b_H = 0.940362+Asp*(0.699466-0.940362)/0.25;                    %Constant Heat Flux Constant 2
        else
            b_H = 0.774133656 - 0.350363736*Asp + 0.198543081*Asp^2;        %Constant Heat Flux Constant 2
        end
        DNusselt_H = a_H.*exp(-b_H.*lnx_star);                              %Constant Heat Flux Nusselt Number for Pr = 0.72
	
        %Correct for Pr
        if (Pr>0.72)	 
            DNurat = 0.6847+0.3153*exp(-1.26544559.*(log(Pr)-log(0.72)));   %Nusselt Number Correction
        else
            DNurat=1.68-0.68.*exp(0.32.*(log(Pr)-log(0.72)));               %Nusselt Number Correction
        end
        
        Nu_T = Nu_T_fd + DNurat.*DNusselt_T;                                %Constant Temperature Nusselt Number for Constant Temperature Boundary
        Nu_H = Nu_H_fd + DNurat.*DNusselt_H;                                %Constant Heat Flux Nusselt Number for Constant Heat Flux Boundary
        
elseif Re_bar > 3000                                                        %Turbulent Flow
    %Friction Factor
    f_fd=(-0.001570232./log(Re)+0.394203137./log(Re).^2+...                 %From Li, Seem, and Li, "A New Explicity Equation for Accurate Friction Factor Calculation for Smooth Tubes, 2011
        2.534153311./log(Re).^3).*4; 
	if relRough>1e-5 
        f_fd = (-2.*log10(relRough./3.71-1.975./Re.*log((relRough/3.93)...  %Offor and Alabi, Advances in Chemical Engineering and Science, 2016, 6, 237-245
            .^1.092+7.627./(Re+395.9)))).^(-2);
    end
    f = f_fd.*(1+(1/LoverD)^0.7);                                           %Account for developing flow (approximate)

    %Nusselt Number
    Nusselt_L=((f_fd./8).*(Re-1000).*Pr)/(1+12.7.*sqrt(f_fd./8)...          %Gnielinski, V.,, Int. Chem. Eng., 16, 359, 1976
        .*(Pr.^(2/3)-1));	
    if (Pr<0.5)
		Nusselt_L_lp = 4.8 + 0.0156.*Re.^0.85.*Pr.^0.93;                    %Notter and Sleicher, Chem. Eng. Sci., Vol. 27, 1972}
		if (Pr<0.1)
			Nusselt_L = Nusselt_L_lp;
        else
			Nusselt_L=Nusselt_L_lp+(Pr-0.1).*(Nusselt_L-Nusselt_L_lp)./0.4;
        end
    end
	Nu = Nusselt_L.*(1+(1/LoverD)^0.7);                                     %Nusselt Number Accounting for Developing Flow (approximate)
    Nu_T = Nu;                                                              %Nusselt Number not a Function of Boundary Conditions for Turbulent Flow
    Nu_H = Nu;                                                              %Nusselt Number not a Function of Boundary Conditions for Turbulent Flow
else
    f = nan;
    Nu_T = nan;
    Nu_H = nan;
end
    
    
    