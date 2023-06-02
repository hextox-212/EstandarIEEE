clc
clear
% DATOS FISICOS DEL CONDUCTOR
D=25.35;
Rl=0.0756;
Rh=0.0937;
e=0.5;
a=0.5;
% DATOS MEDIAMBIENTALES
Ta=30;
V=0.61;
He=2700;
Lat=16.5;
N=270;
h=12;
Tcrange=[]
% COMVECCION
for Tc=(40:80)
    Tcrange=[Tcrange,Tc];
    Uf=(1.458*10^-6*((Tc+Ta)/2+273)^1.5/((Tc+Ta)/2+383.4));
    p=(1.293-1.525*10^-4*He+6.379*10^-9*He^2)/(1+0.00367*(Tc+Ta)/2);
    Kf=2.424*0.01+7.477*10^-5*(Tc+Ta)/2-4.407*10^-9*((Tc+Ta)/2)^2;
    Qc1=(1.01+0.0372*(p*D*V/Uf)^0.52)*Kf*(Tc-Ta);
    Qc2=(0.0119*(D*p*V/Uf)^.6)*Kf*(Tc-Ta);
    qc=max(Qc1,Qc2);
    % RADIACION
    qr=0.0178*D*e*(((273+Tc)/100)^4-((273+Ta)/100)^4);
    % RADIACION SOLAR
    w=15*(h-12);
    d=23.4583*sind((284+N)/365*360);
    Hc=asin(cosd(Lat)*cosd(d)*cosd(w)+sind(Lat)*sind(d))*180/3.14159;
    X=(sind(w)/(sind(Lat)*cosd(w)-cosd(Lat)*tand(d)));
        if X*w<0
            C=0;
        else
            C=180;
        end
    Zc=C+180/pi*atan(X);
    O=acos(cosd(Hc)*cosd(Zc-90))*180/3.14159;
    Ksol=1+1.148*10^-4*He-1.108*10^-8*He^2;
    Qs=-42.2391+63.4044*Hc-1.922*Hc^2+3.46921*0.01*Hc^3-3.61118*10^-4*Hc^4+1.94318*10^-6*Hc^5-4.07608*10^-9*Hc^6;
    Qse=Ksol*Qs;
    A=D/10^3;
    qs=a*Qse*sind(O)*A;
    % EFECTO JOULE
    R=Rl+((Rh-Rl)/(75-20))*(Tc-20);
    % AMPACIDAD
    I(Tc)=((qc+qr-qs)/R*10^3)^.5;
end
Irange=I(1,40:80);
plot(Tcrange,Irange)

