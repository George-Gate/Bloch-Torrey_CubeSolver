function cfs = get_characteristic_functions_weakPump(w_L, OD)
% 计算弱光极限下的各种近似解析式
% (#Translate# Compute various approximate analytic expressions in the weak pumping limit.)

%     w_L = linspace(0,1.0,200);
    F1 = sqrt(2*pi)/3*w_L.^2 .* erf(1/2/sqrt(2)./w_L).*(    12*w_L.*exp(-1/8./w_L.^2) + sqrt(2*pi)*(1-12*w_L.^2).*erf(1/2/sqrt(2)./w_L)    );
    F2 = 2*pi/3*w_L.^2.*erf( 1/2/sqrt(2)./w_L ).^2;
    F3 = (  2*pi*w_L.^2 .*erf(1/2/sqrt(2)./w_L).^2  ).^3;
    F4 = 2*(  2*pi*w_L.^2 .*erf(1/2/sqrt(2)./w_L).^2  ).^2.*(pi*w_L.^2 .* exp(-4*pi^2*w_L.^2).*( erfz(1/2/sqrt(2)./w_L + 1i*sqrt(2)*pi*w_L) + erfz(1/2/sqrt(2)./w_L - 1i*sqrt(2)*pi*w_L)  ).^2);

%     OD = linspace(0,15,1000);
    G1 = (1-exp(-OD))./OD;
    G2 = -(  (12+OD.*OD).*(1-exp(-OD)) - 6*OD.*(1+exp(-OD))    )./(OD.^3);

    I00z1 = Iz1_fun(0,0); I00z1 = I00z1(OD); % I00z1 = ((1-exp(-OD))./(OD));
    % I0pz1_2 = @(p) (  sqrt(2)*OD.*(1+(-1)^(p+1)*exp(-OD))./(OD.^2 + p^2*pi^2)  ).^2;
    G3 = 0;  
    G4 = (I00z1).^2/16;
    for p = 1:100
        I0pz1 = Iz1_fun(0,p);
        G3 = G3 + I0pz1(OD).^2/p^4;
        G4 = G4 + I0pz1(OD).^2/(p^2+4)^2;
    %     G1 = G1 + I0pz1_2(p)/p^4;
    %     G2 = G2 + I0pz1_2(p)/(p^2+4)^2;
    end
    G3 = G3.*I00z1;
    G4 = G4.*I00z1;

    G5 = 0;  G6 = 0;
    for p1 = 0:100
        I0p1z1 = Iz1_fun(0,p1);  I0p1z1 = I0p1z1(OD);
        for p2 = 1:100
            I0p2z1 = Iz1_fun(0,p2);   I0p2z1 = I0p2z1(OD);
            Ip1p2z1 = Iz1_fun(p1,p2);  Ip1p2z1 = Ip1p2z1(OD);
            if p1 > 0
                G5 = G5 + I0p1z1.*Ip1p2z1.*I0p2z1/(p1^2*p2^2);
            end
            G6 = G6 + I0p1z1.*Ip1p2z1.*I0p2z1/((4+p1^2)*p2^2);
        end
    end
    G6 = 2*G6;

    G7 = G3 - G5;
    G8 = G4 - G6;
    
    cfs = struct();
    cfs.F1 = F1;  cfs.F2 = F2;   cfs.F3 = F3;  cfs.F4 = F4;
    cfs.G1 = G1;  cfs.G2 = G2;   cfs.G3 = G3;  cfs.G4 = G4;
    cfs.G5 = G5;  cfs.G6 = G6;   cfs.G7 = G7;  cfs.G8 = G8;
end

% 给定 (m,p), 返回 近似解析表达式 函数 Impz1(OD)。近似到 lambda 的零阶项
% (#Translate# Given (m, p), it returns the approximate analytic expression as the function Impz1(OD), accurate to the zeroth order in lambda.)
function fun = Iz1_fun(m,p)
    if m==0 && p==0
        fun = @(OD)(1-exp(-OD))./OD;
    elseif m==0 || p==0
        if p==0
            p=m;
        end
        fun = @(OD)sqrt(2)*OD.*(1+(-1)^(p+1).*exp(-OD))./(OD.*OD+p^2*pi^2);
    else
        fun = @(OD)2*OD.*( 1+(-1)^(m+p+1).*exp(-OD) ).*(  OD.*OD + (m^2+p^2)*pi^2  )./( OD.*OD+(m-p)^2*pi^2 )./( OD.*OD+(m+p)^2*pi^2 );
    end
end