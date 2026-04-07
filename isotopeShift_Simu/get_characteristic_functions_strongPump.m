function cfs = get_characteristic_functions_strongPump(w_L)
% 计算强光极限下的各种近似解析式
% (#Translate# Compute various approximate analytic expressions in the strong pumping limit.)

%     w_L = linspace(0,1.0,200);
    tmp_erfi = Erfi(1/2/sqrt(2)./w_L);
    F5 = -sqrt(2*pi)/3*w_L.^2  .*  tmp_erfi.*(  -12*w_L.*exp(1/8./w_L.^2)  + sqrt(2*pi)*(1+12*w_L.^2).*tmp_erfi    );


    I00x2 = Ix2_fun(0,0); I00x2 = I00x2(w_L);
    % sum (n)
    summer = 0;
    for n = 1:20
        I0nx2 = Ix2_fun(0,n);  I0nx2 = I0nx2(w_L);
        summer = summer + I0nx2.^2/n^4;
    end
    F6 = -2*I00x2.^4 .* summer;
    % sum (m,n)
    summer = 0;
    for m = 1:20
        I0mx2 = Ix2_fun(0,m);  I0mx2 = I0mx2(w_L);
        for n = 1:20
            I0nx2 = Ix2_fun(0,n);  I0nx2 = I0nx2(w_L);
            summer = summer + (I0mx2.*I0nx2).^2/(m^2+n^2)^2;
        end
    end
    F6 = F6 - I00x2.^2 .* summer;
    % sum (n,n')
    summer = 0;
    for m = 1:20
        I0mx2 = Ix2_fun(0,m);  I0mx2 = I0mx2(w_L);
        for n = 1:20
            I0nx2 = Ix2_fun(0,n);  I0nx2 = I0nx2(w_L);
            Imnx2 = Ix2_fun(m,n);  Imnx2 = Imnx2(w_L);
            summer = summer + (I0mx2.*I0nx2).*(I00x2.*Imnx2 + I0nx2.*I0mx2)/(m^2*n^2);
        end
    end
    F6 = F6 + 2*I00x2.^2 .* summer;
    % sum (m,n,n')
    summer = 0;
    for n1 = 1:20
        I0n1x2 = Ix2_fun(0,n1);  I0n1x2 = I0n1x2(w_L);
        for n = 1:20
            I0nx2 = Ix2_fun(0,n);  I0nx2 = I0nx2(w_L);
            In1nx2 = Ix2_fun(n1,n);  In1nx2 = In1nx2(w_L);
            for m = 1:20
                I0mx2 = Ix2_fun(0,m);  I0mx2 = I0mx2(w_L);
                summer = summer + I0nx2.*I0n1x2.* I0mx2.^2 .* In1nx2 / ((n^2+m^2)*n1^2);
            end
        end
    end
    F6 = F6 + 4*I00x2 .* summer;
    % sum (m,n,m',n')
    summer = 0;
    for n1 = 1:12
        I0n1x2 = Ix2_fun(0,n1);  I0n1x2 = I0n1x2(w_L);
        for n = 1:12
            I0nx2 = Ix2_fun(0,n);  I0nx2 = I0nx2(w_L);
            In1nx2 = Ix2_fun(n1,n);  In1nx2 = In1nx2(w_L);
            for m1 = 1:12
                I0m1x2 = Ix2_fun(0,m1);  I0m1x2 = I0m1x2(w_L);
                for m = 1:12
                    I0mx2 = Ix2_fun(0,m);  I0mx2 = I0mx2(w_L);
                    Im1mx2 = Ix2_fun(m1,m);  Im1mx2 = Im1mx2(w_L);
                    summer = summer + I0nx2.*I0n1x2.* I0mx2.*I0m1x2 .* In1nx2.* Im1mx2 / ((n^2+m^2)*(n1^2+m1^2));
                end
            end
        end
    end
    F6 = F6 + summer;
    
    I02x2 = Ix2_fun(0,2); I02x2 = I02x2(w_L);
    I22x2 = Ix2_fun(2,2); I22x2 = I22x2(w_L);
    F6_approx = I00x2.^6.*(   1/8*(I22x2./I00x2 - 1).*(I02x2./I00x2).^2  +1/64*(7 + 8*I22x2./I00x2 + (I22x2./I00x2).^2).*(I02x2./I00x2).^4   );
    
    cfs = struct();
    cfs.F5 = F5;
    cfs.F6 = F6;
    cfs.F6_approx = F6_approx;
end

% 给定 (m,p), 返回 近似解析表达式 函数 Impx2(w/L)。近似到 lambda 的零阶项
% (#Translate# Given (m, p), returns the approximate analytic expression function Impx2(w/L). Approximation is to zeroth order in lambda.)
function fun = Ix2_fun(m,p)
    flag_use_faddeeva = true;
    if m==0 && p==0
        fun = @(w_L)sqrt(2*pi)*w_L.*Erfi(1/2/sqrt(2)./w_L);
    elseif m==0 || p==0
        if p==0
            p=m;
        end
        if flag_use_faddeeva
            % use Faddeeva_w() function to simplify numerical calculation
            fun = @(w_L)2*sqrt(pi)*w_L.*exp(1/8./w_L.^2).* real( 1i*exp(-0.5i*p*pi) .* Faddeeva_w(-1/2/sqrt(2)./w_L + 1i*p*pi*w_L/sqrt(2)) ).*cos(p*pi/2);
        else
            fun = @(w_L)sqrt(pi)*w_L.*exp(+p^2*pi^2*w_L.^2/2).*(  Erfi(1/2/sqrt(2)./w_L - 1i*p*pi*w_L/sqrt(2)) + Erfi(1/2/sqrt(2)./w_L + 1i*p*pi*w_L/sqrt(2))   ).*cos(p*pi/2);
        end
    else
        if flag_use_faddeeva
            fun = @(w_L)0;
            for mp = [m+p, m-p]
                fun = @(w_L)fun(w_L) + cos(pi/2*mp).*real( 1i*  exp( -0.5i*pi*abs(mp) ).*Faddeeva_w( - 1./(2*sqrt(2)*w_L) + 1i*pi*abs(mp)*w_L/sqrt(2))  );
            end
            fun = @(w_L)fun(w_L).*sqrt(2*pi).*w_L.*exp(1/8./w_L.^2);
        else
            fun = @(w_L)0;
            for s1 = [-1,+1]
                for s2 = [-1,+1]
                    fun = @(w_L)fun(w_L) + cos(pi/2*(s1*m+s2*p)).*exp((s1*m+s2*p).^2*pi^2*w_L.^2/2).*Erfi(1./(2*sqrt(2)*w_L) + 1i*pi*(s1*m+s2*p)*w_L/sqrt(2));
                end
            end
            fun = @(w_L)fun(w_L).*sqrt(pi/2).*w_L;
        end
    end
end

function val = Erfi(z)
    val = erfz(1i*z)/1i;
end