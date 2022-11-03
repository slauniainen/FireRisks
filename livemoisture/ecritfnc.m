function E=ecritfnc(P0,L,a,bT,s,k0,Z50,kexp,h0,hs)

% P0 is base water pressure in MPa
% E is maximum transpiration in mmol/(mL^2-s) for given parameters and P0

% H is tree height in m
% L is path length in m
% a is vulnerability shape factor in 1/MPa

% assumes vulnerability curve f(P)=1/[1+exp(-a(P-b(z)))]
% b=bT+s*(L-z)
% bT is P50 at top of stem
% s is slope of P50 vs z in MPa/m
%   positive s means less negative P50 at bottom

% includes gravity

% kstar is hydraulic conductivity in kg/(mS-s-MPa)
% kstar=k0/[1+(z/(Z50*L))^kexp]
% k0 is kstar at bottom
% Z50 is fraction from bottom where kstar=.5*k0
% kexp is a sensitivity factor

% h is Huber ratio in cmS^2/mL^2
% h=h0+hs*z
% h0 is Huber ratio at bottom
% hs is positive slope
%   negative hs means Huber ratio decreases with height

% GL, 2018/02/24

% Y=1/flc is the primary dependent variable for the calculations

E1 = 1;
E2 = 3;
tol = 0.01;

z50 = Z50*L;
b0 = b(0); 
Y0 = exp(a*(P0-b0))+1;

E = E1;
[~,Y] = ode45(@yprime,[-L,0],1);
y1 = Y(end)-Y0;
E = E2;
[~,Y] = ode45(@yprime,[-L,0],1);
y2 = Y(end)-Y0;
if abs(y1)<abs(y2)
    E = E1;
    y = y1;
    Eold = E2;
    yold = y2;
else
    E = E2;
    y = y2;
    Eold = E1;
    yold = y1;
end    
iter=1;
while (abs(y)>tol) && iter<2000
    m = (E-Eold)/(y-yold);
    Eold = E;
    yold = y;
    E = Eold-m*yold;
    [~,Y] = ode45(@yprime,[-L,0],1);
    y = Y(end)-Y0;
    iter=iter+1;
end

% if iter==2000
%     E=-1;
% end

    function f=yprime(x,Y)
        f = -a*(B(-x)*(1-Y)-E*R(-x)*Y);
    end

    function y=R(z)
        y = .18./(h(z)*kstar(z));
    end

    function y=B(z)
        y = -s+.01;
    end

    function y=b(z)
        y = bT+s*(L-z);
    end

    function y=h(z)
        y = h0-hs*z;
    end

    function y=kstar(z)
        if kexp==0
            y = k0;
        else
            y = k0./(1+(z/z50)^kexp);
        end
    end
end



