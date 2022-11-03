
ETmm=5; %transipiration rate, per unit ground area (mm/day) - from SpaFy
P0ref=-1; %water potential at the bottom of the canopy MPa
LAI=3; %m2leaf/m2ground
Lref=25; %m canopy height
daylength=14; %hrs (can be substituted with a DOY-specific value)

a=1.07; %1.07 in Table 3 Couvrer; 1.6: value for Pinus Ponderosa in Couvrer
bTs=[-5.79,-0.0192];  %Couvrer Table 3, does not match exactly the value for either Psudotsuga or Pinus in Table 2 though
bs=0.08; %MPa/m: average value for conifer in Couvrer
h0s=[4.3,-0.050]; %dependence of Huber value on L; Couvrer Table 3; return Huber in c2sap/m2leaf
k0s=[1.4,0.11]; %dependece of saturated hydraulic conductivity at the base of the tree (z = 0) on L; Couvrer, Table 3; for Pinus, [1.1,0.11]

load("moisture_release_fit.mat")

%pressure-volume curve dataset
j=10; %Picea needle
%pressure
[P,z]=tree_water_potential(ETmm, daylength, P0ref, Lref, LAI, k0s, h0s, a, bTs, bs);
figure(9)
subplot(2,1,1)
plot(P,z,'-k')
xlabel('P (MPa)');ylabel('z (m)')
%RWC
RWC=(afitInv(j,1)-P.*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P);
    figure(9)
    subplot(2,1,2)
    plot(RWC,z,'-');hold on
figure(9)
subplot(2,1,2)
hold off
xlabel('RWC');ylabel('z (m)')


pause
%%
%sensitivity analysis on pressure - volume curve
[P,z]=tree_water_potential(ETmm, daylength, P0ref, Lref, LAI, k0s, h0s, a, bTs, bs);
figure(10)
subplot(2,1,1)
plot(P,z,'-k')
xlabel('P (MPa)');ylabel('z (m)')

% Relative water content (fraction)
for j=1:length(afitInv(:,1))
    RWC=(afitInv(j,1)-P.*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P);
    figure(10)
    subplot(2,1,2)
    plot(RWC,z,'-');hold on
end
figure(10)
subplot(2,1,2)
hold off
xlabel('RWC');ylabel('z (m)')


%%
%sensitivity analysis on ET
ets=[1:2:8]; vals=ets;
for i=1:length(vals)
    et=vals(i)
    [P,z]=tree_water_potential(et, daylength, P0ref, Lref, LAI, k0s, h0s, a, bTs, bs);

    figure(11)
    subplot(3,3,1)
    plot(P,z,'Color',[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)]);hold on
    xlabel('P (MPa)');ylabel('z (m)')
    title('ET increases with darkening gray')
    for j=1:length(afitInv(:,1))
        col=[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)];
        if length(find(needles==j))>0
            symb='^';
        else
            symb='o';
        end
        subplot(3,3,4)
        RWC=(afitInv(j,1)-P(end).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(end));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on%'Color',[1-j/length(afitInv(:,1)) 1-j/length(afitInv(:,1)) 1-j/length(afitInv(:,1))]);
        title('triangles - needles')
        xlabel('dataset #')
        ylabel('RWC (h_c)')
        subplot(3,3,7)
        RWC=(afitInv(j,1)-P(round(length(z)/2)).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(round(length(z)/2)));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on
        xlabel('dataset #')
        ylabel('RWC (h_c/2)')
        if length(find(piceas==j))>0
            subplot(3,3,4)
            plot([j,j],[0,1],':r')
            subplot(3,3,7)
            plot([j,j],[0,1],':r')
        end
    end
end

%sensitivity analysis on P0
P0s=[-1.51:0.2:-0.03]; vals=P0s;
for i=1:length(vals)
    P0=vals(i)
    [P,z]=tree_water_potential(ETmm, daylength,P0, Lref, LAI, k0s, h0s, a, bTs, bs);

    figure(11)
    subplot(3,3,2)
    plot(P,z,'Color',[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)]);hold on
    xlabel('P (MPa)');ylabel('z (m)')
    title('P0 increases with darkening gray')
    for j=1:length(afitInv(:,1))
        col=[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)];
        if length(find(needles==j))>0
            symb='^';
        else
            symb='o';
        end
        subplot(3,3,5)
        RWC=(afitInv(j,1)-P(end).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(end));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on%'Color',[1-j/length(afitInv(:,1)) 1-j/length(afitInv(:,1)) 1-j/length(afitInv(:,1))]);
        title('red lines: Picea')
        xlabel('dataset #')
        ylabel('RWC (h_c)')
        subplot(3,3,8)
        RWC=(afitInv(j,1)-P(round(length(z)/2)).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(round(length(z)/2)));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on
        xlabel('dataset #')
        ylabel('RWC (h_c/2)')
        if length(find(piceas==j))>0
            subplot(3,3,5)
            plot([j,j],[0,1],':r')
            subplot(3,3,8)
            plot([j,j],[0,1],':r')
        end
    end
end

%sensitivity analysis on P0
Ls=[5:10:45]; vals=Ls;
for i=1:length(vals)
    L=vals(i)
    [P,z]=tree_water_potential(ETmm, daylength, P0ref, L, LAI, k0s, h0s, a, bTs, bs);

    figure(11)
    subplot(3,3,3)
    plot(P,z,'Color',[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)]);hold on
    xlabel('P (MPa)');ylabel('z (m)')
    title('L increases with darkening gray')
    for j=1:length(afitInv(:,1))
        col=[1-i/length(vals) 1-i/length(vals) 1-i/length(vals)];
        if length(find(needles==j))>0
            symb='^';
        else
            symb='o';
        end
        subplot(3,3,6)
        RWC=(afitInv(j,1)-P(end).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(end));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on;
        xlabel('dataset #')
        ylabel('RWC (h_c)')
        subplot(3,3,9)
        RWC=(afitInv(j,1)-P(round(length(z)/2)).*(1+afitInv(j,2)))./(afitInv(j,1)-afitInv(j,2).*P(round(length(z)/2)));
        plot(j,RWC,symb,'MarkerEdgeColor',col);hold on
        xlabel('dataset #')
        ylabel('RWC (h_c/2)')
        if length(find(piceas==j))>0
            subplot(3,3,6)
            plot([j,j],[0,1],':r')
            subplot(3,3,9)
            plot([j,j],[0,1],':r')
        end
    end
end