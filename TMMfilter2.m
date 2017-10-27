% Transfer Matrix Methods to calculate R T A for filter
% Hao Wang @ SUTD
% 2017 October 26


%inputs:
%angle of incidence: thetai
%wavelength of incident light: lambda
%thicknesses of the layers: h
%refractive index of the layers: n (may be absorbing and dispersive, but this may require a subfunction)
%polarization: s or p (need one calculation for each for unpolarized light)
clc
clear all
close all

tic
% thetai=0:0.1:90; %angle of incidence (degrees)
thetai=0; %angle of incidence (degrees)
lambda=300:1:1000; %vacuum wavelength (nm)

central_wave=532;
ha=central_wave*(1+1/4)/1.55;
hb=central_wave*(1+1/4)/1;

%% the layer numbers
% the first and last layers are environment, and the layer numbers should
% be 2+material layer number
lnum=9;
h=zeros(1,lnum);
h(1)=NaN;
h(lnum)=NaN;
for gg=1:(length(h)-2)/2
    h(2*gg)=ha;
    h(2*gg+1)=hb;
end
h(lnum-1)=ha;
h((lnum+1)/2)=2*hb;

pol=1; %polarization, 1 for p and 0 for s

%% Materials

% % for dispersion materials
% filename='F:\Code\Plasmonic\material_data/Ag.txt';
% L1=dlmread(filename,'',1);
% n_interp=interp1(L1(:,1), L1(:,2), lambda, 'linear', 'extrap');
% k_interp=interp1(L1(:,1), L1(:,3), lambda, 'linear', 'extrap');
% L11=real((n_interp+i*k_interp).^(1/2))+i*imag((n_interp+i*k_interp).^(1/2));
% size(L11);
% 
% filename2='F:\Code\Plasmonic\material_data/Al2O3.txt';
% L2=dlmread(filename2,'',1);
% n_interp2=interp1(L2(:,1), L2(:,2), lambda, 'linear', 'extrap');
% k_interp2=interp1(L2(:,1), L2(:,3), lambda, 'linear', 'extrap');
% L12=real((n_interp2+i*k_interp2).^(1/2))+i*imag((n_interp2+i*k_interp2).^(1/2));
% size(L12);




n=zeros(length(h),length(lambda));
[nn,m]=size(n);

% substrate refractive index, light incidences from the substrate
% n=1.46 for nanoscribe glass, n=1.52 for ordinarily used fused silica
% n=1.55 for nanoscribe Ip-Dip

n(1,:)=1.46;
n(length(h),:)=1;
for kk=1:(length(h)-2)/2
    % for dispersion materials
%     n(2*kk,:)=L11;
%     n(2*kk+1,:)=L12;
    % for non-dispersion materials
    n(2*kk,:)=1.55;
    n(2*kk+1,:)=1;
end
n(length(h)-1,:)=1.55;

%% TMM

for b=1:length(lambda)
    
%Snell's law:
    theta(1)=thetai*pi/180;

    for a=1:nn-1
    theta(a+1)=real(asin(n(a,b)/n(a+1,b)*sin(theta(a))))-1i*abs(imag(asin(n(a,b)/n(a+1,b)*sin(theta(a)))));
    end
    

%Fresnel coefficients:
    if pol==0 %formulas for s polarization
    
        for a=1:nn-1
            n(a,b);
            Fr(a)=(n(a,b)*cos(theta(a))-n(a+1,b)*cos(theta(a+1)))/(n(a,b)*cos(theta(a))+n(a+1,b)*cos(theta(a+1)));
            Ft(a)=2*n(a,b)*cos(theta(a))/(n(a,b)*cos(theta(a))+n(a+1,b)*cos(theta(a+1)));
    
        end
%phase shift factors:

    delta=zeros((nn-2),length(lambda));
% for a=1:length(n)-2
%     for b=1:length(lambda)
%         delta(a,b)=2*pi*h(a+1)/lambda(b)*n(a+1)*cos(theta(a+1))
%     end
% end

        for a=1:nn-2
            delta(a,b)=2*pi*h(a+1)/lambda(b)*n(a+1,b)*cos(theta(a+1));
        end
% figure
% plot(lambda, delta(1,:))


%build up transfer matrix:


        M=[1,0;0,1]; %start with unity matrix
        for a=1:nn-2
            M=M*1/Ft(a)*[1,Fr(a);Fr(a),1]*[exp(-1i*delta(a,b)),0;0,exp(1i*delta(a,b))];
        end
    
        M=M*1/Ft(nn-1)*[1,Fr(nn-1);Fr(nn-1),1];
    
    
elseif pol==1 %formulas for p polarization
    
    for a=1:nn-1
        Fr(a)=(n(a,b)*cos(theta(a+1))-n(a+1,b)*cos(theta(a)))/(n(a,b)*cos(theta(a+1))+n(a+1,b)*cos(theta(a)));
        Ft(a)=2*n(a,b)*cos(theta(a))/(n(a,b)*cos(theta(a+1))+n(a+1,b)*cos(theta(a)));
    end

%phase shift factors:

    delta=zeros((nn-2),length(lambda));
% for a=1:length(n)-2
%     for b=1:length(lambda)
%         delta(a,b)=2*pi*h(a+1)/lambda(b)*n(a+1)*cos(theta(a+1))
%     end
% end

    for a=1:nn-2
        delta(a,b)=2*pi*h(a+1)/lambda(b)*n(a+1,b)*cos(theta(a+1));
    end
% figure
% plot(lambda, delta(1,:))


%build up transfer matrix:


    M=[1,0;0,1]; %start with unity matrix
    for a=1:nn-2
        M=M*1/Ft(a)*[1,Fr(a);Fr(a),1]*[exp(-1i*delta(a,b)),0;0,exp(1i*delta(a,b))];
    end
    
    M=M*1/Ft(nn-1)*[1,Fr(nn-1);Fr(nn-1),1];
end



%total Fresnel coefficients:
Frtot(b)=M(2,1)/M(1,1);
Fttot(b)=1/M(1,1);

%special case of single interface:
    if length(n)==2
        Frtot(b)=Fr(1);
        Fttot(b)=Ft(1);
    end

%total Fresnel coefficients in intensity:
FR(b)=(abs(Frtot(b)))^2;
FT(b)=(abs(Fttot(b)))^2*real(n(nn,b)*cos(theta(nn)))/real(n(1,b)*cos(theta(1)));
FA(b)=1-FR(b)-FT(b);
end
toc

% % Plot parameters
lw=3; font_axis=25; font_label=25; font_legend=18; ms=8;
colorred=[244 140 132]./255;
colorblue=[0 230 230]./255;

% 
% figure(1)
% clf;
% set(gcf,'units','inches','position',[0.25 2.5 11.3 7.9]);
% % set(gcf,'color','w','Renderer','painters','units','inches','position',[0.25 2.5 8.0 6.0]);
% 
% 
% plot(lambda,n_interp,'b','linewidth',lw);
% hold on
% plot(lambda,k_interp,'r','linewidth',lw);
% hold on
% plot(L1(:,1), L1(:,2),'ob','MarkerSize',ms,'MarkerEdgeColor','b','MarkerFaceColor',colorblue);
% hold on
% plot(L1(:,1), L1(:,3),'or','MarkerSize',ms,'MarkerEdgeColor','r','MarkerFaceColor',colorred);
% 
% xlim([lambda(1),lambda(length(lambda))])
% % lg=legend('interp Re(\epsilon)','interp Im(\epsilon)','J&C Re(\epsilon)','J&C Im(\epsilon)');
% lg=legend('interp Re(\epsilon) of Al','interp Im(\epsilon) of Al','Weaver & Frederikse Re(\epsilon) of Al','Weaver & Frederikse Im(\epsilon) of Al');
% set(lg,'FontSize',font_legend,'box','off');
% set(gca,'FontSize',font_axis,'LineWidth',lw-1,'XMinorTick','off','YMinorTick','off');
% xlabel('Wavelength (nm)','fontsize',font_label);
% ylabel('Permittivity','fontsize',font_label);

%plot results:
figure(2)
clf;
clf;
set(gcf,'units','inches','position',[0.25 2.5 11.3 7.9]);
hold on
plot(lambda,FT,'g','linewidth',lw)
plot(lambda,FR,'b','linewidth',lw)
plot(lambda,FA,'r','linewidth',lw)
xlim([lambda(1),lambda(length(lambda))])
ylim([0 1])

% title('Fresnel coefficients for transmission (blue), reflection (red) and absorption (green)')
title('RTA')
lg=legend('T', 'R','A');
set(lg,'FontSize',font_legend,'box','off');
set(gca,'FontSize',font_axis,'LineWidth',lw-1,'XMinorTick','off','YMinorTick','off');
xlabel('Wavelength (nm)','fontsize',font_label);
ylabel('RTA','fontsize',font_label);
box on

%% Phase
