
clear all vars 
clc
close all force
%load('D:\NIH brainprobe\DATA_Gap1xspan0.2yspan0.5Number8mm.mat');
%% 3D rayleigh integral phase array repurposed from Ultrasonic course project(Vortex beam simulation)
%% Innitialize plot setting
W_bar = waitbar(0,'Please wait...');
set(groot,'DefaultAxesFontSize', 20)
set(groot,'DefaultLineLineWidth', 2)
set(groot,'DefaultAxesLineWidth', 2)

%% Innitialize parameters 
F = 5e-3;         %Focal length
f = 2.5e6;        %Frequency of Ultrasonic wave
rho0 = 1000;        %Density of water
c0 = 1500;          %sound speed in water
K = 2*pi*f/c0;      %wavenumber
lambda = 2*pi/K;    %wavelength
F_offset_x=0;
F_offset_y=0e-2;
h.r=lambda/20;   % mesh size in the r direction
r=(0:h.r/2:6e-3);
u0 = 40e-9*2*pi*f;             %Particle velocity amplitude on the surface of the transducer in m/s 
Ray_dist=0.5*K*(1.7e-3/2)^2;
Material.rho0=rho0;
Material.c0=c0;
Material.K=K;
Material.u0=u0;
loss=2.5; %loss in the brain is 1dB/cm/MHz
NT=24;
a = 10e-3;          %Rate of Roll-off -r^2/a^2
ar=3*a;               %Radius of transducer
xspan_T=0.35e-3;
yspan_T=0.35e-3;
Gap=0.05e-3;


h.XT=lambda/12;       %mesh size for the transducer in the X direction
h.YT=lambda/12;      %mesh size for the transducer in the Y direction

h.X1=lambda/15;       %mesh size for the interrogation plane in the X direction
h.Y1=lambda/15;      %mesh size for the interrogation plane in the Y direction
r_y1=h.Y1;
r_x1=h.X1;

D_X1=12e-3;         %x span of the interrogation plane
D_Y1=40e-3;        %y span of the plane of interrogation

D_XT=5e-3;      % x span of the transducer
D_YT=40e-3;      % y span of the 
if mod(ceil(D_YT/Gap),2)==1
    D_YT=Gap*(ceil(D_YT/Gap)+1);
end





ka=sqrt(xspan_T*yspan_T/pi)*2*pi*f/c0
z0 = F;             %Distance between the transducer and the focal plane

%Rayleigh Integral from P_T to P_1 and add phase plate at P_1
%W_bar = waitbar(0,'Please wait...');
rho0=Material.rho0;
c0=Material.c0;
K=Material.K;
u0 = Material.u0; 

size_XT=floor(D_XT/h.XT);     %Size of the mesh in the X direction on the plane of integration
size_YT=floor(D_YT/h.YT)+1;




P_T=zeros(size_YT,size_XT); %Initialize the coordinate mesh on the surface of the transducer->Then fill in the position of each mesh grid
XcoorT=repmat(linspace(-D_XT/2,D_XT/2,size_XT),length(P_T(:,1)),1);% Setting the coordinate of the transducer
YcoorT=repmat(linspace(D_YT/2,-D_YT/2,size_YT).',1,length(P_T(1,:)));
for ii=1:NT
%     P_T=P_T+exp(-sqrt(-1)*K*sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-ii)-F_offset_y)^2))*(abs(XcoorT)<=xspan_T/2).*((YcoorT>=(Gap/2+(Gap+yspan_T)*(NT/2-ii))).*(YcoorT<=(Gap/2+yspan_T+(Gap+yspan_T)*(NT/2-ii)))).*...
%         sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-ii)-F_offset_y)^2)/sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-NT)-F_offset_y)^2);
    
      P_T=P_T+exp(-sqrt(-1)*K*sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-ii)-F_offset_y)^2))*(abs(XcoorT)<=xspan_T/2).*((YcoorT>=(Gap/2+(Gap+yspan_T)*(NT/2-ii))).*(YcoorT<=(Gap/2+yspan_T+(Gap+yspan_T)*(NT/2-ii))));%.*...
        %sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-ii)-F_offset_y)^2)/sqrt(F^2+F_offset_x^2+((Gap+yspan_T)*(0.5+NT/2-NT)-F_offset_y)^2);
end

figure
imagesc(1e3*linspace(-D_XT/2,D_XT/2,size_XT),1e3*linspace(D_YT/2,-D_YT/2,size_YT),abs(P_T));
xlabel('X [mm]');
ylabel('Y [mm]');
colormap('jet');
colorbar
%axis square



size_X1=floor(D_X1/h.X1);     %Size of the mesh in the X direction on the transducer
size_Y1=floor(D_Y1/h.Y1);



P_1 = zeros(size_Y1,size_X1); %Initialize the pressure pattern on phase plate
Count = 0;                  %Initialize counter

%% This part of the code is to increase the computation efficiency of rayleigh integral P_T->P_1;
Xcoor1=repmat(linspace(-D_X1/2,D_X1/2,size_X1),length(P_1(:,1)),1);
Ycoor1=repmat(linspace(D_Y1/2,-D_Y1/2,size_Y1).',1,length(P_1(1,:)));

S_T=(D_XT*D_YT)/(size_XT*size_YT);


P_r=zeros(1,length(r));
Count=0;
for i=1:length(P_r(:,1)) %Y index
    tic
    for j=1:length(P_r(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-0).^2+(YcoorT-0).^2+(r(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_r(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_r(:,1),time,W_bar)
end
P_Tn=zeros(size_YT,size_XT);
P_1n = zeros(size_Y1,size_X1); %Initialize the pressure pattern on phase plate

for ii=1:NT
    P_Tn=P_Tn+1*(abs(XcoorT)<=xspan_T/2).*((YcoorT>=Gap/2+(Gap+yspan_T)*(NT/2-ii)).*(YcoorT<=Gap/2+yspan_T+(Gap+yspan_T)*(NT/2-ii)));
end

Xcoor1=repmat(linspace(-D_X1/2,D_X1/2,size_X1),length(P_1(:,1)),1);
Ycoor1=repmat(linspace(D_Y1/2,-D_Y1/2,size_Y1).',1,length(P_1(1,:)));

S_T=(D_XT*D_YT)/(size_XT*size_YT);
P_rn=zeros(1,length(r));
Count=0;
for i=1:length(P_rn(:,1)) %Y index
    tic
    for j=1:length(P_rn(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-0).^2+(YcoorT-0).^2+(r(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_Tn.*(10.^(-loss*1e2*R_1/20));
        
        P_rn(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_rn(:,1),time,W_bar)
end

figure
%plot(r*1e3,abs(P_rn))
%hold on
plot(r*1e3,abs(P_r))
hold on

hold off
%legend('without phasing','with phasing');
legend('without phasing','F=5mm','F=10mm','F=20mm','F=40mm');
xlabel('r [mm]');
ylabel('On axis pressure [pa]');
pedit
title(strcat('On axis pressure with focal point at',{' '},num2str(F*1e3),'mm'));

for i=1:length(P_1(:,1)) %Y index
    tic
    for j=1:length(P_1(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-Xcoor1(i,j)).^2+(YcoorT-Ycoor1(i,j)).^2+(z0).^2);
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
        %Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_1(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_1(:,1),time,W_bar)
end
% for uz calculation in the xy plane
P_1_dz=zeros(size_Y1,size_X1);
clear R_1 Q
dz=lambda/10;
for i=1:length(P_1(:,1)) %Y index
    tic
    for j=1:length(P_1(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-Xcoor1(i,j)).^2+(YcoorT-Ycoor1(i,j)).^2+(z0+dz).^2);
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
        %Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_1_dz(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_1_dz(:,1),time,W_bar)
end

P_1_dz2=zeros(size_Y1,size_X1);
clear R_1 Q
dz=lambda/10;
for i=1:length(P_1(:,1)) %Y index
    tic
    for j=1:length(P_1(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-Xcoor1(i,j)).^2+(YcoorT-Ycoor1(i,j)).^2+(z0-dz).^2);
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
        %Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_1_dz2(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_1_dz2(:,1),time,W_bar)
end
uxy_z=-j/(rho0*2*pi*f)*((P_1_dz-P_1)/dz+(P_1-P_1_dz2)/dz)/2;
% for uy calculation in the xy plane
P_1_y1=zeros(length(P_1(:,1))+1,length(P_1(1,:)));
P_1_y2=P_1_y1;% to increase the accuracy by calculating the average of the left gradient and right gradient resulting in N array becoming N-2
P_1_y1(1:end-1,:)=P_1;
P_1_y2(2:end,:)=P_1;
P_1_dy1=(P_1_y2-P_1_y1)/(r_y1);
P_1_y3=zeros(length(P_1(:,1))-1,length(P_1(1,:)));
P_1_y4=P_1_y3;
P_1_y3=P_1(1:end-1,:);
P_1_y4(2:end,:)=P_1(1:end-2,:);
P_1_dy2=(P_1_y4-P_1_y3)/r_y1;
P_1_dyave=zeros(length(P_1(:,1)),length(P_1(1,:)));
P_1_dyave(2:end-1,:)=(P_1_dy2(2:end,:)+P_1_dy1(3:end-1,:))/2;
P_1_dyave(1,:)=P_1_dyave(2,:);
P_1_dyave(end,:)=P_1_dyave(end-1,:);
uxy_y=-j/(rho0*2*pi*f)*P_1_dyave;
% for ux calculation in the xy plane
P_1_x1=zeros(length(P_1(:,1)),length(P_1(1,:))+1);
P_1_x2=P_1_x1;% to increase the accuracy by calculating the average of the left gradient and right gradient resulting in N array becoming N-2
P_1_x1(:,2:end)=P_1;
P_1_x2(:,1:end-1)=P_1;
P_1_dx1=(P_1_x2-P_1_x1)/(r_x1);
P_1_x3=zeros(length(P_1(:,1)),length(P_1(1,:))-1);
P_1_x4=P_1_x3;
P_1_x3=P_1(:,2:end);
P_1_x4(:,1:end-1)=P_1(:,3:end);
P_1_dx2=(P_1_x4-P_1_x3)/r_x1;
P_1_dxave=zeros(length(P_1(:,1)),length(P_1(1,:)));
P_1_dxave(:,2:end-1)=(P_1_dx2(:,1:end-1)+P_1_dx1(:,2:end-2))/2;
P_1_dxave(:,1)=P_1_dxave(:,2);
P_1_dxave(:,end)=P_1_dxave(:,end-1);
uxy_x=-j/(rho0*2*pi*f)*P_1_dxave;

Ixy_x=1/2*real(P_1.*conj(uxy_x));
Ixy_y=1/2*real(P_1.*conj(uxy_y));
Ixy_z=1/2*real(P_1.*conj(uxy_z));

Ixy=abs(Ixy_x+Ixy_y+Ixy_z);

figure
subplot(211)
imagesc(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),1e-4*1/2*abs(P_1).^2/(rho0*c0)) 
xlabel('X [mm]')
ylabel('Y [mm]');
colormap('jet')
colorbar
title('Isppa(W/cm^2) approx Zo');
subplot(212)
imagesc(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),1e-4*Ixy) 
xlabel('X [mm]')
ylabel('Y [mm]');
colormap('jet')
colorbar
title('Isppa(W/cm^2) use uo');


clear R_1 Q

figure
imagesc(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1)) 
xlabel('X [mm]')
ylabel('Y [mm]');
colormap('jet')
colorbar
%axis square
title(strcat('Magnitude at focal point',{' '},num2str(F*1e3),'mm'));
figure
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1(:,round(length(P_1(1,:))/2))));
xlabel('Y [mm]')
ylabel('Pressure [Pa]');

disp(strcat('surface pressure is',{' '},num2str(rho0*c0*u0),'pa, ','focal point pressure is',{' '},num2str(interp1(r*1e3,abs(P_r),F*1e3))))
figure
subplot(1,2,1)
h = surf(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1)/(rho0*c0*u0));
h.EdgeAlpha=0; %hides mesh
title('magnitude at focal point','Interpreter','latex');
axis square
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0[mm]$','interpreter','latex');
ylabel('$y_0[mm]$','interpreter','latex');
view(2)
colorbar
colormap('jet')

subplot(1,2,2)
g = surf(abs(angle(P_1)));
g.EdgeAlpha=0; %hides mesh
%axis square
%axis off
title('phase incident on phase plate','Interpreter','latex')
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)
colorbar
colormap('jet')
clear R_1 Q
P_yz=zeros(round(D_Y1/r_y1),length(r));
Count=0;
Zcoor_yz=repmat(r,length(P_yz(:,1)),1);
Ycoor_yz=repmat(linspace(D_Y1/2,-D_Y1/2,length(P_yz(:,1))).',1,length(P_yz(1,:)));
for i=1:length(P_yz(:,1)) %Z index
    tic
    for j=1:length(P_yz(1,:)) %Y index
        Count=Count+1;
        
        R_1=sqrt(XcoorT.^2+(Ycoor_yz(i,j)-YcoorT).^2+(Zcoor_yz(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_yz(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_yz(:,1),time,W_bar)
end


% % for uy calculation in the xy plane
% P_yz_y1=zeros(length(P_yz(:,1))+1,length(P_yz(1,:)));
% P_yz_y2=P_yz_y1;% to increase the accuracy by calculating the average of the left gradient and right gradient resulting in N array becoming N-2
% P_yz_y1(2:end,:)=P_yz;
% P_yz_y2(1:end-1,:)=P_yz;
% P_yz_dy1=(P_yz_y2-P_yz_y1)/(r_y1);
% P_yz_y3=zeros(length(P_yz(:,1))-1,length(P_yz(1,:)));
% P_yz_y4=P_yz_y3;
% P_yz_y3=P_yz(2:end,:);
% P_yz_y4(1:end-1,:)=P_yz(3:end,:);
% P_yz_dy2=(P_yz_y4-P_yz_y3)/r_y1;
% P_yz_dyave=zeros(length(P_yz(:,1)),length(P_yz(1,:)));
% P_yz_dyave(2:end-1,:)=(P_yz_dy2(1:end-1,:)+P_yz_dy1(2:end-2,:))/2;
% P_yz_dyave(1,:)=P_yz_dyave(2,:);
% P_yz_dyave(end,:)=P_yz_dyave(end-1,:);
% uyz_y=-j/(rho0*2*pi*f)*P_yz_dyave;
% % for ux calculation in the xy plane
% P_yz_x1=zeros(length(P_yz(:,1)),length(P_yz(1,:))+1);
% P_yz_x2=P_yz_x1;% to increase the accuracy by calculating the average of the left gradient and right gradient resulting in N array becoming N-2
% P_yz_x1(:,2:end)=P_yz;
% P_yz_x2(:,1:end-1)=P_yz;
% P_yz_dx1=(P_yz_x2-P_yz_x1)/(r_x1);
% P_yz_x3=zeros(length(P_yz(:,1)),length(P_yz(1,:))-1);
% P_yz_x4=P_yz_x3;
% P_yz_x3=P_yz(:,2:end);
% P_yz_x4(:,1:end-1)=P_yz(:,3:end);
% P_yz_dx2=(P_yz_x4-P_yz_x3)/r_x1;
% P_yz_dxave=zeros(length(P_yz(:,1)),length(P_yz(1,:)));
% P_yz_dxave(:,2:end-1)=(P_yz_dx2(:,1:end-1)+P_yz_dx1(:,2:end-2))/2;
% P_yz_dxave(:,1)=P_yz_dxave(:,2);
% P_yz_dxave(:,end)=P_yz_dxave(:,end-1);
% uyz_x=-j/(rho0*2*pi*f)*P_yz_dxave;

%% no focusing

P_Tn=zeros(size_YT,size_XT);
P_1n = zeros(size_Y1,size_X1); %Initialize the pressure pattern on phase plate

for ii=1:NT
    P_Tn=P_Tn+1*(abs(XcoorT)<=xspan_T/2).*((YcoorT>=Gap/2+(Gap+yspan_T)*(NT/2-ii)).*(YcoorT<=Gap/2+yspan_T+(Gap+yspan_T)*(NT/2-ii)));
end

Xcoor1=repmat(linspace(-D_X1/2,D_X1/2,size_X1),length(P_1(:,1)),1);
Ycoor1=repmat(linspace(D_Y1/2,-D_Y1/2,size_Y1).',1,length(P_1(1,:)));

S_T=(D_XT*D_YT)/(size_XT*size_YT);


for i=1:length(P_1n(:,1)) %Y index
    tic
    for j=1:length(P_1n(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-Xcoor1(i,j)).^2+(YcoorT-Ycoor1(i,j)).^2+(z0).^2);
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
        %Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_Tn.*(10.^(-loss*1e2*R_1/20));
        
        P_1n(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_1n(:,1),time,W_bar)
end

clear R_1 Q
%r=linspace(0,10e-2,4e3);
P_rn=zeros(1,length(r));
Count=0;
for i=1:length(P_rn(:,1)) %Y index
    tic
    for j=1:length(P_rn(1,:)) %X index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-0).^2+(YcoorT-0).^2+(r(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_Tn.*(10.^(-loss*1e2*R_1/20));
        
        P_rn(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_rn(:,1),time,W_bar)
end



clear R_1 Q
P_yzn=zeros(round(D_Y1/(1/50*c0/f)),length(r));
Count=0;
Zcoor_yz=repmat(r,length(P_yzn(:,1)),1);
Ycoor_yz=repmat(linspace(D_Y1/2,-D_Y1/2,length(P_yzn(:,1))).',1,length(P_yzn(1,:)));
for i=1:length(P_yzn(:,1)) %Z index
    tic
    for j=1:length(P_yzn(1,:)) %Y index
        Count=Count+1;
        
        R_1=sqrt(XcoorT.^2+(Ycoor_yz(i,j)-YcoorT).^2+(Zcoor_yz(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_Tn.*(10.^(-loss*1e2*R_1/20));
        
        P_yzn(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_yzn(:,1),time,W_bar)
end

clear R_1 Q
P_xz=zeros(round(D_X1/(0.25*c0/f)),length(r));
Count=0;
Zcoor_xz=repmat(r,length(P_xz(:,1)),1);
Xcoor_xz=repmat(linspace(D_X1/2,-D_X1/2,length(P_xz(:,1))).',1,length(P_xz(1,:)));
for i=1:length(P_xz(:,1)) %X index
    tic
    for j=1:length(P_xz(1,:)) %Z index
        Count=Count+1;
        
        R_1=sqrt((XcoorT-Xcoor_xz(i,j)).^2+(YcoorT).^2+(Zcoor_xz(i,j).^2));
        Q=-sqrt(-1)*rho0*K*c0*u0*S_T*(exp(sqrt(-1)*K*R_1)./(2*pi*R_1));
       % Q=Q.*exp(-(XcoorT.^2 + YcoorT.^2)/((a)^2));
        Q=Q.*P_T.*(10.^(-loss*1e2*R_1/20));
        
        P_xz(i,j) = sum(Q,'all');%exp(-(XcoorT.^2 + YcoorT.^2)/((a/5)^2)).*P_T.*S_T*exp(-sqrt(-1)*K*R_1)./(2*pi*R_1)),'all');
    end
    time=toc;
        waitbar_show(i,P_xz(:,1),time,W_bar)
end
%% Plots are generated here
figure
imagesc(1e3*r,1e3*linspace(D_X1/2,-D_X1/2,length(P_xz(:,1))),(abs(P_xz)) )
colormap('jet')
colorbar
%axis square
xlabel('Z [mm]')
ylabel('X [mm]');
title(strcat('X-Z Magnitude at focal point ',{' '},num2str(F*1e3),'mm'));

figure
imagesc(1e3*r,1e3*linspace(D_X1/2,-D_X1/2,length(P_xz(:,1))),(abs(P_xz)) )
colormap('jet')
colorbar
%axis square
xlabel('Z [mm]')
ylabel('X [mm]');
title(strcat('X-Z Magnitude at focal point ',{' '},num2str(F*1e3),'mm'));


figure

imagesc(1e3*r(1:end),1e3*linspace(D_Y1/2,-D_Y1/2,length(P_yzn(:,1))),(abs(P_yzn(:,1:end))))
colormap('jet')
colorbar
%axis square
xlabel('Z [mm]')
ylabel('Y [mm]');
title(strcat('Y-Z Magnitude at focal point (no focusing) '));

figure

imagesc(1e3*r(1:end),1e3*linspace(D_Y1/2,-D_Y1/2,length(P_yz(:,1))),(abs(P_yz(:,1:end))))
colormap('jet')
colorbar
%axis square
xlabel('Z [mm]')
ylabel('Y [mm]');
title(strcat('Y-Z Magnitude at focal point ',{' '},num2str(F*1e3),'mm'));

figure
imagesc(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1n)) 
colormap('jet')
colorbar
%axis square
xlabel('[mm]')
ylabel('[mm]');
title(strcat('Magnitude at focal point '));
figure
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1n(:,round(length(P_1n(1,:))/2))));
xlabel('Y [mm]')
ylabel('Pressure [Pa] (no focusing)');
figure
plot(r*1e3,abs(P_rn))
xlabel('r [mm]');
ylabel('On axis pressure [pa]');
pedit
title(strcat('On axis pressure (no focusing)'));

figure
subplot(1,2,1)
h = surf(1e3*linspace(-D_X1/2,D_X1/2,size_X1),1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1n)/(rho0*c0*u0));
h.EdgeAlpha=0; %hides mesh
title('magnitude at focal point(no focusing)','Interpreter','latex');
%axis square
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0[mm]$','interpreter','latex');
ylabel('$y_0[mm]$','interpreter','latex');
view(2)
colorbar
colormap('jet')

subplot(1,2,2)
g = surf(abs(angle(P_1n)));
g.EdgeAlpha=0; %hides mesh
%axis square
%axis off
title('phase incident on phase plate(no focusing)','Interpreter','latex')
set(gca,'xtick',[],'ytick',[])
xlabel('$x_0$','interpreter','latex');
ylabel('$y_0$','interpreter','latex');
view(2)
colorbar
colormap('jet')

figure
plot(r*1e3,abs(P_r))
xlabel('r [mm]');
ylabel('On axis pressure [pa]');
hold on
pedit
title(strcat('On axis pressure with focal point at ',{' '},num2str(F*1e3),'mm'));
plot(r*1e3,abs(P_rn))
hold off
legend('With focusing','Without focusing');

figure
subplot(121)
P_1ny=abs((abs(P_1n(:,round(length(P_1n(1,:))/2)))));
P_1nx=abs((abs(P_1n(round(length(P_1n(:,1))/2),:))));
P_1y=abs((abs(P_1(:,round(length(P_1(1,:))/2)))));
P_1x=abs((abs(P_1(round(length(P_1(:,1))/2),:))));
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1n(:,round(length(P_1n(1,:))/2))));
xlabel('y [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),abs(P_1(:,round(length(P_1(1,:))/2))));
hold off
pedit
legend('No focusing','With focusing');
subplot(122)

plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1),abs(P_1n(round(length(P_1n(:,1))/2),:)));
xlabel('x [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1),abs(P_1(round(length(P_1(:,1))/2),:)));
hold off
pedit
legend('No focusing','With focusing');
P_1y_filled=interp1(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),P_1y,1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1*10));
P_1ny_filled=interp1(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),P_1ny,1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1*10));
P_1x_filled=interp1(1e3*linspace(D_X1/2,-D_X1/2,size_X1),P_1x,1e3*linspace(D_X1/2,-D_X1/2,size_X1*10));
P_1nx_filled=interp1(1e3*linspace(D_X1/2,-D_X1/2,size_X1),P_1nx,1e3*linspace(D_X1/2,-D_X1/2,size_X1*10));
figure
subplot(121)
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1*10),P_1ny_filled);
xlabel('y [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1*10),P_1y_filled);
hold off
pedit
legend('No focusing','With focusing');
subplot(122)
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1*10),P_1nx_filled);
xlabel('x [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1*10),P_1x_filled);
hold off
pedit
legend('No focusing','With focusing');


figure
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1*10),0.5e-4*abs(P_1y_filled).^2/(rho0*c0));
hold on
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1*10),0.5e-4*abs(P_1x_filled).^2/(rho0*c0));
hold off
pedit
set(gca, 'FontSize', 26);
xlabel('Position [mm]', 'FontSize', 26);
ylabel('Isppa [W/cm^2]', 'FontSize', 26);
legend('Y crosssection(X=0)','X crosssection(Y=0)');




figure
subplot(121)
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),P_1ny);
xlabel('y [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(D_Y1/2,-D_Y1/2,size_Y1),P_1y);
hold off
pedit
legend('No focusing','With focusing');
subplot(122)
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1),P_1nx);
xlabel('x [mm]')
ylabel('Pressure [Pa]');
hold on 
plot(1e3*linspace(-D_X1/2,D_X1/2,size_X1),P_1x);
hold off
pedit
legend('No focusing','With focusing');
title('Hilbert Transform');
figure
subplot(211)
imagesc(1e3*linspace(D_X1/2,-D_X1/2,size_X1),1e3*linspace(-D_Y1/2,D_Y1/2,size_Y1),abs(P_1n)) 
colormap('jet')
colorbar
%axis square
xlabel('X [mm]')
ylabel('Y [mm]');
title(strcat('Magnitude at focal point(no focus)'));
subplot(212)
imagesc(1e3*linspace(D_X1/2,-D_X1/2,size_X1),1e3*linspace(-D_Y1/2,D_Y1/2,size_Y1),abs(P_1)) 
colormap('jet')
colorbar
%axis square
xlabel('X [mm]')
ylabel('Y [mm]');
title(strcat('Magnitude at focal point ',{' '},num2str(F*1e3),'mm'));
filename=strcat('D:\NIH brainprobe\','DATA_Gap',num2str(Gap*1e3),'xspan',num2str(xspan_T*1e3),'yspan',num2str(yspan_T*1e3),'Number',num2str(NT),'mm.mat');
%save(filename);