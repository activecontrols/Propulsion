% This function contains the core characteristic equation calculations and
% plotting code.
function mainfunc(hObject, eventdata, s1)
global a b cstar dcdMR m mf mo MR f fmax pc pod pfd
global tau_f tau_o theta_g xmax zmax Zf Zo df
global LoO CmO LoF CmF
global h2
format('long')
% Initialize variables (These must be initialized mainly to carry values
% through to the save function, since the save function can use either 2 of
% these parameters or 4 of these parameters)
dpf_pc=0;
dpo_pc=0;
dpf_pc1=0;
dpf_pc2=0;
dpo_pc1=0;
dpo_pc2=0;
% get screensize and define figure layout to better layout figures
scrsz=get(0,'ScreenSize'); % get screen size
figvert=420; % vertical figure dimension
bordwidth=10; % figure border width
titlebar=77; % title bar width
spw=20; % space width
%%% Numbers to Strings for putting text on plots later
str1=num2str(pfd/pc);
str2=num2str(pod/pc);
str3=num2str(tau_o);
str4=num2str(tau_f);
str6=num2str(fmax); %str5 is defined below
%%% General Characteristic Equation
omega=2*pi*f; % Circular Frequency
s=i*omega; % Stab. space curve is when Re. part of s=a+i*omega is 0
% Impedance for plotting Nyquist
if get(h2,'Value')==1 % If Feed System On
    % Real and imaginary parts from manifold impedance formulation
    RRo=real(-LoO*s-1./(CmO*s));
    RRf=real(-LoF*s-1./(CmF*s));
    Io=imag(-LoO*s-1./(CmO*s));
    Ifu=imag(-LoF*s-1./(CmF*s));
    % impedance formulation
    Zo=-pod/mo/a+RRo+i*Io;
    Zf=-pfd/mf/b+RRf+i*Ifu;
else
    Zo=-pod/mo/a;
    Zf=-pfd/mf/b;
end
X=(1+(1+MR)*dcdMR/cstar)*(pc/m);
Y=(1-MR*(1+MR)*dcdMR/cstar)*(pc/m);
H=1-1./(1+theta_g*s).*( X./Zo.*exp(-s*tau_o) + Y./Zf.*exp(-s*tau_f) );
C=H-1;
%%%------------------------------------------------------------------------
%%% This Section Plots the Nyquist Diagram
%%%------------------------------------------------------------------------
%%% Nyquist Stability Diagram for H(s)-1 plane
% % Nyquist Plot
xm=2; % x-axis maximum range
% figure border width is 10 and title bar, etc is 77
figure('Name','Nyquist Diagram','NumberTitle','off','Position',...
    [bordwidth+spw scrsz(4)-figvert-titlebar 4*figvert/3 figvert])
% Plot in H(s)-1 plane
plot(real(C),imag(C),'LineWidth',2)
title('Nyquist Diagram, H(s) - 1 plane')
xlabel('Real Axis')
ylabel('Imaginary Axis')
axis equal
axis([-xm xm -xm xm])
% Hold so Nyquist Diagram is plotted with Polar Grid
hold on
% % Include Stability Origin Point
plot([-1],[0],'bo','MarkerFaceColor','b')
% % Include Specified DP/Pc Information written on plot
text(-5/6*xm,5/6*xm, ['\DeltaP_f/P_c= ',str1], 'HorizontalAlignment', 'left')
text(-5/6*xm,2.1/3*xm, ['\DeltaP_o/P_c= ',str2], 'HorizontalAlignment', 'left')
text(-5/6*xm+xm,5/6*xm, ['\tau_o (sec) = ',str3], 'HorizontalAlignment', 'left')
text(-5/6*xm+xm,2.1/3*xm, ['\tau_f (sec) = ',str4], 'HorizontalAlignment', 'left')
text(-5/6*xm,-5/6*xm, ['Freq. Range (Hz) = -',str6, ' to ', str6], 'HorizontalAlignment', 'left')
% % Polar Grid Defined with Next Series of Plots
xr=(-xm:xm/100:xm); % Full x range
% % Polar Grid Defined with Next Series of Plots
xr=(-xm:xm/100:xm); % Full x range
% Radial Circles
plot(xr,real(sqrt((xm*1/3)^2-xr.^2)),':','Color',[.6 .6 .6])
plot(xr,real(-sqrt((xm*1/3)^2-xr.^2)),':','Color',[.6 .6 .6])
plot(xr,real(sqrt((xm*2/3)^2-xr.^2)),':','Color',[.6 .6 .6])
plot(xr,real(-sqrt((xm*2/3)^2-xr.^2)),':','Color',[.6 .6 .6])
plot(xr,real(sqrt(xm^2-xr.^2)),':','Color',[.6 .6 .6])
plot(xr,real(-sqrt(xm^2-xr.^2)),':','Color',[.6 .6 .6])
% Spokes
plot(sqrt(3)/2*xr,sqrt(3)/2*xr*sqrt(3)/3,':','Color',[.6 .6 .6])
plot(0.5*xr,0.5*xr*sqrt(3),':','Color',[.6 .6 .6])
plot(sqrt(3)/2*xr,-sqrt(3)/2*xr*sqrt(3)/3,':','Color',[.6 .6 .6])
plot(0.5*xr,-0.5*xr*sqrt(3),':','Color',[.6 .6 .6])
line([0,0],[-xm,xm],'LineStyle',':','Color',[.6 .6 .6])
line([-xm,xm],[0,0],'LineStyle',':','Color',[.6 .6 .6])
% Unit Circle
plot(xr/xm,sqrt(1-(xr/xm).^2),'k','LineWidth',2)
plot(xr/xm,-sqrt(1-(xr/xm).^2),'k','LineWidth',2)
hold off
% Redefine Frequency to only the positive range for the rest of the analysis
f=0:df:fmax;
omega=2*pi*f; % Circular Frequency
s=i*omega; % Stab. space curve is when Re. part of s=a+i*omega is 0
%%%--------------------------------------------------------------------
%%% This Section is the Main Calcaultioon
%%%--------------------------------------------------------------------
%%% Decide which analysis method to use
if get(h2,'Value')==1 % Feed System On
    % Redefine real and imag parts for only the positive range of
    % frequencies for the rest of the analysis
    RRo=real(-LoO*s-1./(CmO*s));
    RRf=real(-LoF*s-1./(CmF*s));
    Io=imag(-LoO*s-1./(CmO*s));
    Ifu=imag(-LoF*s-1./(CmF*s));
    % Obtain DPf/Pc as a function of DPo/Pc and Frequency
    % Real parts of Zo and Zf are solved from the characteristic equation
    % The real parts can be related back to Dpf/Pc and Dpo/Pc
    A=(2*pi*theta_g*f.*cos(2*pi*f*tau_f)+sin(2*pi*f*tau_f))...
        +Ifu.*(1+(2*pi*f*theta_g).^2)/Y;
    B=X*sin(2*pi*f*(tau_o-tau_f))...
        +2*X/Y*Ifu.*(2*pi*theta_g*f.*sin(2*pi*f*tau_o)-cos(2*pi*f*tau_o));
    C=Io.^2.*(2*pi*theta_g*f.*cos(2*pi*f*tau_f)+sin(2*pi*f*tau_f))...
        +X*Io.*cos(2*pi*f*(tau_o-tau_f))...
        +2*X/Y*Io.*Ifu.*(2*pi*theta_g*f.*cos(2*pi*f*tau_o)+sin(2*pi*f*tau_o))...
        +Ifu.*Io.^2.*(1+(2*pi*f*theta_g).^2)/Y+X^2/Y*Ifu;
    % Real part of Zo
    Ro1=(-B+sqrt(B.^2-4*A.*C))./(2*A);
    Ro2=(-B-sqrt(B.^2-4*A.*C))./(2*A);
    % Real part of Zf
    Rf1=(-Y*Ro1.*sin(2*pi*f*tau_f)+Y*Io.*cos(2*pi*f*tau_f)...
        +X*Ifu.*cos(2*pi*f*tau_o)+2*pi*theta_g*f.*Io.*Ifu-Ifu.*Ro1)...
        ./(Ro1.*f*2*pi*theta_g+Io+X*sin(2*pi*f*tau_o));
        Rf2=(-Y*Ro2.*sin(2*pi*f*tau_f)+Y*Io.*cos(2*pi*f*tau_f)...
        +X*Ifu.*cos(2*pi*f*tau_o)+2*pi*theta_g*f.*Io.*Ifu-Ifu.*Ro2)...
        ./(Ro2.*f*2*pi*theta_g+Io+X*sin(2*pi*f*tau_o));
    % Relate Ro to Dpo/Pc
    dpo_pc1=(RRo*mo*a/pc-Ro1*mo*a/pc);
    dpo_pc2=(RRo*mo*a/pc-Ro2*mo*a/pc);
    dpf_pc1=(RRf*mf*b/pc-Rf1*mf*b/pc);
    dpf_pc2=(RRf*mf*b/pc-Rf2*mf*b/pc);
    % Elliminate Values that are less than 0 since Matlab still Plots
    % Negative Values when a Positive Range is Defined. This clutters the
    % plot. Convert Negatives to NaN.
    [row,col]=find(dpf_pc1<0);
    %negs=[ ];
    negs=zeros(size(row));
    for k=1:length(row)
    negs(k)=[dpf_pc1(row(k),col(k))];
    dpf_pc1(row(k),col(k))=NaN;
    end
    [row,col]=find(dpo_pc1<0);
    negs=zeros(size(row));
    for k=1:length(row)
    negs(k)=[dpo_pc1(row(k),col(k))];
    dpo_pc1(row(k),col(k))=NaN;
    end
    [row,col]=find(dpf_pc2<0);
    negs=zeros(size(row));
    for k=1:length(row)
        negs(k)=[dpf_pc2(row(k),col(k))];
        dpf_pc2(row(k),col(k))=NaN;
    end
    [row,col]=find(dpo_pc2<0);
    negs=zeros(size(row));
    for k=1:length(row)
        negs(k)=[dpo_pc2(row(k),col(k))];
        dpo_pc2(row(k),col(k))=NaN;
    end
    % Eliminate Values that are complex since Matlab still Plots
    % Real Part. Convert complex numbers to NaN.
    [row,col]=find(imag(dpf_pc1)~=0);
    comp=zeros(size(row));
    for k=1:length(row)
        comp(k)=[dpf_pc1(row(k),col(k))];
        dpf_pc1(row(k),col(k))=NaN;
    end
    [row,col]=find(imag(dpo_pc1)~=0);
    comp=zeros(size(row));
    for k=1:length(row)
        comp(k)=[dpo_pc1(row(k),col(k))];
        dpo_pc1(row(k),col(k))=NaN;
    end
    [row,col]=find(imag(dpf_pc2)~=0);
    comp=zeros(size(row));
    for k=1:length(row)
        comp(k)=[dpf_pc2(row(k),col(k))];
        dpf_pc2(row(k),col(k))=NaN;
    end
    [row,col]=find(imag(dpo_pc2)~=0);
    comp=zeros(size(row));
    for k=1:length(row)
        comp(k)=[dpo_pc2(row(k),col(k))];
        dpo_pc2(row(k),col(k))=NaN;
    end
    %%%--------------------------------------------------------------------
    %%% This Section Plots the Stability SpaceCurve for Feed System On
    %%%--------------------------------------------------------------------
    % % Plot the spaceCurve Dpf/Pc(Dpo/pc,f)
    figure('Name','Space Curve','NumberTitle','off','Position',...
        [scrsz(3)-4*figvert/3-spw scrsz(4)-figvert-titlebar 4*figvert/3 figvert])
    F1=plot3(dpo_pc1,f,dpf_pc1,'b-',dpo_pc2,f,dpf_pc2,'b-');
    set(F1,'LineWidth',2)
    axis([0 xmax 0 fmax 0 zmax]);
    view([15,30]);
    title('Space Curve')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,fmax,0.8*zmax,...
        ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    text(0.6*xmax,fmax,0.65*zmax,...
        ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    %%%--------------------------------------------------------------------
    %%% This Section Plots the Wenzel and Szuch Figures for Feed System On
    %%%--------------------------------------------------------------------
    % % Plot the curve Dpf/Pc vs Dpo/pc
    figure('Name','Pressure Drop - Stability Plot','NumberTitle','off','Position',...
        [bordwidth+spw bordwidth+spw 4*figvert/3 figvert])
    hold on
    F2a=plot3(dpo_pc1,f,dpf_pc1,'b-',dpo_pc2,f,dpf_pc2,'b-');
    F2b=plot3(pod/pc,0,pfd/pc,'bo','MarkerFaceColor','b');
    hold off
    set(F2a,'LineWidth',2)
    axis([0 xmax 0 fmax 0 zmax]);
    view([0,0]);
    title('Pressure Drop - Stability Plot')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,0,0.85*zmax,...
        ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    text(0.6*xmax,0,0.75*zmax,...
        ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    % % Plot the curve Dpf/Pc vs f
    figure('Name','Frequency at Neutral Stability','NumberTitle','off','Position',...
        [scrsz(3)-4*figvert/3-spw bordwidth+spw 4*figvert/3 figvert])
        F3=plot3(dpo_pc1,f,dpf_pc1,'b-',dpo_pc2,f,dpf_pc2,'b-');
    set(F3,'LineWidth',2)
    axis([0 xmax 0 fmax 0 zmax]);
    view([0,90]);
    title('Frequency at Neutral Stability')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,0.85*fmax,0,...
        ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    text(0.6*xmax,0.75*fmax,0,...
        ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
savefile(dpf_pc,dpo_pc,dpf_pc1,dpf_pc2,dpo_pc1,dpo_pc2,f,h2)
else
    % % Obtain DPf/Pc as a function of DPo/Pc and Frequency
    % % Need to put constants (X and Y) in correct form
    % % Equation in TN D-3080 is in a slightly different form
    % (MR/(MR+1)+MR/cstar*dcdMR)*(pc/mo) = (1+(1+MR)/cstar*dcdMR)*(pc/m)
    % (1/(MR+1)-MR/cstar*dcdMR)*(pc/mf) = (1-MR*(1+MR)/cstar*dcdMR)*(pc/m)
    X=X*a*mo/pc; % Include oxidizer exponent also in X
    Y=Y*b*mf/pc; % Include fuel exponent also in Y
    % % Obtain DPf/Pc as a function of DPo/Pc and Frequency
    % Also consider case when tau_f = tau_o; general solution becomes
    % trivial so the single time lag model is used in the conditional
    % statement. Single time lag model uses a solver to obtain frequency.
    if tau_o==tau_f % single time lag model (Equations derived in Maple)
        tau=tau_o;
        dpo_pc=0:xmax/100:xmax;
        
        % Newton Rhapson for frequency calc
        f_new=1/(4*tau);% initial guess shown as a bound from Maple
        eps=1; % initialize
        count=1; % initialize

        while (eps>=0.00000001) && (count<100000)
            f_old=f_new;
            g_old=2*theta_g*pi*f_old*cos(2*pi*f_old*tau)...
                +sin(2*pi*f_old*tau);
            gp_old=(2*theta_g*pi+2*pi*tau)*cos(2*pi*f_old*tau)...
                -4*theta_g*pi^2*f_old*sin(2*pi*f_old*tau);
            f_new=f_old-g_old/gp_old;
            eps=abs((f_new-f_old)/f_old);
            count=count+1;
        end
        if count==100000
            disp('check convergence')
        end
        % Solved 2*theta_g*pi*f*cos(2*pi*f*tau)+sin(2*pi*f*tau)=0
        freqsolve=f_new;
        str5=num2str(freqsolve,4);
        f=freqsolve*ones(size(dpo_pc));
        dpf_pc=Y.*dpo_pc./(dpo_pc.*sqrt(1+(2*pi*f*theta_g).^2)-X);
    else % use the double time lag model
        dpo_pc=(sin(2*pi*f*(tau_o-tau_f))*X)./...
            (sin(2*pi*tau_f*f)+theta_g*2*pi*f.*cos(2*pi*tau_f*f));
        dpf_pc=Y*sin(2*pi*f*(tau_f-tau_o))./...
            (sin(2*pi*tau_o*f)+theta_g*2*pi*f.*cos(2*pi*f*tau_o));
    end
    % Remove Values that are Less than 0 since Matlab still Plots
    % Negative Values when a Positive Range is Defined. This clutters the
    % plot. Convert Negatives to NaN.
    [row,col]=find(dpf_pc<0);
    negs=zeros(size(row));
    for k=1:length(row)
    negs(k)=[dpf_pc(row(k),col(k))];
    dpf_pc(row(k),col(k))=NaN;
    end
    [row,col]=find(dpo_pc<0);
    negs=zeros(size(row));
    for k=1:length(row)
    negs(k)=[dpo_pc(row(k),col(k))];
    dpo_pc(row(k),col(k))=NaN;
    end
    %%%--------------------------------------------------------------------
    %%% This Section Plots the Stability SpaceCurve for Feed System Off
    %%%--------------------------------------------------------------------
    % % Plot the spaceCurve Dpf/Pc(Dpo/pc,f)
    figure('Name','Space Curve','NumberTitle','off','Position',...
        [scrsz(3)-4*figvert/3-spw scrsz(4)-figvert-titlebar 4*figvert/3 figvert])
    plot3(dpo_pc,f,dpf_pc,'LineWidth',2)
    axis([0 xmax 0 fmax 0 zmax]);
    view([15,30]);
    title('Space Curve')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,fmax,0.8*zmax,...
        ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    text(0.6*xmax,fmax,0.65*zmax,...
        ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    %%%--------------------------------------------------------------------
    %%% This Section Plots the Wenzel and Szuch Figures for Feed System Off
    %%%--------------------------------------------------------------------
    % % Plot the curve Dpf/Pc vs Dpo/pc
    figure('Name','Pressure Drop - Stability Plot','NumberTitle','off','Position',...
    [bordwidth+spw bordwidth+spw 4*figvert/3 figvert])
    hold on
    plot3(dpo_pc,f,dpf_pc,'LineWidth',2)
    plot3(pod/pc,0,pfd/pc,'bo','MarkerFaceColor','b')
    hold off
    axis([0 xmax 0 fmax 0 zmax]);
    view([0,0]);
    title('Pressure Drop - Stability Plot')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,0,0.85*zmax,...
        ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    text(0.6*xmax,0,0.75*zmax,...
        ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
        'FontSize', 16)
    % Show frequency on plot if single time lag model is used
    if (tau_o==tau_f)
        text(0.6*xmax,0,0.65*zmax,...
            ['freq= ',str5,' Hz'], 'HorizontalAlignment', 'left',...
            'FontSize', 16)
    end
    % % Plot the curve Dpf/Pc vs f
    figure('Name','Frequency at Neutral Stability','NumberTitle','off','Position',...
        [scrsz(3)-4*figvert/3-spw bordwidth+spw 4*figvert/3 figvert])
    plot3(dpo_pc,f,dpf_pc,'LineWidth',2)
    axis([0 xmax 0 fmax 0 zmax]);
    view([0,90]);
    title('Frequency at Neutral Stability')
    xlabel('\DeltaP_o / P_c')
    zlabel('\DeltaP_f / P_c')
    ylabel('Frequency (Hz)')
    grid on
    % Show Time lags on Plot
    text(0.6*xmax,0.85*fmax,0,...
    ['\tau_o = ',str3,' sec'], 'HorizontalAlignment', 'left',...
    'FontSize', 16)
    text(0.6*xmax,0.75*fmax,0,...
    ['\tau_f = ',str4,' sec'], 'HorizontalAlignment', 'left',...
    'FontSize', 16)
savefile(dpf_pc,dpo_pc,dpf_pc1,dpf_pc2,dpo_pc1,dpo_pc2,f,h2)
end