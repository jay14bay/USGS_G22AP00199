clear all;close all;clc

%% define the grid of stations X and Y, in km
SiteX=-80:1:180;
SiteY=-80:1:180;

%% define the rupture

% the rupture is described as required for use with GC2 (Spudich and Chiou, 2015)
% ftraces is the description of the rupture trace coordinates, strike angles, and segment lengths
% ftraces is a structure with length equal to the number of fault strands, where the strands do not need to be connected at their endpoints.
% for each strand i, ftraces(i).trace lists the X, Y coordinates of the segments which define the strand. this must be (n+1) by 2, where n is the number of segments in the ith strand.
% for each strand i, ftraces(i).strike is a 1xn vector of segment strike angles
% for each strand i, ftraces(i).l is a 1xn vector of segment lengths
% see further descriptions in the accompanying function GC2.m

% to illustrate, these three examples below all define the same rupture in a different way.  Two of them are commented out.

% Option 1: a single strand with a single segment 80 km in length
    clear ftraces
    ftraces(1).trace=[0 0;
                      0 80];
    ftraces(1).strike=[0];
    ftraces(1).l= [80];

% Option 2: two strands each with one segment 40 km in length
%     clear ftraces
%     ftraces(1).trace=[0 0;
%                       0 40];
%     ftraces(1).strike=[0];
%     ftraces(1).l= [40]; 
% 
%     ftraces(2).trace=[0 40;
%                       0 80];
%     ftraces(2).strike=[0];
%     ftraces(2).l= [40];

% Option 3: one strand with two segments, each 40 km in length
%     clear ftraces
%     ftraces(1).trace=[0 0;
%                       0 40
%                       0 80];
%     ftraces(1).strike=[0 0];
%     ftraces(1).l= [40 40]; 

nt=length(ftraces);

M=7.2; % moment magnitude
% determinine the model version. 1->simulation-based. 2->NGA-W2 data-based
Version=1; 
% select the period at which to show the effect
Tdo=3;

% characteristic rupture parameters
Rake=0; % rake in deg
Ztor=0; % Ztor, must be positive, in km

% specify the coordinates of the epicenter and GC2 origin, po
type.epi=[0 10]; % X, Y
type.po=[0 10]; % in this case, the same as the epicenter
    
%% call the Spudich and Chiou (2015) GC2 function
type.str='JB'; 
discordant=false;
gridflag=true;
[T,U,W,reference_axis,p_origin,nominal_strike,Upo]=GC2(ftraces,SiteX,SiteY,type,discordant,gridflag);

% calculate the maximum value of S in each direction for this hypocenter; it is U calculated at the nominal strike ends
[~,Uend,~,~,~,~,~,~]=GC2(ftraces,nominal_strike.a(1,1),nominal_strike.a(1,2),type,discordant,gridflag);
[~,Uend2,~,~,~,~,~,~]=GC2(ftraces,nominal_strike.a(2,1),nominal_strike.a(2,2),type,discordant,gridflag);
Smax1=min(Uend,Uend2);
Smax2=max(Uend,Uend2); 

%% call the directivity model

fDi=zeros(size(U));
for ii=1:size(U,2)
    [fD,fDi(:,ii),PhiRed,PhiRedi,PredicFuncs,Other]=Bea24(M,U(:,ii),T(:,ii),Smax1,Smax2,Ztor,Rake,Tdo,Version);
    
    S2(:,ii)=Other.S2;
    fs2(:,ii)=PredicFuncs.fs2;
    ftheta(:,ii)=PredicFuncs.ftheta;
    fdist(:,ii)=PredicFuncs.fdist;
    fGprime(:,ii)=PredicFuncs.fGprime;
    
end  

%% plot U T contours
figure;  set(gcf,'position',[311   188    747 391 ]); 
subplot(1,2,1)
    Z=[fliplr(0:-5:round(min(min(T)))) 5:5:round(max(max(T)))]; % contour interval
    V=[fliplr(0:-20:round(min(min(T)))) 20:20:round(max(max(T)))]; % label interval
    [c,h]=contour(SiteX,SiteY,T,Z); hold on
    clabel(c,h,V)
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis square
    title('GC2, T Coordinate')
    xlabel('Easting (km)')
    ylabel('Northing (km)')

subplot(1,2,2)
    Z=[fliplr(0:-5:round(min(min(U)))) 5:5:round(max(max(U)))]; % contour interval
    V=[fliplr(0:-20:round(min(min(U)))) 20:20:round(max(max(U)))]; % label interval
    [c,h]=contour(SiteX,SiteY,U,Z); hold on
    clabel(c,h,V)
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis square
    title('GC2, U Coordinate')
    xlabel('Easting (km)')
    ylabel('Northing (km)')

    
 %%  plot the directivity model
 
figure; set(gcf,'position',[454 166 996 758]);
subplot(2,3,1)
    contourf(SiteX,SiteY,S2,'linestyle','none'); hold on
    colorbar; 
    title('\itS2')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis equal
    ylabel('Northing (km)')
    axis([-80 80 -80 170])
subplot(2,3,2)
    contourf(SiteX,SiteY,fs2,'linestyle','none'); hold on
    colorbar; 
    title('\itf_{S2}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis equal
    axis([-80 80 -80 170])
subplot(2,3,3)
    contourf(SiteX,SiteY,ftheta,'linestyle','none'); hold on
    colorbar;
    title('\itf_{\theta}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    axis equal
    axis([-80 80 -80 170])
subplot(2,3,4)
    contourf(SiteX,SiteY,fGprime,'linestyle','none'); hold on
    colorbar; 
    title('\itf_G''')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    ylabel('Northing (km)')
    axis equal
    axis([-80 80 -80 170])
subplot(2,3,5)
    contourf(SiteX,SiteY,fdist,'linestyle','none'); hold on
    colorbar; 
    title('\itf_{dist}')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    axis equal
    axis([-80 80 -80 170])
subplot(2,3,6)
    contourf(SiteX,SiteY,exp(fDi),'linestyle','none'); hold on
    colorbar; 
    title('Amplification, T=3 sec')
    colormap(othercolor('BuDRd_18',256))
    for ii=1:nt
        plot(ftraces(ii).trace(:,1),ftraces(ii).trace(:,2),'k','linewidth',2)
    end
    plot(type.epi(1),type.epi(2),'kp','markerfacecolor','r','markersize',12)
    xlabel('Easting (km)')
    axis equal
    axis([-80 80 -80 170])
    clim([.6 1.8])

    
disp('Finished with example script')