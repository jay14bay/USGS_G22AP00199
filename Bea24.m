%%  Matlab function for the Bayless et al. (2024) Directivity Model: Bea24
%   as described in USGS External Grants Report G22AP00199
%
%	Jeff Bayless (jeff.bayless@aecom.com) 
%	Created: Feb 2024
%	Copyright (c) 2024, Jeff Bayless, covered by GNU General Public License
%	All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% INPUT: M,U,T,Smax1,Smax2,Ztor,Rake,Period,Version
% M:            Moment magnitude, 6<=M<=8. 1x1 double.
% U,T:          the GC2 coordinates in km. Must both be nX1 doubles where n is the number of locations at which the model provides a prediction. Must be columnar.
% Smax1:        the maximum S in the antistrike direction for the scenario in km (defined to be a negative value) . 1x1 double.
% Smax2:        the maximum S in the strike direction for the scenario in km (a positive value). 1x1 double.
% Ztor:         the depth to top of rupture in km (positive value). 1x1 double.
% Rake:         the characteristic rupture rake angle, in degrees. 1x1 double. For strike slip ruptures the Rake is between -180 to -150 degrees, -30 to 30 degrees, or 150 to 180 
% Period:       the spectral period for which fD and PhiRed are requested, in sec. 0.01<=Period<=10. 1x1 double.
% Version:      a flag for determining the model version. 1->simulation-based. 2->NGA-W2 data-based. 1x1 double.

%% OUTPUT: fD,fDi,PhiRed,PhiRedi,PredicFuncs,Other
% fD:           the directivity adjustment in ln units. nx1000 double at 1000 log-spaced periods between 0.01 and 10 sec. The periods are provided in Other.Per
% fDi:          the directivity adjustment in ln units at user provided 'Period'. nx1 double.
% PhiRed:       the phi reduction. nx1000 double.
% PhiRedi:      the phi reduction at user provided 'Period'. nx1 double.
% PredicFuncs:  a struct with eight fields:
    % fG:       the period independent (uncentered) geometric directivity predictor. nx1 double.
    % fGprime:  the period independent centered geometric directivity predictor (includes distance and Ztor tapers). nx1 double.
    % fGbar:    the directivity predictor centering term. nx1 double.
    % fdist:    the distance taper. nx1 double.
    % fztor:    the Ztor taper. nx1 double.
    % ftheta:   the azimuthal component of the directivity predictor. nx1 double.
    % fs2:      the rupture travel distance component of the directivity predictor. nx1 double.
    % A:        the period- and mag-dependent lower and upper bound of fD. 1x1000 double.
% Other:        a struct with eight fields:
    % Per:      the periods at which fD and PhiRed are provided. 1x1000 double with 1000 log-spaced periods between 0.01 and 10 sec.
    % Rmax:     the maximum distance of the distance taper. 1x1 double.
    % Footprint:the index of sites within the footprint of the directivity effect (those with nonzero distance taper). nx1 logical.
    % Tpeak:    the peak period of the directivity effect. 1x1 double.
    % k:        the logistic function slope (model coefficient). 1x1 double.
    % Amax:     the limiting upper and lower bound of A (model coefficient). 1x1 double.
    % Rst:      the distance from the surface trace in km. nx1 double.
    % Ry0:      the Ry0 distance in km. nx1 double.
    % S2:       the generalized rupture travel distance parameter. nx1 double.
    % theta:    the angle theta in degrees. nx1 double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [fD,fDi,PhiRed,PhiRedi,PredicFuncs,Other]=Bea24(M,U,T,Smax1,Smax2,Ztor,Rake,Period,Version)

    % impose the limits on M; Table 4-3
    if M>8; error('Upper M limit is 8.0'); end
    if M<6; error('Lower M limit is 6.0'); end

%% (1) determine constants based on the model version

    if Version==1 % simulation-based
        Amax=0.54; k=1.58; SigG=0.38;
        % phi reduction model coeff
        PhiPer=[0.01	0.3	    0.4	     0.5	 0.75	 1	     1.5	 2	     3	     4	     5	     7.5	 10];
        e1=    [0.00    0.000   0.0003   0.011   0.038   0.072   0.107   0.143   0.172   0.189   0.195   0.206   0.200];
    elseif Version==2 % NGA-W2 event based
        Amax=0.34; k=1.58; SigG=0.26; 
        % phi reduction model coeff
        PhiPer=[0.01	0.3	    0.4	    0.5	    0.75	1	    1.5	    2	    3	    4	    5	    7.5	    10];
        e1=    [0.00	0.000	0.0024	0.0074	0.024	0.041	0.064	0.076	0.091	0.110	0.124	0.145	0.157];
    end
    Per=logspace(-2,1,1000); % define period array
    Tpeak=10.^(-2.15+0.404.*M); % [Eq. 4a]
    x=log10(Per./Tpeak); % narrow band gaussian centered on Tpeak
    A=Amax.*exp(-x.^2./(2*SigG^2)); % [Eq. 4]
    e1interp=interp1(log(PhiPer),e1,log(Per)); %interpolate e1 to periods in Per

%% (2) Calculate the period-independent predictors fG, fs2, ftheta, tapers

% (2a) Convert U to S
    S=zeros(size(U)); 
    uneg=U<0; upos=U>=0; % negative and positive indices of U, S is positive in the direction of strike and negative in the opposite
    S(uneg)=-min(abs(U(uneg)),abs(Smax1)); % converts U to S
    S(upos)=min(abs(U(upos)),abs(Smax2));  % converts U to S

% (2b) convert U to Ry0
    Ry=zeros(size(U));
    Ry(upos)=U(upos)-Smax2;
    Ry(uneg)=abs(U(uneg))-abs(Smax1);
    utween=U<=Smax2 & U>=Smax1; % in between, Ry0=0
    Ry(utween)=0;

% (2c) Calculate S2
    Srake=S.*cosd(Rake); % note this can change the sign of S depending on the rake angle, the resulting negative Srake values are opposite the direction of rake
    Dmin=3; % fixed value, km
    D=Dmin;
    S2=sqrt(D.^2+Srake.^2);  % positive everywhere. [Eq. 3b]

% (2d) predictor variable fs2
    fs2=log(S2);
    % apply the cap to Fs, at 465 km or approx. L for a M8, about 6.14 ln units
    fsCap=log(465);
    fs2(fs2>fsCap)=fsCap;

% (2e) angular predictor variables
    theta=abs(atan(T./U));  % [Eq. 3c] 
    theta(isnan(theta))=0;  % set to 0 when located exactly on the trace and atan is undefined
    ftheta=abs(cos(2.*theta));

% (2f) Distance taper
    R=sqrt(T.^2+Ry.^2+Ztor.^2);   % R is the distance from the surface trace  
    Rmax=-60+20.*M; % Rmax is 60km for M6, 80 km for M7 and higher, linear between [Eq. 3f]
    if M>7
        Rmax=80;    
    end
    Footprint=R<=Rmax;  % logical index of sites within the footprint of the directivity effect (nonzero distance taper)
    fdist=zeros(size(R));
    fdist(Footprint)=1-exp(-4*Rmax./R(Footprint) + 4); % [Eq. 3e]

% (2g) Ztor taper [Eq. 3d]
    if abs(Ztor)<20
        fztor=repmat(1-abs(Ztor)/20,size(R));
    else
        fztor=zeros(size(R));
    end

% (2h) fG, fGprime
    fG=fs2.*ftheta; % [Eq. 3a]
    fGbar=zeros(length(R),1);
    for ii=1:length(R)
        fGbar(ii,1)=centerfunc(Smax2,abs(Smax1),R(ii),D,Rake);
    end
    fGprime=(fG-fGbar).*fdist.*fztor; % [Eq. 3]
        
%% (3) Calculate fD

    fD=A.*( 2./ (1 + exp(-k.*fGprime) ) -1 ); % logistic function [Eq. 2] 
    [~,ti]=min(abs(Per-Period)); % get fD at the user requested period
    fDi=fD(:,ti);  
        
%% (4) Calculate PhiReduction 

    %  only apply reduction within the directivity model footprint
    PhiRed=repmat(e1interp,size(fD,1),1);
    PhiRed(not(Footprint),:)=0;

    % PhiRed at user requested period
    PhiRedi=PhiRed(:,ti);
    
%% (5) Format Output    

% period-independent predictor functions
    PredicFuncs.fG=fG;
    PredicFuncs.fGprime=fGprime;
    PredicFuncs.fGbar=fGbar;
    PredicFuncs.fdist=fdist;
    PredicFuncs.fztor=fztor;
    PredicFuncs.ftheta=ftheta;
    PredicFuncs.fs2=fs2;
    PredicFuncs.A=A;

% other 
    Other.Per=Per;
    Other.Rmax=Rmax;
    Other.Footprint=Footprint;
    Other.Tpeak=Tpeak;
    Other.k=k;
    Other.Amax=Amax;
    Other.Rst=R;
    Other.Ry0=Ry;
    Other.S2=S2;
    Other.theta=theta;

end

% this function calculates the centering term, fGbar, for a given rupture dimension, hypocenter location, rake angle, and distance
function [fGbar]=centerfunc(L1,L2,R,D,Rake)
    if R<0.1;  R=0.1; end

    dx=0.1;
    % I1 - between the ends of the fault, strike direction
    x=0:dx:L1; xr=x.*cosd(Rake);
    s=sqrt(xr.^2+D^2);
    y1=log(s).*abs(cos(2.*atan(R./x)));

    % I2 - between the ends of the fault, anti-strike direction
    x=0:dx:L2; xr=x.*cosd(Rake);
    s=sqrt(xr.^2+D^2);
    y2=log(s).*abs(cos(2.*atan(R./x)));

    % I3 - off the end of the fault, strike direction
    l=L1;  lr=l.*cosd(Rake);
    s=sqrt(lr^2+D^2);
    x=l+dx:dx:l+R; 
    r=sqrt(R^2-(x-l).^2);
    y3=log(s).*abs(cos(2.*atan(r./x)));

    % I4 off the end, anti-strike dir
    l=L2;   lr=l.*cosd(Rake);
    s=sqrt(lr^2+D^2);
    x=l+dx:dx:l+R;
    r=sqrt(R^2-(x-l).^2);
    y4=log(s).*abs(cos(2.*atan(r./x)));

    % total
    fGbar=mean([y1 y2 y3 y4]);

end