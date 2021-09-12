function [vx,vy]=cornerflow(X,Y,DYN,GEOM)


%% parameters
%analyticalPARAMS.Ts = 273;
%analyticalPARAMS.Tp = 1573;
%analyticalPARAMS.slabAge = 1;
%analyticalPARAMS.slabAge_s = analyticalPARAMS.slabAge*1e6*365.25*24*3600; % slab age in seconds
%analyticalPARAMS.kappa = PARAMS.kref/PARAMS.rhoref/PARAMS.cpref;
%analyticalPARAMS.H_slab = 0; %erfinv(0.999)*2*sqrt(analyticalPARAMS.kappa*analyticalPARAMS.slabAge_s); % make sure slab is ~thick enough for plate model temperature

% analyticalPARAMS.xSize = 1e5;
% xx = linspace(-analyticalPARAMS.xSize,analyticalPARAMS.xSize,100);
% yy = linspace(0,analyticalPARAMS.xSize,100);
% [X,Y]=meshgrid(xx,yy);


%analyticalPARAMS.slabVel = analyticalPARAMS.slabVel/100/365.25/24/3600; % slab velocity in m/sec

analyticalPARAMS.n = 1; % 1 for newtownian analytical solution (McKenzie), 3 for power law (Tobisch)

V = DYN.v_slab;
alpha = deg2rad(GEOM.slab_dip);
omega = pi-alpha;


%% analytical corner flow 
syms TH_0 TH_1 PI_0 th_0 th
    
    % upper plate side

    streamFxn = PI_0/2*th.*sin(th+th_0) + TH_0*sin(th) + TH_1*cos(th);
    streamFxn1 = diff(streamFxn, th);
    streamFxnN = -PI_0/2*th.*sin(th+th_0) - TH_0*sin(th) - TH_1*cos(th);
    
    bc1 = subs(streamFxn1, th, 0) == 0;
    bc2 = subs(streamFxnN, th, 0) == 0;
    bc3 = subs(streamFxn1, th, alpha) == V;
    bc4 = subs(streamFxnN, th, alpha) == 0;

    S1 = solve([bc1 bc2 bc3 bc4], [TH_0, TH_1, PI_0, th_0], 'Real', true);
    
    TH_0arc = eval(S1.TH_0(2));
    TH_1arc = eval(S1.TH_1(2));
    th_0arc = eval(S1.th_0(2));
    PI_0arc = eval(S1.PI_0(2)); 
    
    EQTNS.arcVth_fxn = @(th) -(PI_0arc/2*th.*sin(th + th_0arc) + TH_0arc*sin(th) + TH_1arc*cos(th));
    EQTNS.arcVr_fxn  = @(th) PI_0arc/2*(sin(th+th_0arc) + th.*cos(th + th_0arc)) + TH_0arc*cos(th) - TH_1arc*sin(th);
    
    % repeat for oceanic wedge
    
    bc1 = subs(streamFxn1, th, 0) == -V;
    bc2 = subs(streamFxnN, th, 0) == 0;
    bc3 = subs(streamFxn1, th, omega) == V;
    bc4 = subs(streamFxnN, th, omega) == 0;
    
    
    S2 = solve([bc1 bc2 bc3 bc4], [TH_0, TH_1, PI_0, th_0], 'Real', true, 'PrincipalValue', true);
    
    TH_0oc = eval(S2.TH_0);
    TH_1oc = eval(S2.TH_1);
    th_0oc = eval(S2.th_0);
    PI_0oc = eval(S2.PI_0);    
    
    EQTNS.oceanVth_fxn = @(th) -(PI_0oc/2*th.*sin(th + th_0oc) + TH_0oc*sin(th) + TH_1oc*cos(th));
    EQTNS.oceanVr_fxn  = @(th) PI_0oc/2*(sin(th+th_0oc) + th.*cos(th + th_0oc)) + TH_0oc*cos(th) - TH_1oc*sin(th);

%% flow field


[vx, vy] = wedgeAnalyticalSolution2(X, Y, analyticalPARAMS, EQTNS, DYN, GEOM);