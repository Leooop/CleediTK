function [vx, vy] = wedgeAnalyticalSolution2(X, Y, analyticalPARAMS, EQTNS, DYN, GEOM)

%% define values and the coordinate system

delta = deg2rad(GEOM.slab_dip); % slab dip in radians

a = GEOM.lplate_thickness; % downgoing plate thickness in m
Zw = GEOM.uplate_thickness; % overriding plate thickness in m
V = DYN.v_slab; % plate velocity in m/yr

%maxX = analyticalPARAMS.xSize;

% cartesian coordinates
%x = -600:600; % horizontal dimension in km
%x = x*1e3; % horizontal dimention in m
%y = 0:800; % vertical dimension (pos down) in km
%y = y*1e3; % vertical dimension in m

Xorig = X;
Yorig = Y;

X = X-Zw/tan(delta); % put wedge corner in left third of box
Y = Yorig - Zw; % set Y indexing to start at base of upper plate

%[XX, YY] = meshgrid(x, y);

%YY = a-YY;

% wedge polar coordinates
%[theta, ~] = cart2pol(X,Y);
theta = atan2(Y,X);

%theta(theta>delta) = nan;

% define slab top cartesian coordinates

%YYp = X*cos(delta) + Y*sin(delta);
%XXp = X*sin(delta) - Y*cos(delta);

% under subducting plate polar coordinates
Y2 = Y-a+Zw; % (Y=0 at base of downgoing plate)
X2 = -X-a/sin(delta)-(Zw-a)/tan(delta);
%[theta2, ~] = cart2pol(X2, Y2);
theta2 = atan2(Y2,X2);

%% calculate velocity

% in slab/plate, boundaries

omega = pi-delta;

%inOverRidingPlate = Y<=0 & theta<delta;
%Ux(inOverRidingPlate) = 0;
%Uy(inOverRidingPlate) = 0;

%inHorizontalSlab = Y<=0 & abs(theta)>omega;
%Ux(inHorizontalSlab) = V;
%Uy(inHorizontalSlab) = 0;
 
%inDippingSlab = Y>0 & abs(theta)>-omega;
%Ux(inDippingSlab) = cos(delta)*V;
%Uy(inDippingSlab) = sin(delta)*V;

% beneath subducting plate
%theta2(theta2<0) = nan;
%theta2(theta2>(pi-delta)) = nan;

% define some additional useful variables
%theta2_2pi = theta2;
%theta2_2pi(theta2_2pi<=0) = theta2_2pi(theta2_2pi<=0) + 2*pi;

%theta2pi = theta;
%theta2pi(theta2pi<=0) = theta2pi(theta2pi<=0) + 2*pi;

%% creat masks for different zones and assign vx, vy piecewise

vx = nan(size(X));
vy = nan(size(X));

% in the overriding plate
overridingPlate = theta>-pi+delta & theta<pi & Y<=0;
vx(overridingPlate) = 0;
vy(overridingPlate) = 0;

% in the dipping slab
dippingSlab = Xorig<=Yorig/tan(delta);% | omega-theta2<pi;
vx(dippingSlab) = cos(delta)*V;
vy(dippingSlab) = sin(delta)*V;

% in the arc corner
arcCorner = theta<=delta & theta>=0;
arcTheta = theta(arcCorner);
%Utheta = V/(delta^2 - sin(delta)^2)*((delta-arcTheta).*sin(arcTheta)*sin(delta) - delta*arcTheta.*sin(delta-arcTheta));
%Ur = -V/(delta^2 - sin(delta)^2)*((arcTheta-delta).*cos(arcTheta)*sin(delta)+sin(arcTheta)*sin(delta)-delta*sin(arcTheta-delta)-delta*arcTheta.*cos(arcTheta-delta));

if analyticalPARAMS.n>0
    Utheta = EQTNS.arcVth_fxn(arcTheta);
    Ur = EQTNS.arcVr_fxn(arcTheta);
       
    vx(arcCorner) = (Ur.*cos(arcTheta) - Utheta.*sin(arcTheta));
    vy(arcCorner) = (Utheta.*cos(arcTheta) + Ur.*sin(arcTheta));
else
    vx(arcCorner) = 0;
    vy(arcCorner) = 0;
end

origin = Xorig == Zw & Y == 0;
vx(origin) = cos(delta)*V;
vy(origin) = sin(delta)*V;

%in the oceanic side corner

slabCorner = theta2>=0 & theta2<=omega;
slabTheta2 = theta2(slabCorner);
%Utheta2 = V/(omega+sin(omega))*((slabTheta2-omega).*sin(slabTheta2)+slabTheta2.*sin(slabTheta2-omega));  % positive here (i.e. extra negative) because theta2 goes in opposite direction(?)
%Ur2 = V/(omega+sin(omega))*((slabTheta2-omega).*cos(slabTheta2)+sin(slabTheta2)+sin(slabTheta2-omega)+slabTheta2.*cos(slabTheta2-omega));

if analyticalPARAMS.n>0
    Utheta2 = EQTNS.oceanVth_fxn(slabTheta2);
    Ur2 = EQTNS.oceanVr_fxn(slabTheta2);

    vx(slabCorner) = -(Ur2.*cos(slabTheta2) - Utheta2.*sin(slabTheta2));
    vy(slabCorner) = (Utheta2.*cos(slabTheta2) + Ur2.*sin(slabTheta2));
else
    vx(slabCorner) = 0;
    vy(slabCorner) = 0;
end

% at the origin (HOW SHOULD THIS BE DEFINED?)
% origin = X == 0 & Y == 0;
% vx(origin) = cos(delta)*V;
% vy(origin) = sin(delta)*V;


% % in wedge
% Utheta = V/(delta^2 - sin(delta)^2)*((delta-theta).*sin(theta)*sin(delta) - delta*theta.*sin(delta-theta));
% Ur = -V/(delta^2 - sin(delta)^2)*((theta-delta).*cos(theta)*sin(delta)+sin(theta)*sin(delta)-delta*sin(theta-delta)-delta*theta.*cos(theta-delta));
% 
% %Ux = (Ur-tan(theta).*Utheta)./(cos(theta).*(1+tan(theta).^2));
% %Uy = -(Utheta-tan(theta).*Ur)./(cos(theta).*(1+tan(theta).^2));
% 
% Ux = (Ur.*cos(theta) - Utheta.*sin(theta));
% Uy = -(Utheta.*cos(theta) - Ur.*sin(theta));
% 
% Utheta2 = V/(omega+sin(omega))*((theta2-omega).*sin(theta2)+theta2.*sin(theta2-omega));  % positive here (i.e. extra negative) because theta2 goes in opposite direction(?)
% Ur2 = V/(omega+sin(omega))*((theta2-omega).*cos(theta2)+sin(theta2)+sin(theta2-omega)+theta2.*cos(theta2-omega));
% 
% Ux2 = -(Ur2.*cos(2*pi-theta2) - Utheta2.*sin(2*pi-theta2));
% Uy2 = -(Utheta2.*cos(2*pi-theta2) + Ur2.*sin(2*pi-theta2));

%Ux2 = -(Ur2-tan(theta2).*Utheta2)./(cos(theta2).*(1+tan(theta2).^2));
%Uy2 = -(Utheta2-tan(theta2).*Ur2)./(cos(theta2).*(1+tan(theta2).^2));


% %% assign vx & vy from the different components
% 
% vx = nan(size(X));
% vy = nan(size(X));
% 
% % define some additional useful variables
% theta2_2pi = theta2;
% theta2_2pi(theta2_2pi<=0) = theta2_2pi(theta2_2pi<=0) + 2*pi;
% 
% theta2pi = theta;
% theta2pi(theta2pi<=0) = theta2pi(theta2pi<=0) + 2*pi;
% 
% % in the arc corner
% arcCorner = theta<delta & theta>0;
% vx(arcCorner) = Ux(arcCorner);
% vy(arcCorner) = Uy(arcCorner);
% 
% % in the oceanic side corner
% slabCorner = theta2>0 & theta2<omega;
% vx(slabCorner) = Ux2(slabCorner);
% vy(slabCorner) = Uy2(slabCorner);
% 
% % in the overriding plate
% overridingPlate = theta2pi>pi+delta & theta<=0;
% vx(overridingPlate) = 0;
% vy(overridingPlate) = 0;
% 
% % in the horizontal slab
% horizontalSlab = theta2_2pi>(2*pi+omega)/2;
% vx(horizontalSlab) = V;
% vy(horizontalSlab) = 0;
% 
% % in the dipping slab
% dippingSlab = theta2_2pi>omega & theta2_2pi<=(2*pi+omega)/2 & theta2pi>=delta & theta2pi<=pi+delta;% | omega-theta2<pi;
% vx(dippingSlab) = cos(delta)*V;
% vy(dippingSlab) = sin(delta)*V;

if any(isnan(vx(:))) || any(isnan(vy(:)))
    error('there''s a hole in the velocity field!')
end

end