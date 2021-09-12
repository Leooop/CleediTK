function [P,T,X_mat]=thermokin(OPT,DYN,GEOM,KIN,THERMAL,BC)

%%%%%% model design %%%%%%%%%%
x      = linspace(GEOM.x0(1),GEOM.x0(2),GEOM.res(1)); 
y      = linspace(GEOM.y0(1),GEOM.y0(2),GEOM.res(2));
[X,Y]  = meshgrid(x,y) ;

[VX,VY]=cornerflow(X,Y,DYN,GEOM);

Y  = flip(-Y);
X  = flip(X);
vx = flip(VX);
vy =flip(-VY);



%physics
T2D    = init_geotherm(X,Y,vx,vy,THERMAL,GEOM,BC);
[verts,averts] = streamslice(X,Y,vx,vy);

% user input
f1 = figure(1);
hp = pcolor(X,Y,T2D); set(hp,'edgecolor','none'); axis equal tight; colorbar;
sl = streamline(verts);  
[STARTX,STARTY]=ginput(1);
close(f1)

% compute streamline and extract P-T
XY = stream2(X,Y,vx,vy,STARTX,STARTY);
iXY = interpstreamspeed(X,Y,vx,vy,XY,DYN.dt);
F_temp = griddedInterpolant(X',Y',T2D');
T      = F_temp(iXY{1});
P = -iXY{1}(:,2)*3300*9.81;

% kinetic
[X_mat, id_vec, t_vec] = kinetic(P,T,KIN,DYN);

% plot kinetic on 
f2 = figure(2);
c = X_mat(:,OPT.color_phase);
sz = 50;
if OPT.hold.fig2 == 0
    close(f2)
    f2 = figure(2);
    hp = contour(X,Y,T2D,'r','ShowText','on'); axis equal tight; 
    sl = streamline(verts); hold on;
end
h=scatter(iXY{1}(:,1),iXY{1}(:,2),sz,c,'.'); colorbar;
caxis([0 max(c)])
%set(h,'Color',X_mat(:,9),'Marker','.')

% vizualisation
f3 = figure(3);
if OPT.hold.fig3 == 0
    close(f3)
    f3 = figure(3);
end
contour(KIN.Tdiag,KIN.Pdiag,KIN.phases_diag,'k')
hold on
[nt,nc] = size(c);
scatter(T(1:nt),P(1:nt),sz,c,'.')
caxis([0 max(c)])
axis ij
%plot(T_vec,1e6*T_vec)
%plot(T_vec,1e6*T_vec + 0.3e9)

xlabel('T (°C)')
ylabel('P (Pa)')
colorbar
grid on

%
phases = unique(id_vec);
cphases = lines(KIN.nphases);
f4 = figure(4);
if OPT.hold.fig4 == 0
    close(f4)
    f4 = figure(4);
    %set(gca, 'ColorOrder', hsv(length(phases)), 'NextPlot', 'replacechildren')
end
hold on
for i = 1:length(phases)
    id = phases(i);
    plot(t_vec,X_mat(:,id),'LineWidth',3,'DisplayName',strcat('phase',' ',num2str(id)),'Color',cphases(id,:))
    %scatter(t_vec,X_mat(:,id),4,'.')
end
legend()

% write PT path on file if needed
if OPT.write_PT_path == 1
    fid = fopen('PTpath.dat','w');
    fprintf(fid,'%e\t%e\n',[P, T]') % output in K, bar
    fclose(fid)
end

% figure(4);plot(T,-P) ;hold on;
% for i=1:size(verts,2)
% T      = F_temp(verts{i});
% P = verts{i}(:,2)*3300*10;
% PT{i}=[T,-P];
% figure(4);plot(T,-P) ;hold on;
% end


function T2D=init_geotherm(X,Y,VX,VY,THERMAL,GEOM,BC)

bctype.top   = 1; % T = const
bctype.bot   = 3; % T imposé en inflow et no flux conductif en outlow ( gradient orthogonal au bord nul)
bctype.left  = 3; % T selon half-space cooling model avec Tmax = bcval.bot
bctype.right = 3; % T imposé en inflow avec un gradient linéaire entre 0 et bcval.right et no flux conductif en outlow ( gradient orthogonal au bord nul)
bcval.top    = 0;
bcval.left   = 0;
bcval.right  = 900; % T à la limite inflow / outflow
bcval.bot    = 1300; 
bcval.ageMa   =80;

% diffusivity field
chi2D  = THERMAL.diffusivity*ones(size(X));

% heat production field
is_below_slab = (X <= 0 | Y > tand(GEOM.slab_dip) .* X);
is_upperplate = (~is_below_slab & Y <= GEOM.radiogenic_uplate_thickness);
is_lowerplate = (is_below_slab & (Y <= GEOM.radiogenic_lplate_thickness | Y < (tand(GEOM.slab_dip).*X + GEOM.radiogenic_lplate_thickness)));
is_radiogenic = (is_upperplate | is_lowerplate);
P2D           = THERMAL.heat_prod .* is_radiogenic; %zeros(size(X)); % heat production term degree/sec -> H / rho*Cp ~ 1e-11, T/s
 

[T2D]= compute_stationary_geotherm(X,Y,VX,VY,chi2D,P2D,BC,THERMAL);


% U=1e-9; 
% [VXa,VYa,VXb,VYb]=cornerflow(U,X,Y);
% 
% figure(4);quiver(X,Y,VXb,VYb);
% figure(5);quiver(X,Y,VXa,VYa); 
% 
% function [u_a,v_a,u_s,v_s]=cornerflow(U,X,Y)
% 
% A= 0;
% B= 0;
% C= -pi*U*sqrt(2)/2/(2-pi^2/4);
% D= -U*sqrt(2)*(2-pi/2)/(2-pi^2/4);
% 
% u_a = -D*atan(Y./X)+(C*X+D*Y).*(-X./(X.^2+Y.^2));
% v_a = -C*atan(Y./X)+(C*X+D*Y).*(-Y./(X.^2+Y.^2));
% 
% C = U/(9*pi^2/4-2)*(2-sqrt(2)/(1+3*pi/2)*(3*pi/2+9*pi^2/4));
% D = U/(9*pi^2/4-2)*(sqrt(2)*(2+3*pi/2)-2*(1+3*pi/2));
% 
% u_s = -D*atan(Y./X)+(C*X+D*Y).*(-X./(X.^2+Y.^2));
% v_s = -C*atan(Y./X)+(C*X+D*Y).*(-Y./(X.^2+Y.^2));





function [T2D]= compute_stationary_geotherm(X,Y,VX,VY,chi2D,P2D,BC,THERMAL)
%numerical parameters
nn         = size(X); %node number 
dX         = [Y(2,1)-Y(1,1),X(1,2)-X(1,1)] ; %node spacing 
nnrm       = reshape(1:prod(nn),nn(1),nn(2)); %  nodes numbering
ncoef      = nnrm(1:nn(1),1:nn(2)); % reordering for vectorisation
node_left  = nnrm(:,1); 
node_right = nnrm(:,end);
node_bot   = nnrm(1,:);
node_top   = nnrm(end,:);
%allocation of matrix and RHS
    A    = spalloc(prod(nn),prod(nn),5*prod(nn));
    R    = zeros(prod(nn),1);
% compute derivative of diffusivity
[chidx, chidy]=gradient(chi2D,dX(2),dX(1));
% vector of coeffients
kdx = chidx(ncoef(:));
kdy = chidy(ncoef(:));
vy1  = VY(ncoef(:));
vx1  = VX(ncoef(:));
e = chi2D(ncoef(:));
% matrix assembly
A2 =spdiags(-[e/dX(1)^2+kdy/2/dX(1)-vy1/2/dX(1), ...
              e/dX(1)^2-kdy/2/dX(1)+vy1/2/dX(1), ...
              -2*e/dX(1)^2-2*e/dX(2)^2, ...
              e/dX(2)^2+kdx/2/dX(2)-vx1/2/dX(2), ...
              e/dX(2)^2-kdx/2/dX(2)+vx1/2/dX(2)],...
              [-1 1 0 -nn(1) nn(1)], prod(nn), prod(nn));
          A2=A2';
R  = P2D(:); 
% boudary conditions
A2(node_right,:)   = 0; 
A2(node_right,node_right)   = eye(length(node_right));
R(node_right)                     = 0; 
switch BC.right.type
    case 1 % const temperature
    R(node_right)=BC.right.val;
    case 2 %
    A2(node_right,node_right-nn(1)) = A2(node_right,node_right-nn(1))-eye(length(node_right)); 
    case 3
    inflownodes      =  node_right(VX(node_right)<=0);   
    outflownodes     =  node_right(VX(node_right)>0);    
    R(inflownodes)   = Y(inflownodes)*BC.right.val/min(Y(inflownodes)); 
    A2(outflownodes,outflownodes-nn(1)) = A2(outflownodes,outflownodes-nn(1))-eye(length(outflownodes));
end 

      
A2(node_top,:)         = 0;        
A2(node_top,node_top)  =  eye(length(node_top)); 
R(node_top)            = BC.top.val; 
if BC.top.type == 2
   A2(node_top,node_top-1) = A2(node_top,node_top-1) - eye(length(node_top));
end


A2(node_bot,:)           = 0;
A2(node_bot,node_bot)    = eye(length(node_bot));
R(node_bot)              = 0;
 switch BC.bot.type
     case 2
  A2(node_bot,node_bot+1) = A2(node_bot,node_bot+1) - eye(length(node_bot));
     case 3
    inflownodes      =  node_bot(VY(node_bot)>=0);   
    outflownodes     =  node_bot(VY(node_bot)<0);    
    R(inflownodes)   =  BC.bot.val; 
    A2(outflownodes,outflownodes+1) = A2(outflownodes,outflownodes+1)-eye(length(outflownodes));
end

A2(node_left,:)          = 0; 
A2(node_left,node_left)  = eye(length(node_left));
R(node_left)             = 0;

switch BC.left.type
    case 2
      A2(node_left,node_left+nn(1)) = A2(node_left,node_left+nn(1)) - eye(length(node_left));
    case 3
      R(node_left)  =   BC.top.val +(BC.bot.val-BC.top.val)*erf(-Y(node_left)/sqrt(2*THERMAL.diffusivity*BC.age_slab));    
end

T = A2\R;
T2D= reshape(T,size(nnrm)); 

     
