function [VX,VY,X,Y]=main_term()
x0    = [-500e3 1000e3]; %width of box
y0    = [-250e3 0]; %length of box
det   = [100,25]; %resolution of the map in x,y
v_slab = 1; %cm/yr
%%%%%% model design %%%%%%%%%%
x      = linspace(x0(1),x0(2),det(1)); 
y      = linspace(y0(1),y0(2),det(2));
[X,Y]  = meshgrid(x,y) ;
chi2D  = 1e-6*ones(size(X));
P2D    = zeros(size(X));
%VX     = 3e-10*ones(size(X));
%VY     = 0*1.5e-9*ones(size(X)); 
[VX,VY]=cornerflow(X,-Y,v_slab);
VY=-VY;

save('corner_flow3','X','Y','VX','VY')

conti=0
if conti==1

VY = -VY;
VX=VX*0;
T_haut = 0;
T_bas   = 1300; 
[T2D]= compute_stationary_geotherm(X,Y,VX,VY,chi2D,P2D,T_haut,T_bas);

figure()
pcolor(X,Y,T2D);
U=1e-9; 
%[VXa,VYa,VXb,VYb]=cornerflow(U,X,Y); % HERE


% figure(4);quiver(X,Y,VXb,VYb);
figure(5);quiver(X,Y,VX,VY); 
% 
% function [u_a,v_a,u_s,v_s]=cornerflow(U,X,Y)
% 
% A= 0
% B= 0
% C= -pi*U*sqrt(2)/2/(2-pi^2/4)
% D= -U*sqrt(2)*(2-pi/2)/(2-pi^2/4)
% 
% u_a = -D*atan(Y./X)+(C*X+D*Y).*(-X./(X.^2+Y.^2));
% v_a = -C*atan(Y./X)+(C*X+D*Y).*(-Y./(X.^2+Y.^2));
% 
% C = U/(9*pi^2/4-2)*(2-sqrt(2)/(1+3*pi/2)*(3*pi/2+9*pi^2/4))
% D = U/(9*pi^2/4-2)*(sqrt(2)*(2+3*pi/2)-2*(1+3*pi/2));
% 
% u_s = -D*atan(Y./X)+(C*X+D*Y).*(-X./(X.^2+Y.^2));
% v_s = -C*atan(Y./X)+(C*X+D*Y).*(-Y./(X.^2+Y.^2));

end



function [T2D]= compute_stationary_geotherm(X,Y,VX,VY,chi2D,P2D,T_haut,T_bas)
%param?tres num?riques
nn    = size(X); %nombre de noeud 
dX    = [Y(2,1)-Y(1,1),X(1,2)-X(1,1)] ; 
nnrm          = reshape(1:prod(nn),nn(1),nn(2));
noeuds_gauche   = nnrm(:,1); 
noeuds_droite    = nnrm(:,end);
noeuds_bas = nnrm(1,:);
noeuds_haut = nnrm(end,:);
%boucle sur le temps

    A    = spalloc(prod(nn),prod(nn),5*prod(nn));
    R    = zeros(prod(nn),1);

[chidx, chidy]=gradient(chi2D,dX(2),dX(1)); 
ncoef = nnrm(1:nn(1),1:nn(2));
kdx = chidx(ncoef(:));
kdy = chidy(ncoef(:));
vy1  = VY(ncoef(:));
vx1  = VX(ncoef(:));

e = chi2D(ncoef(:));
A2 =spdiags(-[e/dX(1)^2+kdy/2/dX(1)-vy1/2/dX(1)^2, ...
              e/dX(1)^2-kdy/2/dX(1)+vy1/2/dX(1)^2, ...
              -2*e/dX(1)^2-2*e/dX(2)^2, ...
              e/dX(2)^2+kdx/2/dX(2)-vx1/2/dX(2), ...
              e/dX(2)^2-kdx/2/dX(2)+vx1/2/dX(2)],...
              [-1 1 0 -nn(1) nn(1)], prod(nn), prod(nn));
          A2=A2';
%A2 = A2 + spdiags(0.5*[-kdx/dX(1) kdx/dX(1) 0*kdx  kdy/dX(2)  kdy/dX(2)], [-1 1 0 -nn(1) nn(1)], prod(nn), prod(nn));
R  = P2D(:); 
% conditions aux limites
A2(noeuds_droite,:)   = 0; 
A2(noeuds_droite,noeuds_droite)   = eye(length(noeuds_droite)); 
A2(noeuds_droite,noeuds_droite-nn(1)) = A2(noeuds_droite,noeuds_droite-nn(1))-eye(length(noeuds_droite)); 

       % temp?rature constante en haut et en bas
A2(noeuds_haut,:)   = 0;        
A2(noeuds_haut,noeuds_haut)           =  eye(length(noeuds_haut)); 

A2(noeuds_bas,:)   = 0; 
A2(noeuds_bas,noeuds_bas)             = eye(length(noeuds_bas));
A2(noeuds_bas,noeuds_bas)             = A2(noeuds_bas,noeuds_bas+nn(1)) - eye(length(noeuds_bas));



A2(noeuds_gauche,:)   = 0; 
A2(noeuds_gauche,noeuds_gauche)       = eye(length(noeuds_gauche));
A2(noeuds_gauche,noeuds_gauche+nn(1)) = A2(noeuds_gauche,noeuds_gauche+nn(1)) - eye(length(noeuds_gauche));



R(noeuds_droite)                      = 0; 
R(noeuds_haut)                        = T_haut; 
R(noeuds_bas)                         = 0;% T_bas;
R(noeuds_gauche)                      =  0;
T = A2\R;
T2D= reshape(T,size(nnrm)); 

     
