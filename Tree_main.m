clear all;
tic;

%% Parameters
prm = Tree_parameter(2.5847);
prm.LtoR = 10;        % Ratio between length and radius

%% Create Fractal tree
% Arterial tree
[node,segment] = Tree_artery(prm);
% venous tree
[vnode,vsegment] = Tree_vein(prm);
% Tree visualization
GV=graph(vsegment.init,vsegment.end);
GA=graph(segment.init,segment.end);

figure()
plot(GA,'XData',node.XData,'YData',node.YData);
title('Artery');
daspect([1 1 1])
hold on;
plot(GV,'XData',vnode.XData,'YData',vnode.YData);
title('Vein');
daspect([1 1 1])

%% Capillary Domain
[P,T,pvec,Grid,node,vnode] = Tree_Capillary_Darcy(node, segment, vnode, vsegment, prm);

%% Calculate pressure in nodes (Update April 26, 2018)
[node,segment,x] = simpleEquation(node,prm,segment,vnode,vsegment,P);
Co = node.n + vnode.n;
Cv = Co + segment.n;
seg = 1:1:Co;
q = [x([Co+seg(node.in)-1]);-1*x([Cv+seg(vnode.out)-1])];
Pmat = x(Co)+pvec*q;
Pmat = reshape(Pmat,Grid.Nx,Grid.Ny);
V = Flux(Pmat,Grid,T);

figure()
contourf(Pmat,20); colorbar; title('Pressure');
hold on;
quiver(V.qx,V.qy,'r'); %title('Fluxes');
toc

%% Calculate Flux
function V = Flux(Pmat,Grid,T)
   % Fluxes at staggered grid
   Ny=Grid.Ny; Nx=Grid.Nx; N=Ny*Nx;
       
   V.x = zeros(Nx,Ny+1);
   V.y = zeros(Nx+1,Ny);;
   V.x(:,2:Ny) = (Pmat(:,1:Ny-1)-Pmat(:,2:Ny)).*T.y(:,2:Ny);
   V.y(2:Nx,:) = (Pmat(1:Nx-1,:)-Pmat(2:Nx,:)).*T.x(2:Nx,:);

   V.qx = movmean(V.x,2,2,'Endpoints','discard');
   V.qy = movmean(V.y,2,1,'Endpoints','discard');
   V.u = V.qx.^2+V.qy.^2;
   u = reshape(V.u,N,1);
   V.range = max(u)-min(u);
end