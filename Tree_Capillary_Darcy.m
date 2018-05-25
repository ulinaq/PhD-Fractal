function [P,T,pvec,Grid,node,vnode] = Tree_Capillary_Darcy(node, segment, vnode, vsegment, prm)
%% Create darcy domain
in = node.junction==0;
in(1) = 0;
node.in = in;
out = vnode.junction==0;
out(1) = 0;
vnode.out = out;
DomainSize = [120 160]; % x and y
DomainLimit = [-60 60 0 160]; % xmin, xmax, ymin ymax
Reps = 3;   % Epsilon radii for dirac

% Grid for Darcy domain
Grid.Nx=120; Grid.hx=DomainSize(1)/Grid.Nx;
Grid.Ny=160; Grid.hy=DomainSize(2)/Grid.Ny;
Grid.Nz=1; Grid.hz=1;
Grid.K=ones(3,Grid.Nx,Grid.Ny);

X = ceil((node.XData(in)- DomainLimit(1))./Grid.hx);
Y = ceil((node.YData(in)- DomainLimit(3))./Grid.hy);
Xout = floor((vnode.XData(out)-DomainLimit(1))./Grid.hx);
Yout = floor((vnode.YData(out))./Grid.hy);

input = sub2ind([Grid.Nx Grid.Ny],X,Y);
output = sub2ind([Grid.Nx Grid.Ny],Xout,Yout);

%% Darcy Flow
% Grid sizes
Ny=Grid.Ny; Nx=Grid.Nx; N=Ny*Nx;
hx=Grid.hx; hy=Grid.hy; 

L = Grid.K.^(-1);
tx = 2*hy/hx;  ty = 2*hx/hy;

% Cell center
x = hx/2:hx:hx*Nx;
y = hy/2:hy:hy*Ny;
[ycell, xcell]=meshgrid(y,x);

% Transmissibility in x-direction: Nx+1 to give transmissibility at boundaries
Tymat = zeros(Nx,Ny+1); 
Txmat = zeros(Nx+1,Ny);

% Calculate transmissibilities at interfaces (general case- homo, hetero)
Tymat(:,2:Ny) = (ty)./(L(1,:,1:Ny-1)+L(1,:,2:Ny)); % Take harmonic mean
Txmat(2:Nx,:) = (tx)./(L(2,1:Nx-1,:)+L(2,2:Nx,:)); % Take harmonic mean -- No flow boundary condition: 0 for 1 and Nx+1 -> therefore 2:Nx

% Off-diagonal vectors of the system
% by reshaping transmissibility matrices
Tyvec1 = reshape(Tymat(:,1:Ny), Ny*Nx, 1);
Txvec1 = reshape(Txmat(1:Nx,:), Ny*Nx, 1);
Tyvec2 = reshape(Tymat(:,2:Ny+1), Ny*Nx, 1);
Txvec2 = reshape(Txmat(2:Nx+1,:), Ny*Nx, 1);

% Matrix of the system
Amat = spdiags([-Tyvec1,-Txvec1,Tyvec1+Txvec1+Tyvec2+Txvec2,-Txvec2,-Tyvec2],[Nx,1,0,-1,-Nx],Ny*Nx,Ny*Nx);
Amat(N+1,1) = 1;
a = [input,output];
z = zeros(Nx,Ny);
for i=1:length(a)
    % Create circle region for spreading epsilon
    z=Reps+1-(xcell-xcell(a(i))).^2-(ycell-ycell(a(i))).^2;
    area=(xcell-xcell(a(i))).^2+(ycell-ycell(a(i))).^2<=Reps;
    zn=z(area)*hx*hy/sum(z(area)*hx*hy);
    q=zeros(N+1,1);   % without source as default
    q(area)=zn;
    q(N+1) = 0;       % dirichlet BC set 0 at corner (0,0)
    % Solve system equation
    pvec(:,i) = Amat\q;
end

P=pvec([input,output],:);
T.y=Tymat;
T.x=Txmat;
end

