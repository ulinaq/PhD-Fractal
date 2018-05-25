function [node,segment,x] = simpleEquation(node,prm,segment,vnode,vsegment,P)
Pin = prm.Pin; Qin = prm.Qin;
Co = node.n + vnode.n;
Cv = Co+segment.n;
eqCount = 1;
seg = 1:1:node.n+vnode.n;
j=1;
in = node.in;
out = vnode.out;
Nin = nnz(in);
Nout = nnz(out);
%% Create nonlinear system of equation from artery
F{eqCount} = @(x) Pin - x(1) - 8*prm.mu*segment.len(1)*Qin/(pi*segment.rad(1)^4);
F{eqCount+1} = @(x) Qin - x(Co+1);
eqCount = eqCount + 2;
for i = 2 : node.n
    segend = segment.end==i;
    lsg = seg(segend);                      % parent segment
    ch = segment.end(segment.init==i);      % 2 nodes in bifurcation
    sg = seg(segment.init==i);              % 2 segments in bifurcation
    cphi = cos(3*(segment.theta(segment.init==i)-segment.theta(segend))/4);

    if node.junction(i)==1              % equation in bifurcation
            % equation
        F{eqCount} = @(x) x(i-1) - x(ch(1)-1) - 8*prm.mu*segment.len(sg(1))*x(Co+sg(1))/(pi*segment.rad(sg(1))^4)...
            - prm.Rho*x(Co+lsg)^2/(2*pi^2*segment.rad(segend)^4)...
            *(1+x(Co+sg(1))^4*segment.rad(segend)^4/(x(Co+lsg)^2*segment.rad(sg(1))^4)...
            -2*segment.rad(segend)^2*cphi(1)*x(Co+sg(1))/(x(Co+lsg)*segment.rad(sg(1))^2));
        F{eqCount+1} = @(x) x(i-1) - x(ch(2)-1) - 8*prm.mu*segment.len(sg(2))*x(Co+sg(2))/(pi*segment.rad(sg(2))^4)...
            - prm.Rho*x(Co+lsg)^2/(2*pi^2*segment.rad(segend)^4)...
            *(1+x(Co+sg(2))^4*segment.rad(segend)^4/(x(Co+lsg)^2*segment.rad(sg(2))^4)...
            -2*segment.rad(segend)^2*cphi(2)*x(Co+sg(2))/(x(Co+lsg)*segment.rad(sg(2))^2));
        F{eqCount+2} = @(x) x(Co+lsg) - x(Co+sg(1)) - x(Co+sg(2));
        eqCount = eqCount + 3;
        else
            F{eqCount} = @(x) x(i-1) - x(Co) - P(j,1:Nin)*x(Co+seg(in)-1) + P(j,Nin+1:Nin+Nout)*x(Cv+seg(out)-1);
            eqCount = eqCount + 1;
            j=j+1;
        end
end
F{eqCount} = @(x) x(node.n+1) - x(node.n) - 8*prm.mu*vsegment.len(1)*Qin/(pi*vsegment.rad(1)^4);
F{eqCount+1} = @(x) Qin - x(Cv+1);
eqCount = eqCount + 2;
for i = 2 : vnode.n
    segend = vsegment.end==i;
    lsg = seg(segend);                      % parent segment
    ch = vsegment.end(vsegment.init==i);      % 2 nodes in bifurcation
    sg = seg(vsegment.init==i);              % 2 segments in bifurcation
    cphi = cos(3*(vsegment.theta(vsegment.init==i)-vsegment.theta(segend))/4);

    if vnode.junction(i)==1              % equation in bifurcation
            % equation
        F{eqCount} = @(x) x(node.n+i-1) - x(node.n+ch(1)-1) + 8*prm.mu*vsegment.len(sg(1))*x(Cv+sg(1))/(pi*vsegment.rad(sg(1))^4)...
            + prm.Rho*x(Cv+lsg)^2/(2*pi^2*vsegment.rad(segend)^4)...
            *(1+x(Cv+sg(1))^4*vsegment.rad(segend)^4/(x(Cv+lsg)^2*vsegment.rad(sg(1))^4)...
            -2*vsegment.rad(segend)^2*cphi(1)*x(Cv+sg(1))/(x(Cv+lsg)*vsegment.rad(sg(1))^2));
        F{eqCount+1} = @(x) x(node.n+i-1) - x(node.n+ch(2)-1) + 8*prm.mu*vsegment.len(sg(2))*x(Cv+sg(2))/(pi*vsegment.rad(sg(2))^4)...
            + prm.Rho*x(Cv+lsg)^2/(2*pi^2*vsegment.rad(segend)^4)...
            *(1+x(Cv+sg(2))^4*vsegment.rad(segend)^4/(x(Cv+lsg)^2*vsegment.rad(sg(2))^4)...
            -2*vsegment.rad(segend)^2*cphi(2)*x(Cv+sg(2))/(x(Cv+lsg)*vsegment.rad(sg(2))^2));
        F{eqCount+2} = @(x) x(Cv+lsg) - x(Cv+sg(1)) - x(Cv+sg(2));
        eqCount = eqCount + 3;
        else
            F{eqCount} = @(x) x(node.n+i-1) - x(Co) - P(j,1:Nin)*x(Co+seg(in)-1) + P(j,Nin+1:Nin+Nout)*x(Cv+seg(out)-1);
            eqCount = eqCount + 1;
            j=j+1;
        end
end

%% Solve
% Create Linear system for initial guess
guess = LinTreeSolver(prm,node,segment,vnode,vsegment,P);
% guess = ones(2*Co-2,1);
options = optimoptions('fsolve','FunctionTolerance',1e-8,'MaxFunctionEvaluations',200000,'MaxIterations',100000);
[x,fval,exitflag,output,jacobian] = fsolve(@(y) cellfun(@(x) x(y), F), guess, options);
end