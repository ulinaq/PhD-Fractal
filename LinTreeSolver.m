function xlin = LinTreeSolver(prm,node,segment,vnode,vsegment,P)
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
G = zeros(2*Co-2);
b = zeros(2*Co-2,1);
G(eqCount,1) = 1;       b(eqCount) = Pin - 8*prm.mu*segment.len(1)*Qin/(pi*segment.rad(1)^4);
G(eqCount+1,Co+1) = 1;  b(eqCount+1) = Qin;
eqCount = eqCount + 2;
for i = 2 : node.n
    segend = segment.end==i;
    lsg = seg(segend);                      % parent segment
    ch = segment.end(segment.init==i);      % 2 nodes in bifurcation
    sg = seg(segment.init==i);              % 2 segments in bifurcation
    cphi = cos(3*(segment.theta(segment.init==i)-segment.theta(segend))/4);

    if node.junction(i)==1              % equation in bifurcation
        % equation
        G(eqCount,i-1) = 1;         G(eqCount,ch(1)-1) = -1;    G(eqCount,Co+sg(1)) = - 8*prm.mu*segment.len(sg(1))/(pi*segment.rad(sg(1))^4);
        G(eqCount+1,i-1) = 1;       G(eqCount+1,ch(2)-1) = -1;    G(eqCount+1,Co+sg(2)) = - 8*prm.mu*segment.len(sg(2))/(pi*segment.rad(sg(2))^4);
        G(eqCount+2,Co+lsg) = 1;    G(eqCount+2,Co+sg(1)) = -1; G(eqCount+2,Co+sg(2)) = -1; 
        eqCount = eqCount + 3;
    else
        G(eqCount,i-1) = 1;   G(eqCount,Co) = -1;  G(eqCount,Co+seg(in)-1) = - P(j,1:Nin);   G(eqCount,Cv+seg(out)-1) = P(j,Nin+1:Nin+Nout);
        eqCount = eqCount + 1;
        j=j+1;
    end
end
G(eqCount,node.n+1) = 1;    G(eqCount,node.n) = -1;   b(eqCount) = 8*prm.mu*vsegment.len(1)*Qin/(pi*vsegment.rad(1)^4);
G(eqCount+1,Cv+1) = 1;      b(eqCount+1) = Qin;
eqCount = eqCount + 2;
for i = 2 : vnode.n
    segend = vsegment.end==i;
    lsg = seg(segend);                      % parent segment
    ch = vsegment.end(vsegment.init==i);      % 2 nodes in bifurcation
    sg = seg(vsegment.init==i);              % 2 segments in bifurcation
    cphi = cos(3*(vsegment.theta(vsegment.init==i)-vsegment.theta(segend))/4);

    if vnode.junction(i)==1              % equation in bifurcation
        % equation
        G(eqCount,node.n+i-1) = 1;     G(eqCount,node.n+ch(1)-1) = -1;    G(eqCount,Cv+sg(1)) = 8*prm.mu*vsegment.len(sg(1))/(pi*vsegment.rad(sg(1))^4);
        G(eqCount+1,node.n+i-1) = 1;   G(eqCount+1,node.n+ch(2)-1) = -1;    G(eqCount+1,Cv+sg(2)) = 8*prm.mu*vsegment.len(sg(2))/(pi*vsegment.rad(sg(2))^4);
        G(eqCount+2,Cv+lsg) = 1;       G(eqCount+2,Cv+sg(1)) = -1;        G(eqCount+2,Cv+sg(2)) = -1; 
        eqCount = eqCount + 3;
    else
        G(eqCount,node.n+i-1) = 1; G(eqCount,Co) = -1;  G(eqCount,Co+seg(in)-1) = - P(j,1:Nin);   G(eqCount,Cv+seg(out)-1) = P(j,Nin+1:Nin+Nout);
        eqCount = eqCount + 1;
        j=j+1;  
    end
end
xlin = G\b;
end
