function [vnode,vsegment] = Tree_vein(prm)
%% Create fractal tree
% Create first vnode (initial vnode)
vnode.XData(1) = 0;
vnode.YData(1) = 160;
vnode.junction(1) = 0;             % indicate Junction in the vnode

% Create first vsegment (first trunk)
vsegment.init(1) = 1;                 % initial vnode
vsegment.end(1) = 2;                  % end vnode
vsegment.len(1) = prm.LtoR*prm.Rvo;    % length
vsegment.rad(1) = prm.Rvo;             % radius
vsegment.theta(1) = -pi/2;             % angles between each branch of fractal and vertical axes

% Create second vnode
vnode.XData(2) = vnode.XData(1)+vsegment.len*cos(vsegment.theta);
vnode.YData(2) = vnode.YData(1)+vsegment.len*sin(vsegment.theta);
vnodeCount = 2;

% Create remaining vsegments and vnodes
init = 2;
while init <= vnodeCount
    kpi=randperm(prm.N);
    if vsegment.rad(init-1)>prm.Rvmin
        vnode.junction(init) = 1;
        for i=1:1:prm.N
            [vnode,vsegment,vnodeCount] = drawvsegment(vnode,vsegment,prm,init,vnodeCount,i,kpi(i));
        end
    else vnode.junction(init) = 0;
    end
    init = init+1;
end
vnode.n = vnodeCount;
vsegment.n = length(vsegment.len);

clear gamma i init delta beta kpi; 

function [vnode,vsegment,vnodeCount] = drawvsegment(vnode,vsegment,prm,init,vnodeCount,i,num)
% vsegment count always less than vnode count by 1
vsegment.init(vnodeCount) = init;
vsegment.end(vnodeCount) = vnodeCount+1;
vsegment.len(vnodeCount) = vsegment.len(init-1)*prm.ratio(num);
vsegment.rad(vnodeCount) = vsegment.rad(init-1)*prm.ratio(num);
vsegment.theta(vnodeCount) = vsegment.theta(init-1)+(-1)^i*prm.vphi(num);
% create vnode in the end of new vsegment
vnodeCount = vnodeCount+1;
vnode.Parent(vnodeCount) = init;
vnode.XData(vnodeCount) = vnode.XData(init)+vsegment.len(vnodeCount-1)*cos(vsegment.theta(vnodeCount-1));
vnode.YData(vnodeCount) = vnode.YData(init)+vsegment.len(vnodeCount-1)*sin(vsegment.theta(vnodeCount-1));
end

end