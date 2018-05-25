# PhD-Fractal
Starting Project in PhD journey
## Arterial structure
We use fractal geometry structure on this vascular model based on previous study (Arterial Tree Structure). The fractal parameter that we use in the structure are bifurcation exponent \xi=2.7, degree of symmetric (ratio between child vessel radius) $\gamma=0.7$, and a constant length to radius ratio $l/r=20$ \cite{Kassab1995,Takahashi2009,Qureshi2014}.

## Mathematical model
We simplify arterial vessel using fractal tree and output terminal nodes of vessel become point sources in the capillaries tissue. We assumed a certain blood flowrate is injected in the input terminal and the capillary tissue have similar static pressure in the whole domain. Blood flow in the output terminal would give additional (dynamic) pressure point as source in capillary domain. Using 0D model that represent arterial vessel as electrical circuit we can get system of equation for arterial vessel and solve it to obtain pressure and flowrate in output terminals.
