% series solution to 
% EQN c_t = D c_{zz}
% BCs c(0,t) = 0, c(H,t) = k C_{pre} (t)
% IC c(z,0) = c_{init}
D = 1; % diffusion coefficient
T = 1; % final time
H = 1; % height of contact lens
N = 100; % terms in series
nv = 0:N;
c_init = 1; % initial concentration
h_post = 0; % change later
t_steps = 100; % dt = .01 ?
z_nodes = 100; % number of z nodes
zv = linspace(h_post,H+h_post,z_nodes);
tv = linspace(0,T,t_steps);
c = zeros(z_nodes,t_steps);
for j = 1:length(tv)
	c_pre = 1; % ????????????????
	for i = 1:length(zv)
		c(i,j) = ...
			% first term
			c_init* (4 ./ (2*nv + 1)/pi * sin( (2*nv+1)*pi*zv(i)/H ) ...
			*exp( -(2*nv+1).^2 * pi^2 / H^2 ) * D*tv(j) ...
			% second term
			k*c_pre(j) * ( zv(i)/H + 2*(-1).^nv ./ nv / pi ...
							   sin(nv*pi*zv(i)/H) ) ...
			k * trapz( ...%%% (-1)n TYPO IN PAPER????
					c_pre(1:j)* 2*(-1)^nv .* nv*pi*D/H^2 ...
					* sin(nv*pi*zv(j)/H) ...
					* exp(-nv.^2*pi^2/H^2*D*(t) ) %%% ??????)
					);
	end
end
