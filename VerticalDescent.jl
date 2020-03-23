"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct VERTICAL_DESCENT_SAMPLER <: SAMPLER
	z_init::Sampleable
	z_descent::Sampleable
	vx_init::Sampleable
	vz_init::Sampleable
	θ::Sampleable
	t_extra::Sampleable
	t_descent::Sampleable
end

function vertical_descent_sampler(;z_init = Uniform(300,600),
						z_descent = Normal(100,8),
						vx_init = Uniform(30,45),
						vz_init = Uniform(-500fpm2fps,-300fpm2fps),
						θ = Truncated(Normal(6,3),2,15),
						t_extra = Uniform(0.1,0.3), #Normal(10,3),
						t_descent = Normal(30,5))
	return VERTICAL_DESCENT_SAMPLER(z_init, z_descent, vx_init, vz_init, θ, t_extra, t_descent)
end

mutable struct VERTICAL_DESCENT <: UAM_TRAJECTORY
	z_init::Float64
	z_descent::Float64
	vx_init::Float64
	vz_init::Float64
	θ::Float64
	t_extra::Float64
	t_descent::Float64
	sampler::SAMPLER
	dt::Float64
	vz_max::Float64
	az_max::Float64
	ax_max::Float64
	p::Array{Float64,2}
	v::Array{Float64,2}
end

function vertical_descent(;z_init = 0.0,
						z_descent = 0.0,
						vx_init = 0.0,
						vz_init = 0.0,
						θ = 0.0,
						t_extra = 0.0,
						t_descent = 0.0,
						sampler = vertical_descent_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps,
						az_max = 0.3g,
						ax_max = 0.1g,
						p = zeros(1,3),
						v = zeros(1,3))
	return VERTICAL_DESCENT(z_init, z_descent, vx_init, vz_init, θ, t_extra, t_descent, sampler, dt,
							 vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::VERTICAL_DESCENT)
	# Glide portion #################################################
	# Get some needed info
	τ.vz_init = clamp(τ.vz_init, -τ.vz_max, Inf)
	ttot = abs(τ.z_init/τ.vz_init) + τ.t_extra*abs(τ.z_init/τ.vz_init)
	N = convert(Int64, floor(ttot/τ.dt))
	# Difference matrix
	D = zeros(N-1,N)
	for i = 1:(N-1)
	    D[i,i] = -1
	    D[i,i+1] = 1
	end
	D = D./τ.dt

	# Create and solve the optimization problem
	px = Variable(N)
	pz = Variable(N)
	vx = Variable(N)
	vz = Variable(N)

	objec = norm(D*vx) + λ*norm(vz*ft2m + vx*tand(τ.θ))

	constraints = [[px[i+1] == px[i] + vx[i]*τ.dt for i = 1:N-1];
				   [pz[i+1] == pz[i] + vz[i]*τ.dt for i = 1:N-1];
				   vx[1] == τ.vx_init; px[1] == 0;
				   vz[1] == τ.vz_init; pz[1] == τ.z_init;
				   vx[N] == 0; pz[N] == τ.z_descent; vz[N] == 0;
				   vx ≥ 0; pz ≥ 0; vz ≤ 0; vz ≥ -τ.vz_max;
				   D*vx ≤ τ.ax_max; D*vz*ft2m ≤ τ.az_max]

	prob1 = minimize(objec, constraints)
	solve!(prob1, ECOS.Optimizer(verbose=0))

	# Vertical descent portion ######################################
	Nd = convert(Int64, floor(τ.t_descent/τ.dt))
	# Difference matrix
	Dd = zeros(Nd-1,Nd)
	for i = 1:(Nd-1)
	    Dd[i,i] = -1
	    Dd[i,i+1] = 1
	end
	Dd = Dd./τ.dt

	pzd = Variable(Nd)
	vzd = Variable(Nd)

	objec = norm(Dd*vzd)

	constraints = [[pzd[i+1] == pzd[i] + vzd[i]*τ.dt for i = 1:Nd-1];
				   vzd[1] == 0; pzd[1] == τ.z_descent;
				   pzd[Nd] == 0; vzd[Nd] == 0;
				   pzd ≥ 0; vzd ≤ 0; vzd ≥ -τ.vz_max;
				   Dd*vzd*ft2m ≥ -τ.az_max; Dd*vzd*ft2m ≤ τ.az_max]

	prob2 = minimize(objec, constraints)
	solve!(prob2, ECOS.Optimizer(verbose=0))

	# Update the trajectory
	p = zeros(N+Nd,3)
	p[1:N,1] = px.value[:,1]
	p[N+1:end,1] = px.value[end,1]*ones(Nd)
	p[1:N,3] = pz.value[:,1]
	p[N+1:end,3] = pzd.value
	v = zeros(N+Nd,3)
	v[1:N,1] = vx.value[:,1]
	v[1:N,3] = vz.value[:,1]
	v[N+1:end,3] = vzd.value

	τ.p = copy(p)
	τ.v = copy(v)

	return Int(prob1.status) == 1
end