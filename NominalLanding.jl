"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct NOMINAL_LANDING_SAMPLER <: SAMPLER
	z_init::Sampleable
	vx_init::Sampleable
	vz_init::Sampleable
	θ::Sampleable
	t_extra::Sampleable
end

function nominal_landing_sampler(;z_init = Uniform(300,600),
						vx_init = Uniform(30,45),
						vz_init = Uniform(-500fpm2fps,-300fpm2fps),
						θ = Truncated(Normal(6,3),2,15),
						t_extra = Uniform(0.1,0.3)) #Normal(10,3))
	return NOMINAL_LANDING_SAMPLER(z_init, vx_init, vz_init, θ, t_extra)
end

mutable struct NOMINAL_LANDING <: UAM_TRAJECTORY
	z_init::Float64
	vx_init::Float64
	vz_init::Float64
	θ::Float64
	t_extra::Float64
	sampler::SAMPLER
	dt::Float64
	vz_max::Float64
	az_max::Float64
	ax_max::Float64
	p::Array{Float64,2}
	v::Array{Float64,2}
end

function nominal_landing(;z_init = 0.0,
						vx_init = 0.0,
						vz_init = 0.0,
						θ = 0.0,
						t_extra = 0.0,
						sampler = nominal_landing_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps,
						az_max = 0.3g_ft,
						ax_max = 0.1g,
						p = zeros(1,3),
						v = zeros(1,3))
	return NOMINAL_LANDING(z_init, vx_init, vz_init, θ, t_extra, sampler, dt, 
							vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::NOMINAL_LANDING)
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
				   vx[N] == 0; pz[N] == 0; vz[N] == 0;
				   vx ≥ 0; pz ≥ 0; vz ≤ 0; vz ≥ -τ.vz_max;
				   D*vx ≤ τ.ax_max; D*vz*ft2m ≤ τ.az_max]

	prob = minimize(objec, constraints)
	solve!(prob, ECOS.Optimizer(verbose=0))

	# Update the trajectory
	p = zeros(N,3)
	p[:,1] = px.value[:,1]
	p[:,3] = pz.value[:,1]
	v = zeros(N,3)
	v[:,1] = vx.value[:,1]
	v[:,3] = vz.value[:,1]

	τ.p = copy(p)
	τ.v = copy(v)

	return Int(prob.status) == 1
end