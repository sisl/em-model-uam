"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct CURVED_APPROACH_SAMPLER <: SAMPLER
	z_init::Sampleable
	vx_init::Sampleable
	vz_init::Sampleable
	radius::Sampleable
	t_extra::Sampleable
end

function curved_approach_sampler(;z_init = Uniform(300,600),
						vx_init = Uniform(30,45),
						vz_init = Uniform(-500fpm2fps,-300fpm2fps),
						radius = Uniform(200,400),
						t_extra = Normal(20,4))
	return CURVED_APPROACH_SAMPLER(z_init, vx_init, vz_init, radius, t_extra)
end

mutable struct CURVED_APPROACH <: UAM_TRAJECTORY
	z_init::Float64
	vx_init::Float64
	vz_init::Float64
	radius::Float64
	t_extra::Float64
	sampler::SAMPLER
	dt::Float64
	vz_max::Float64
	az_max::Float64
	ax_max::Float64
	p::Array{Float64,2}
	v::Array{Float64,2}
end

function curved_approach(;z_init = 0.0,
						vx_init = 0.0,
						vz_init = 0.0,
						radius = 0.0,
						t_extra = 0.0,
						sampler = curved_approach_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps, # this doesn't actually get used right now
						az_max = 0.3g_ft,
						ax_max = 0.1g, # this doesn't actually get used right now
						p = zeros(1,3),
						v = zeros(1,3))
	return CURVED_APPROACH(z_init, vx_init, vz_init, radius, t_extra, sampler, dt, 
							vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::CURVED_APPROACH)
	# Get some needed info
	ttot = abs(τ.z_init/τ.vz_init) + τ.t_extra
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
	py = Variable(N)
	pz = Variable(N)
	vx = Variable(N)
	vy = Variable(N)
	vz = Variable(N)

	objec = 100*sumsquares(D*vx) + 100*sumsquares(D*vy) + sumsquares(D*vz*ft2m)

	constraints = [[px[i+1] == px[i] + vx[i]*τ.dt for i = 1:N-1];
				   [py[i+1] == py[i] + vy[i]*τ.dt for i = 1:N-1];
				   [pz[i+1] == pz[i] + vz[i]*τ.dt for i = 1:N-1];
				   px[1] == 0; py[1] == 2*τ.radius; pz[1] == τ.z_init;
				   vx[1] == -τ.vx_init; vy[1] == 0; vz[1] == τ.vz_init;
				   px[end] == 0; py[end] == 0; pz[end] == 0;
				   vx[end] == 0.001; vy[end] == 0; vz[end] == 0; # Giving small x speed because I want it to constrain heading
				   D*vz*ft2m ≥ -τ.az_max; D*vz*ft2m ≤ τ.az_max]
   
	prob = minimize(objec, constraints)
	solve!(prob, ECOS.Optimizer(verbose=0))

	# Update the trajectory
	p = zeros(N,3)
	p[:,1] = px.value[:,1]
	p[:,2] = py.value[:,1]
	p[:,3] = pz.value[:,1]
	v = zeros(N,3)
	v[:,1] = vx.value[:,1]
	v[:,2] = vy.value[:,1]
	v[:,3] = vz.value[:,1]

	τ.p = copy(p)
	τ.v = copy(v)

	return Int(prob.status) == 1
end