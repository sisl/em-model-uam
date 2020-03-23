"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct NOMINAL_TAKEOFF_SAMPLER <: SAMPLER
	z_final::Sampleable
	vx_final::Sampleable
	θ::Sampleable
	t_extra::Sampleable
end

function nominal_takeoff_sampler(;z_final = Normal(500,20),
						vx_final = Uniform(30,45),
						θ = Uniform(3,5),
						t_extra = Uniform(0.1,0.3)) #Normal(10,4))
	return NOMINAL_TAKEOFF_SAMPLER(z_final, vx_final, θ, t_extra)
end

mutable struct NOMINAL_TAKEOFF <: UAM_TRAJECTORY
	z_final::Float64
	vx_final::Float64
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

function nominal_takeoff(;z_final = 0.0,
						vx_final = 0.0,
						θ = 0.0,
						t_extra = 0.0,
						sampler = nominal_takeoff_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps,
						az_max = 0.3g_ft,
						ax_max = 0.1g,
						p = zeros(1,3),
						v = zeros(1,3))
	return NOMINAL_TAKEOFF(z_final, vx_final, θ, t_extra, sampler, dt, 
							vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::NOMINAL_TAKEOFF)
	# Get some needed info
	#ttot = abs(τ.z_final/500fpm2fps) + τ.t_extra*abs(τ.z_final/500fpm2fps)
	ttot = abs(τ.z_final/τ.vz_max) + τ.t_extra*abs(τ.z_final/τ.vz_max)
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
				   vx[1] == 0; px[1] == 0;
				   vz[1] == 0; pz[1] == 0;
				   vx[N] == τ.vx_final; pz[N] == τ.z_final; vz[N] == 0;
				   vx ≥ 0; pz ≥ 0; vz ≥ 0; vz ≤ τ.vz_max;
				   D*vx ≤ τ.ax_max; D*vz*ft2m ≤ τ.az_max;
				   D*vx ≥ -τ.ax_max; D*vz*ft2m ≥ -τ.az_max]

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