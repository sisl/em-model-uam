"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct VERTICAL_ASCENT_SAMPLER <: SAMPLER
	z_final::Sampleable
	z_ascent::Sampleable
	vx_final::Sampleable
	θ::Sampleable
	t_extra::Sampleable
	t_ascent::Sampleable
end

function vertical_ascent_sampler(;z_final = Normal(500,20),
						z_ascent = Normal(100,8),
						vx_final = Uniform(30,45),
						θ = Uniform(3,5),
						t_extra = Uniform(0.1,0.3), #Normal(10,4),
						t_ascent = Normal(30,5))
	return VERTICAL_ASCENT_SAMPLER(z_final, z_ascent, vx_final, θ, t_extra, t_ascent)
end

mutable struct VERTICAL_ASCENT <: UAM_TRAJECTORY
	z_final::Float64
	z_ascent::Float64
	vx_final::Float64
	θ::Float64
	t_extra::Float64
	t_ascent::Float64
	sampler::SAMPLER
	dt::Float64
	vz_max::Float64
	az_max::Float64
	ax_max::Float64
	p::Array{Float64,2}
	v::Array{Float64,2}
end

function vertical_ascent(;z_final = 0.0,
						z_ascent = 0.0,
						vx_final = 0.0,
						θ = 0.0,
						t_extra = 0.0,
						t_ascent = 0.0,
						sampler = vertical_ascent_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps,
						az_max = 0.3g_ft,
						ax_max = 0.1g,
						p = zeros(1,3),
						v = zeros(1,3))
	return VERTICAL_ASCENT(z_final, z_ascent, vx_final, θ, t_extra, t_ascent, sampler, dt, 
							vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::VERTICAL_ASCENT)
	# Vertical ascent portion ######################################
	Nd = convert(Int64, floor(τ.t_ascent/τ.dt))
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
				   vzd[1] == 0; pzd[1] == 0;
				   pzd[Nd] == τ.z_ascent; vzd[Nd] == 0;
				   pzd ≥ 0; vzd ≥ 0; vzd ≤ τ.vz_max;
				   Dd*vzd*ft2m ≥ -τ.az_max; Dd*vzd*ft2m ≤ τ.az_max]

	prob2 = minimize(objec, constraints)
	solve!(prob2, ECOS.Optimizer(verbose=0))

	# Glide portion #################################################
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
				   vz[1] == 0; pz[1] == τ.z_ascent;
				   vx[N] == τ.vx_final; pz[N] == τ.z_final; vz[N] == 0;
				   vx ≥ 0; pz ≥ 0; vz ≥ 0; vz ≤ τ.vz_max;
				   D*vx ≤ τ.ax_max; D*vz*ft2m ≤ τ.az_max;
				   D*vx ≥ -τ.ax_max; D*vz*ft2m ≥ -τ.az_max]

	prob1 = minimize(objec, constraints)
	solve!(prob1, ECOS.Optimizer(verbose=0))

	# Update the trajectory
	p = zeros(Nd+N,3)
	p[Nd+1:end,1] = px.value[:,1]
	p[1:Nd,3] = pzd.value[:,1]
	p[Nd+1:end,3] = pz.value
	v = zeros(Nd+N,3)
	v[Nd+1:end,1] = vx.value[:,1]
	v[1:Nd,3] = vzd.value[:,1]
	v[Nd+1:end,3] = vz.value

	τ.p = copy(p)
	τ.v = copy(v)

	return Int(prob1.status) == 1
end