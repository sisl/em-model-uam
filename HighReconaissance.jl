"""
---------------------------------------------------
Types and Constructors
---------------------------------------------------
"""

mutable struct HIGH_RECONAISSANCE_SAMPLER <: SAMPLER
	z::Sampleable
	vx::Sampleable
	radius::Sampleable
	arc_length::Sampleable
end

function high_reconaissance_sampler(;z = Normal(500,10),
						vx = Uniform(30,45),
						radius = Uniform(200,400),
						arc_length = Uniform(90,270))
	return HIGH_RECONAISSANCE_SAMPLER(z, vx, radius, arc_length)
end

mutable struct HIGH_RECONAISSANCE <: UAM_TRAJECTORY
	z::Float64
	vx::Float64
	radius::Float64
	arc_length::Float64
	sampler::SAMPLER
	dt::Float64
	vz_max::Float64
	az_max::Float64
	ax_max::Float64
	p::Array{Float64,2}
	v::Array{Float64,2}
end

function high_reconaissance(;z = 0.0,
						vx = 0.0,
						radius = 0.0,
						arc_length = 0.0,
						sampler = high_reconaissance_sampler(),
						dt = 1.0,
						vz_max = 550fpm2fps, # This doesn't actually get used
						az_max = 0.3g, # This doesn't actually get used
						ax_max = 0.1g, # This doesn't actually get used
						p = zeros(1,3),
						v = zeros(1,3))
	return HIGH_RECONAISSANCE(z, vx, radius, arc_length, sampler, dt, 
								vz_max, az_max, ax_max, p, v)
end

"""
---------------------------------------------------
Functions
---------------------------------------------------
"""

function solve_trajectory!(τ::HIGH_RECONAISSANCE)
	# Get the circular portion first
	θ̇ = rad2deg(τ.vx/τ.radius)
	arc_time = τ.arc_length/θ̇
	arc_steps = convert(Int64, floor(arc_time/τ.dt))
	θ = zeros(arc_steps)
	for i = 2:arc_steps
		θ[i] = θ[i-1] + θ̇*τ.dt
	end

	x = τ.radius*cosd.(θ)
	y = τ.radius*sind.(θ)
	z = τ.z*ones(arc_steps)

	p = hcat(x, y, z)

	N = size(p,1)
	D = zeros(N-1,N)
	for i = 1:(N-1)
	    D[i,i] = -1
	    D[i,i+1] = 1
	end
	D = D./τ.dt

	v = D*p
	v = vcat(v, v[end,:]')

	τ.p = copy(p)
	τ.v = copy(v)

	return true # Because others return if they are optimal
end

	

