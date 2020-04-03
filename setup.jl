using Pkg

if !in("Distributions", keys(Pkg.installed()))
	Pkg.add("Distributions")
end

if !in("Convex", keys(Pkg.installed()))
	Pkg.add("Convex")
end

if !in("ECOS", keys(Pkg.installed()))
	Pkg.add("ECOS")
end

if !in("DataFrames", keys(Pkg.installed()))
	Pkg.add("DataFrames")
end

if !in("CSV", keys(Pkg.installed()))
	Pkg.add("CSV")
end