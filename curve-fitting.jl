### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 393b2f5e-6556-11eb-2119-cf7309ee7392
begin
	using Chain
	using CSV
	using DataFrames
	using LsqFit
	using Measurements
	using Plots
	using PlutoUI
	using Statistics
end

# ╔═╡ 412eedac-658a-11eb-2326-93e5bf3d1a2c
md"# Non-linear curve fitting

This notebook plots a dataset (with error bars for datasets containing replicates) and performs non-linear curve fitting of a model to the data."

# ╔═╡ 227f554a-658a-11eb-0a12-41577ee98192
md"## Load raw data and process into mean ± std"

# ╔═╡ 270cc0cc-660f-11eb-241e-b75746a39cc7
md"Indicate which data file to process (path is relative to the notebook file):"

# ╔═╡ 14baf100-660f-11eb-1380-ebf4a860eed8
dataFile = "datasets/dataset_003.csv"
# Try to also make it work with an array of several file names.
# datasets/dataset_005.csv is a good one to plot with 003.

# ╔═╡ 4233708c-6557-11eb-0581-4f79b924af56
df = @chain dataFile begin
	CSV.File(footerskip=2)
	DataFrame()
	transform(
		# Rename column 1, so we can always call it by the
		# same name regardless of its name in the input file.
		1 => :concentration,
		copycols = false
	)
	@aside cols = ncol(_)
	transform(
		# Calculate mean and stddev of replicates
		# (all columns in input except first one).
		AsTable(2:cols) => ByRow(mean) => :mean,
		AsTable(2:cols) => ByRow(std) => :std
	)
	transform(
		# Mean and stddev together define a measurement
		# (this is only for plotting; fitting uses the two
		# original columns separately).
		[:mean, :std] => ByRow(measurement) => :measurement
	)
end

# ╔═╡ 5c5e0392-658a-11eb-35be-3d940d4504cb
md"## Visualizations"

# ╔═╡ 01b59d8a-6637-11eb-0da0-8d3e314e23af
md"Indicate the receptor concentration (in an FP experiment, this is the concentration of fluorescently labeled probe). This parameter only appears in the quadratic model."

# ╔═╡ 171fcd26-6637-11eb-0deb-3dac3dd858b4
@bind R0 PlutoUI.Slider(0.1:0.1:20.0, default = 5.0, show_value = true)

# ╔═╡ 663a4cae-658a-11eb-382f-cf256c08c9d1
md"## Model functions"

# ╔═╡ 88d941e6-658a-11eb-08a2-0f021e5ae3a4
md"""This is the Hill equation:

$S = S_{min} + (S_{max} - S_{min}) \times \frac{L^h}{{K_D}^h + L^h}$

In which $S$ is the measured signal at a given value of ligand concentration $L$, $S_{min}$ and $S_{max}$ are the minimum and maximum values the observed signal can take, respectively, $K_D$ is the equilibrium dissociation constant and $h$ is the Hill coefficient."""

# ╔═╡ e9fc3d44-6559-11eb-2da7-314e8fc76ee9
@. hill(conc, p) = p[1] + (p[2] - p[1]) * conc^p[4] / (p[3]^p[4] + conc^p[4])

# ╔═╡ 9d1f24cc-6a0f-11eb-3b16-35f89aff5d4a
md"The hyperbolic model is a special case, where $h = 1$:"

# ╔═╡ 78e664d0-6618-11eb-135b-5574bb05ddef
@. hyperbolic(conc, p) = p[1] + (p[2] - p[1]) * conc / (p[3] + conc)

# ╔═╡ b0b17206-6a0f-11eb-2f5e-5fc8fa06cd36
md"""Unlike the Hill and hyperbolic models, the quadratic model does not make the approximation that the concentration of free ligand at equilibrium is equal to the total ligand concentration:

$S = S_{min} + (S_{max} - S_{min}) \times \frac{(K_{D} + R_{tot} + L_{tot}) - \sqrt{(- K_{D} - R_{tot} - L_{tot})^2 - 4 \times R_{tot} \times L_{tot}}}{2 \times R_{tot}}$

Symbols have the same meaning as in the previous equations, except here $L_{tot}$ is the total concentration of ligand, and the parameter $R_{tot}$ is the total concentration of receptor.

$R_{tot}$ could in principle be left as a free parameter to be determined by the fitting procedure, but in general it is known accurately enough from the experimental set up, and one should replicate the same experiment with different concentrations of receptor to check its effect on the results."""

# ╔═╡ 5694f1da-6636-11eb-0fed-9fee5c48b114
@. quadratic(conc, p) = p[1] + (p[2] - p[1]) * ( (p[3] + R0 + conc) - sqrt((-(p[3] + R0 + conc)) ^ 2 - 4 * R0 * conc) ) / (2 * R0)

# ╔═╡ 0f960a8a-6a0f-11eb-04e2-b543192f6354
md"## Parameters and their initial values"

# ╔═╡ 2eb890a4-658c-11eb-1fc4-af645d74109d
md"The following arrays store initial values for the model parameters (in this order): $S_{min}$, $S_{max}$, $K_D$ and $h$ (for the Hill model only).

Initial values for $S_{min}$ and $S_{max}$ are simply taken as the minimal and maximal values found in the data. The initial estimate for $K_D$ is the concentration of the data point that has a signal closest to halfway between $S_{min}$ and $S_{max}$ (if the experiment was properly designed, this is a reasonable estimate and close enough to the true value for the fit to converge). The initial estimate of $h$ is $1.0$, meaning we assume no cooperativity."

# ╔═╡ 90083670-655a-11eb-13e3-ddbc86ab2ac8
hillParams = [
	minimum(df.mean),
	maximum(df.mean),
	df.concentration[
		findmin(abs.((maximum(df.mean) - minimum(df.mean) / 2) .- df.mean))[2]
	],
	1.0
]

# ╔═╡ 3c5132f6-6619-11eb-285a-03d30ac28819
hyperbolicParams = [
	minimum(df.mean),
	maximum(df.mean),
	df.concentration[
		findmin(abs.((maximum(df.mean) - minimum(df.mean) / 2) .- df.mean))[2]
	]
]

# ╔═╡ 852d017c-6636-11eb-0112-2f7621211717
quadraticParams = [
	minimum(df.mean),
	maximum(df.mean),
	df.concentration[
		findmin(abs.((maximum(df.mean) - minimum(df.mean) / 2) .- df.mean))[2]
	]
]

# ╔═╡ 213e8ffa-6a0f-11eb-357e-638146193c5d
md"## Fitting"

# ╔═╡ c91ec0aa-655a-11eb-1916-a70d97224aeb
fitHill = curve_fit(hill, df.concentration, df.mean, df.std, hillParams)

# ╔═╡ 1f0384de-659b-11eb-043e-5b86fcdd36e6
dof(fitHill)

# ╔═╡ a74998b4-659c-11eb-354d-09ff62710b87
stderror(fitHill)

# ╔═╡ e4b34962-6618-11eb-1847-035583a15e8e
fitHyp = curve_fit(hyperbolic, df.concentration, df.mean, df.std, hyperbolicParams)

# ╔═╡ f23da6b8-6618-11eb-1718-7971ade4c485
dof(fitHyp)

# ╔═╡ f66336fe-6618-11eb-0406-a90c1a1359da
stderror(fitHyp)

# ╔═╡ 26f4a4f6-6637-11eb-32eb-8d8bf82667ea
fitQuad = curve_fit(quadratic, df.concentration, df.mean, df.std, quadraticParams)

# ╔═╡ e27cb090-6558-11eb-0fae-178b14e7fa8c
begin
	scatter(
		df.concentration,
		df.measurement,
		xscale = :log10,
		label = "Data"
	)
	plot!(
		df.concentration,
		hill(df.concentration, fitHill.param),
		label = "Hill fit"
	)
	plot!(df.concentration,
		hyperbolic(df.concentration, fitHyp.param),
		label = "Hyperbolic fit"
	)
	plot!(df.concentration,
		quadratic(df.concentration, fitQuad.param),
		label = "Quadratic fit (R0 = $R0)"
	)
	plot!(
		title = dataFile,
		xlabel = "Concentration",
		ylabel = "Signal",
		legend = :topleft
	)
end

# ╔═╡ a5f01744-655c-11eb-1248-23eaa89fcf09
begin
	plot(
		title = dataFile,
		xlabel = "Concentration",
		ylabel = "Fit residual",
		legend = :topleft
	)
	scatter!(
		df.concentration,
		fitHill.resid,
		xscale = :log10,
		label = "Hill fit"
	)
	scatter!(
		df.concentration,
		fitHyp.resid,
		label = "Hyperbolic fit"
	)
	scatter!(
		df.concentration,
		fitQuad.resid,
		label = "Quadratic fit (R0 = $R0)"
	)
	
end

# ╔═╡ be17b97e-663a-11eb-2158-a381c19ece3f
md"## Results

Hill:

- Kd = $(round(fitHill.param[3])) ± $(round(stderror(fitHill)[3]))
- h = $(round(fitHill.param[4])) ± $(round(stderror(fitHill)[4]))

Hyperbolic: Kd = $(round(fitHyp.param[3])) ± $(round(stderror(fitHyp)[3]))

Quadratic: Kd = $(round(fitQuad.param[3])) ± $(round(stderror(fitQuad)[3]))"

# ╔═╡ 88d13a40-6637-11eb-2963-3be613d0e154
dof(fitQuad)

# ╔═╡ 8e71d808-6637-11eb-24d4-cf8862fe58b9
stderror(fitQuad)

# ╔═╡ Cell order:
# ╟─412eedac-658a-11eb-2326-93e5bf3d1a2c
# ╠═393b2f5e-6556-11eb-2119-cf7309ee7392
# ╟─227f554a-658a-11eb-0a12-41577ee98192
# ╟─270cc0cc-660f-11eb-241e-b75746a39cc7
# ╠═14baf100-660f-11eb-1380-ebf4a860eed8
# ╠═4233708c-6557-11eb-0581-4f79b924af56
# ╟─5c5e0392-658a-11eb-35be-3d940d4504cb
# ╟─01b59d8a-6637-11eb-0da0-8d3e314e23af
# ╠═171fcd26-6637-11eb-0deb-3dac3dd858b4
# ╠═e27cb090-6558-11eb-0fae-178b14e7fa8c
# ╠═a5f01744-655c-11eb-1248-23eaa89fcf09
# ╟─be17b97e-663a-11eb-2158-a381c19ece3f
# ╟─663a4cae-658a-11eb-382f-cf256c08c9d1
# ╟─88d941e6-658a-11eb-08a2-0f021e5ae3a4
# ╠═e9fc3d44-6559-11eb-2da7-314e8fc76ee9
# ╟─9d1f24cc-6a0f-11eb-3b16-35f89aff5d4a
# ╠═78e664d0-6618-11eb-135b-5574bb05ddef
# ╟─b0b17206-6a0f-11eb-2f5e-5fc8fa06cd36
# ╠═5694f1da-6636-11eb-0fed-9fee5c48b114
# ╟─0f960a8a-6a0f-11eb-04e2-b543192f6354
# ╟─2eb890a4-658c-11eb-1fc4-af645d74109d
# ╠═90083670-655a-11eb-13e3-ddbc86ab2ac8
# ╠═3c5132f6-6619-11eb-285a-03d30ac28819
# ╠═852d017c-6636-11eb-0112-2f7621211717
# ╟─213e8ffa-6a0f-11eb-357e-638146193c5d
# ╠═c91ec0aa-655a-11eb-1916-a70d97224aeb
# ╠═1f0384de-659b-11eb-043e-5b86fcdd36e6
# ╠═a74998b4-659c-11eb-354d-09ff62710b87
# ╠═e4b34962-6618-11eb-1847-035583a15e8e
# ╠═f23da6b8-6618-11eb-1718-7971ade4c485
# ╠═f66336fe-6618-11eb-0406-a90c1a1359da
# ╠═26f4a4f6-6637-11eb-32eb-8d8bf82667ea
# ╠═88d13a40-6637-11eb-2963-3be613d0e154
# ╠═8e71d808-6637-11eb-24d4-cf8862fe58b9
