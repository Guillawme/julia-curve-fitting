### A Pluto.jl notebook ###
# v0.19.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 393b2f5e-6556-11eb-2119-cf7309ee7392
begin
	using Chain
	using CSV
	using DataFrames
	using HTTP
	using LsqFit
	using Measurements
	using Plots
	using PlutoUI
	using Statistics
	using Unitful
	using URIs
end

# ╔═╡ 412eedac-658a-11eb-2326-93e5bf3d1a2c
md"""
# Fitting of equilibrium binding data

This notebook plots an equilibrium binding dataset (with error bars, for datasets containing replicates) and performs non-linear curve fitting of a model to the data.

The following reference explains very well the theory of equilibrium binding experiments, as well as many important practical considerations:

> Jarmoskaite I, AlSadhan I, Vaidyanathan PP & Herschlag D (2020) How to measure and evaluate binding affinities. *eLife* **9**: e57264 <https://doi.org/10.7554/eLife.57264>
"""

# ╔═╡ 270cc0cc-660f-11eb-241e-b75746a39cc7
md"""
## Load data

The data file must be in CSV format. The first row is assumed to contain column names. The first column is assumed to be the $X$ values, all other columns are assumed to be replicate $Y$ values that will be averaged (fitting will be done against the mean values, weighted by their standard deviations). In addition, there must not be any row with an $X = 0$ value (this would result in an error when attempting to plot with a logarithmic scale). I always perform two measurements at $X = 0$, and I always sort rows by descending $X$ values, so this notebook automatically skips the last two rows of the CSV file; adjust accordingly if you don't measure at $X = 0$ and want to keep all rows (see section [Data processing](#1884912a-6aeb-11eb-2b4a-d14d4a321dc5) below). The data does *not* need to be scaled such that $Y$ takes values between $0$ and $1$: the binding models can account for arbitrary minimum and maximum $Y$ values (see section [Model functions](#663a4cae-658a-11eb-382f-cf256c08c9d1) below). Scaling data will hide differences in signal change between datasets, while these differences may tell you something about the system under study, so scaling should never be done "blindly"; always look at the raw data.

Indicate below which data files to process:

- if the path given is not absolute, it is assumed to be relative to the notebook file (wherever the notebook is located)
- a list of files can be provided with one file path per line and separated by commas
- files can be located by local path or URL, but each type of location should be in the dedicated list
"""

# ╔═╡ a2c02fcf-9402-4938-bc3d-819b45a66afa
dataURLs = [
	"https://raw.githubusercontent.com/Guillawme/julia-curve-fitting/main/datasets/dataset_006.csv",
	"https://raw.githubusercontent.com/Guillawme/julia-curve-fitting/main/datasets/dataset_007.csv"
]

# ╔═╡ 14baf100-660f-11eb-1380-ebf4a860eed8
dataFiles = [
	"datasets/dataset_008.csv",
	"datasets/dataset_009.csv"
]

# ╔═╡ 2ce72e97-0133-4f15-bf1d-7fd04ccf3102
md"""
**Number of rows to ignore at the end of files:** $(@bind footerRows PlutoUI.NumberField(0:30, default = 2))
"""

# ╔═╡ 77316a92-425a-4902-9828-52a7c4a74f27
md"""
**Unit of your concentration values ($X$ axis):** $(
@bind chosenConcUnit PlutoUI.Select(
	[
		"pM" => "pM",
		"nM" => "nM",
		"μM" => "μM",
		"mM" => "mM"
	],
	default = "nM"
))
"""

# ╔═╡ 214acce6-6ae8-11eb-3abf-492e50140317
md"""
Your data should appear below shortly, check that it looks normal. In addition to the columns present in your CSV file, you should see four columns named `mean`, `std`, `measurement` and `conc` (these values will be used for fitting and plotting).
"""

# ╔═╡ d5b0e2a1-865c-489c-9d0d-c4ae043828fb
# By defaults, use file names and URLs to identify datasets.
datasetNames = vcat(dataURLs, dataFiles)
# But one can also use custom names.
#datasetNames = ["Monday", "Tuesday", "Wednesday", "Thursday"]

# ╔═╡ 5c5e0392-658a-11eb-35be-3d940d4504cb
md"""
## Visualizations

Your data and fit should appear below shortly. Take a good look at the [data and fit](#3dd72c58-6b2c-11eb-210f-0b547bf38ebe), make sure you check the [residuals](#4f4000b4-6b2c-11eb-015f-d76a0adda0a0). Once you're happy with it, check the [numerical results](#be17b97e-663a-11eb-2158-a381c19ece3f).
"""

# ╔═╡ 3dd72c58-6b2c-11eb-210f-0b547bf38ebe
md"""
### Data and fit

Select binding model:
"""

# ╔═╡ 3da83f72-6a11-11eb-1a74-49b66eb39c96
@bind chosenModel PlutoUI.Radio(
	[
		"Hill" => :Hill,
		"Hyperbolic" => :Hyperbolic,
		"Quadratic" => :Quadratic
	],
	default = "Hill"
)

# ╔═╡ d15cba72-6aeb-11eb-2c80-65702b48e859
md"""
Show fit line with initial parameters?
$@bind showInitialFit PlutoUI.CheckBox(default = false)
"""

# ╔═╡ 01b59d8a-6637-11eb-0da0-8d3e314e23af
md"""
For the quadratic model, indicate receptor concentration (the receptor is the binding partner kept at constant, low concentration across the titration series).
Parameter $R_0 =$
$@bind R0 PlutoUI.Slider(0.01:0.1:500.0, default = 5.0, show_value = true)
"""

# ╔═╡ 4f4000b4-6b2c-11eb-015f-d76a0adda0a0
md"""
### Residuals
"""

# ╔═╡ c50cf18c-6b11-11eb-07d3-0b8e332ec5bc
md"""
The fit residuals should follow a random normal distribution around $0$. If they show a systematic trend, it means the fit systematically deviates from your data, and therefore the model you chose might not be justified (but be careful when considering alternative models: introducing more free parameters will likely get the fit line closer to the data points and yield a lower [sum of squared residuals](#124c4f94-6b99-11eb-2921-d7c2cd00b893), but this is not helpful if these additional parameters don't contribute to explaining the physical phenomenon being modeled). Another possibility is a problem with your data. The most common problems are:

- the data does not cover the proper concentration range
- the concentration of receptor is too high relative to the $K_D$

In either case, your best option is to design a new experiment and collect new data.
"""

# ╔═╡ 5a36fc3f-ce74-42c8-8284-19321e0d687f
md"""
#### Scatter plot
"""

# ╔═╡ 5392d99b-70f9-48cd-90b4-58cba5fc9681
md"""
#### Histogram
"""

# ╔═╡ be17b97e-663a-11eb-2158-a381c19ece3f
md"""
## Numerical results

### Model parameters

The parameter $K_D$ is in unit of $(chosenConcUnit).
"""

# ╔═╡ 124c4f94-6b99-11eb-2921-d7c2cd00b893
md"""
### Sum of squared residuals
"""

# ╔═╡ 7e7a9dc4-6ae8-11eb-128d-83544f01b78b
md"""
## Code

The code doing the actual work is in this section. Do not edit unless you know what you are doing.
"""

# ╔═╡ 512e3028-6ae9-11eb-31b4-1bc9fc66b322
md"### Necessary packages and notebook setup"

# ╔═╡ abc03f64-6a11-11eb-0319-ed7cea455cb5
PlutoUI.TableOfContents()

# ╔═╡ 1884912a-6aeb-11eb-2b4a-d14d4a321dc5
md"### Data processing"

# ╔═╡ c94ef72d-bd12-4434-935e-01e94d5a4588
md"""
#### Concentration unit selection
"""

# ╔═╡ 54272681-dfeb-4034-b127-7e68c19fd576
md"""
First, we need an alias of M (molar) for mol/L, since this is the notation most widely used in the field:
"""

# ╔═╡ 731492c6-95c7-449f-8c19-53e22ab438b8
begin
	Unitful.register(@__MODULE__)
	@unit M "M" Molar 1u"mol/L" true
end

# ╔═╡ 93ef0431-643b-4f6f-8c09-beabff59e0c6
md"""
Check that our alias works:
"""

# ╔═╡ 1c10231d-3bea-4468-a2d8-886c05c6474c
typeof(M)

# ╔═╡ 730b693d-cca8-46de-8382-c151b5f63352
1u"nM" == 1u"nmol/l"

# ╔═╡ 24152a84-523b-4027-9b5a-7e3524b9c659
dimension(1u"nM") == dimension(1u"mol/l")

# ╔═╡ 97f73ee4-3db7-43ba-93eb-17025b485f4f
md"""
This dictionary maps dropdown menu options (in section [Load data](#270cc0cc-660f-11eb-241e-b75746a39cc7) above) to their corresponding unit:
"""

# ╔═╡ 8f2d959e-5e61-48a9-bbf0-538bc1d478a8
concUnits = Dict(
	"pM" => u"pM",
	"nM" => u"nM",
	"μM" => u"μM",
	"mM" => u"mM"
)

# ╔═╡ c866d213-4680-4b08-8ecc-1faefc8661a4
md"""
The following cells simply check which unit is selected in section [Load data](#008f4f8c-9d21-11eb-0fea-f3b2e58957d1) above.
"""

# ╔═╡ 8e9cd30f-1722-4f2f-a26b-2f558805d4a1
chosenConcUnit

# ╔═╡ 13417555-be95-4662-ad11-eca4dece81b5
concUnits[chosenConcUnit]

# ╔═╡ 9c7e922e-c82f-41a4-9513-4462a0559c3f
md"""
#### Processing
"""

# ╔═╡ 4f4b580d-507c-4ad0-b1d5-5967c8ed829e
md"""
The `commonProcessing()` function computes the mean and standard deviation of replicates, defines measurements as mean ± std, and returns a DataFrame containing all the data. It is used by all methods of the following `processData()` function, which handle various inputs (path to a local file, URL to a remote file, loaded CSV file, loaded DataFrame).
"""

# ╔═╡ 5eb607c7-172b-4a6c-a815-acbc195108f0
function commonProcessing(data::DataFrame)
	df = @chain data begin
		# Rename column 1, so we can always call it by the
		# same name regardless of its name in the input file.
		rename(1 => :concentration)
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
		transform(
			# Assign unit to concentration column, store as
			# a new column callec conc (only for plotting).
			:concentration => ByRow(x -> x * concUnits[chosenConcUnit]) => :conc
		)
	end
	return df
end

# ╔═╡ 992c86a2-6b13-11eb-1e00-95bdff2736d0
md"""
The `processData()` function loads one data file, computes the mean and standard deviation of replicates, defines measurements as mean ± std, and returns a DataFrame containing all the data.
"""

# ╔═╡ 1fe4d112-6a11-11eb-37a6-bf95fbe032b1
function processData(dataFile::String)
	df = @chain dataFile begin
		CSV.read(footerskip=footerRows, DataFrame)
		commonProcessing()
	end
	return df
end

# ╔═╡ 5cc812af-4f8e-444b-9652-fb063cc6be06
md"""
This functions should also work if passed an URL to a remote file:
"""

# ╔═╡ 7c2af794-a0ed-4827-bc8e-1f15ad205eca
function processData(dataFile::URI)
	df = @chain dataFile begin
		HTTP.get(_).body
		CSV.read(footerskip=footerRows, DataFrame)
		commonProcessing()
	end
	return df
end

# ╔═╡ a08d1b29-fa5e-4312-977c-39b050de4516
md"""
This functions should also work if passed an already loaded CSV file:
"""

# ╔═╡ 0181a0b1-ec6c-4325-8b3a-851a3fe33846
function processData(dataFile::CSV.File)
	df = @chain dataFile begin
		DataFrame()
		commonProcessing()
	end
	return df
end

# ╔═╡ 247c2416-6c67-11eb-01df-8dac01cdbf8f
md"""
This functions should also work if passed an already loaded data frame (for example, if the user wants to load data and pre-process it in a different way before averaging replicates):
"""

# ╔═╡ 36ebe112-6c66-11eb-11f0-7fdb946865e4
function processData(data::DataFrame)
	return commonProcessing(data)
end

# ╔═╡ d904fd76-6af1-11eb-2352-837e03072137
begin
	allData = vcat(URI.(dataURLs), dataFiles)
	dfs = [ processData(df) for df in allData ]
end

# ╔═╡ 0e8af3be-7ae7-4ec2-8d7a-670878cd52ee
md"""
We will need to keep track of dataset names. You can edit them here if you want to use something more meaningful than the file name or URL; changes will propagate to plot legends. This list **must** contain the same number of elements as you have datasets: **$(length(allData))** in the present case.
"""

# ╔═╡ 8e105fae-6aec-11eb-1471-83aebb776241
md"### Plotting"

# ╔═╡ 97c5020c-6aec-11eb-024b-513b1e603d98
md"The `initMainPlot()` function initializes a plot, the `plotOneDataset!()` function plots one dataset (call it repeatedly to plot more datasets on the same axes)."

# ╔═╡ caf4a4a2-6aec-11eb-2765-49d67afa47dd
function initMainPlot()
	plot(
		xlabel = "Concentration",
		ylabel = "Signal",
		xscale = :log10,
		legend = :topleft
	)
end

# ╔═╡ fc194672-6aed-11eb-0a06-2d967ec094b1
md"The `initResidualPlot()` function initializes a plot, the `plotOneResiduals!()` function plots the fit residuals from one dataset (call it repeatedly to plot more datasets on the same axes)."

# ╔═╡ 14db987c-6aee-11eb-06cf-a11987b98f1e
function initResidualPlot()
	plot(
		xlabel = "Concentration",
		ylabel = "Fit residual",
		xscale = :log10,
		legend = :topleft
	)
	hline!([0], label = nothing, color = :red)
end

# ╔═╡ 7f83b838-6a11-11eb-3652-bdff24f3473e
function plotOneResiduals!(plt, df, fit, filePath)
	title = split(filePath, "/")[end]
	scatter!(
		plt,
		df.concentration,
		fit.resid,
		label = "$title: $chosenModel fit residual"
	)
end

# ╔═╡ 9020fa5d-7408-4161-a52f-df37b3c2e6f5
md"The `initResidualHistogram()` function initializes a histogram, the `plotOneResidualsHistogram!()` function plots a histogram of the fit residuals from one dataset (call it repeatedly to plot more datasets on the same axes)."

# ╔═╡ 3184d209-1cc9-40ed-a9ec-9f094c5c94b5
function initResidualHistogram()
	histogram(
		xlabel = "Fit residual",
		ylabel = "Count",
		legend = :topleft
	)
	vline!([0], label = nothing, color = :red)
end

# ╔═╡ b38cd229-64e2-4ca4-a78b-8881ec166b09
function plotOneResidualsHistogram!(plt, df, fit, filePath)
	title = split(filePath, "/")[end]
	histogram!(
		plt,
		fit.resid,
		bins = length(fit.resid),
		label = "$title: $chosenModel fit residual"
	)
end

# ╔═╡ 663a4cae-658a-11eb-382f-cf256c08c9d1
md"### Model functions"

# ╔═╡ 594e7534-6aeb-11eb-1254-3b92b71877ed
md"#### Model selection"

# ╔═╡ 32dce844-6aee-11eb-3cf2-3ba420d311d3
md"""
This dictionary maps radio button options (in section [Visualizations](#5c5e0392-658a-11eb-35be-3d940d4504cb) above) to their corresponding model function:
"""

# ╔═╡ a1f56b0a-6aeb-11eb-0a44-556fad58f368
md"""
The remaining cells in this section are only meant to check that the model selection buttons work. This first cell should return the name of the selected binding model (corresponding to the active radio button in section [Visualizations](#5c5e0392-658a-11eb-35be-3d940d4504cb) above):
"""

# ╔═╡ 5be2e5d2-6a11-11eb-1421-492f5af16f9c
chosenModel

# ╔═╡ 58617378-6aee-11eb-23e8-c13d89b4c57f
md"""
This other cell should return the model function corresponding to the selected binding model (the active radio button in section [Visualizations](#5c5e0392-658a-11eb-35be-3d940d4504cb) above):
"""

# ╔═╡ 88d941e6-658a-11eb-08a2-0f021e5ae3a4
md"""
#### Hill model

This is the [Hill equation](https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)):

$S = S_{min} + (S_{max} - S_{min}) \times \frac{L^h}{{K_D}^h + L^h}$

In which $S$ is the measured signal ($Y$ value) at a given value of ligand concentration $L$ ($X$ value), $S_{min}$ and $S_{max}$ are the minimum and maximum values the observed signal can take, respectively, $K_D$ is the equilibrium dissociation constant and $h$ is the Hill coefficient."""

# ╔═╡ e9fc3d44-6559-11eb-2da7-314e8fc76ee9
@. hill(conc, p) = p[1] + (p[2] - p[1]) * conc^p[4] / (p[3]^p[4] + conc^p[4])

# ╔═╡ 9d1f24cc-6a0f-11eb-3b16-35f89aff5d4a
md"""
#### Hyperbolic model

The hyperbolic equation is a special case of the Hill equation, in which $h = 1$:
"""

# ╔═╡ 78e664d0-6618-11eb-135b-5574bb05ddef
@. hyperbolic(conc, p) = p[1] + (p[2] - p[1]) * conc / (p[3] + conc)

# ╔═╡ b0b17206-6a0f-11eb-2f5e-5fc8fa06cd36
md"""
#### Quadratic model

Unlike the Hill and hyperbolic models, the quadratic model does not make the approximation that the concentration of free ligand at equilibrium is equal to the total ligand concentration:

$S = S_{min} + (S_{max} - S_{min}) \times \frac{(K_{D} + R_{tot} + L_{tot}) - \sqrt{(- K_{D} - R_{tot} - L_{tot})^2 - 4 \times R_{tot} \times L_{tot}}}{2 \times R_{tot}}$

Symbols have the same meaning as in the previous equations, except here $L_{tot}$ is the total concentration of ligand, not the concentration of free ligand at equilibrium. $R_{tot}$ is the total concentration of receptor.

In principle, $R_{tot}$ could be left as a free parameter to be determined by the fitting procedure, but in general it is known accurately enough from the experimental set up, and one should replicate the same experiment with different concentrations of receptor to check its effect on the results. $R_{tot}$ should be set in the experiment to be smaller than $K_D$, ideally, or at least of the same order of magnitude than $K_D$. It might take a couple experiments to obtain an estimate of $K_D$ before one can determine an adequately small concentration of receptor at which to perform a definite experiment.
"""

# ╔═╡ 5694f1da-6636-11eb-0fed-9fee5c48b114
@. quadratic(conc, p) = p[1] + (p[2] - p[1]) * ( (p[3] + R0 + conc) - sqrt((- p[3] - R0 - conc) ^ 2 - 4 * R0 * conc) ) / (2 * R0)

# ╔═╡ 605f06ae-6a11-11eb-0b2f-eb81f6526829
bindingModels = Dict(
	"Hill" => hill,
	"Hyperbolic" => hyperbolic,
	"Quadratic" => quadratic
)

# ╔═╡ 7c03fcbe-6a11-11eb-1b7b-cbad863156a6
function plotOneDataset!(plt, df, fit, filePath, showInitial = false, initialValues = nothing)
	title = split(filePath, "/")[end]
	scatter!(
		plt,
		df.conc,
		df.measurement,
		label = "$title: data"
	)
	if showInitial
		plot!(
			plt,
			df.conc,
			bindingModels[chosenModel](df.concentration, initialValues),
			label = "$title: $chosenModel fit (initial)",
			color = :grey
		)
		plot!(
			plt,
			df.conc,
			bindingModels[chosenModel](df.concentration, fit.param),
			label = "$title: $chosenModel fit (converged)",
			color = :red
		)
	else
		plot!(
			plt,
			df.conc,
			bindingModels[chosenModel](df.concentration, fit.param),
			label = "$title: $chosenModel fit",
			color = :red
		)
	end
end

# ╔═╡ 6668c49a-6a11-11eb-2abf-5feecaee8972
bindingModels[chosenModel]

# ╔═╡ 0f960a8a-6a0f-11eb-04e2-b543192f6354
md"### Parameters and their initial values"

# ╔═╡ 2eb890a4-658c-11eb-1fc4-af645d74109d
md"""
The `findInitialValues()` function takes the measured data and returns an array containing initial values for the model parameters (in this order): $S_{min}$, $S_{max}$, $K_D$ and $h$ (for the Hill model only, so the function needs to know which model was selected).

Initial values for $S_{min}$ and $S_{max}$ are simply taken as the minimal and maximal values found in the data. The initial estimate for $K_D$ is the concentration of the data point that has a signal closest to halfway between $S_{min}$ and $S_{max}$ (if the experiment was properly designed, this is a reasonable estimate and close enough to the true value for the fit to converge). The initial estimate of $h$ is $1.0$, meaning we assume no cooperativity.
"""

# ╔═╡ ca2d2f12-6a1a-11eb-13ca-1f93df2b8e4a
function findInitialValues(df, model)
	halfSignal = minimum(df.mean) + (maximum(df.mean) - minimum(df.mean)) / 2
	if model == "Hill"
		params = [
			minimum(df.mean),
			maximum(df.mean),
			df.concentration[findmin(abs.(halfSignal .- df.mean))[2]],
			1.0
		]
	else
		params = [
			minimum(df.mean),
			maximum(df.mean),
			df.concentration[findmin(abs.(halfSignal .- df.mean))[2]]
		]
	end
end

# ╔═╡ 9be41e32-6af0-11eb-0904-d1cf3c288cab
md"Determine initial values of the selected model's parameters from the currently loaded datasets:"

# ╔═╡ 67b538f6-6a1b-11eb-3004-2d89c2f941e8
initialParams = [ findInitialValues(df, chosenModel) for df in dfs ]

# ╔═╡ 213e8ffa-6a0f-11eb-357e-638146193c5d
md"### Fitting"

# ╔═╡ babcb896-6af0-11eb-194a-15922bc2df83
md"Perform fit of the selected model to the measurements' mean values using initial values for the model parameters determined previously. If the dataset contains replicates, the fit will be weighted by the measurements' standard deviations."

# ╔═╡ 47426056-6af2-11eb-17f8-6d27d35003ca
begin
	fits = Vector{LsqFit.LsqFitResult}(undef, length(allData))
	for (df, initialValues, i) in zip(dfs, initialParams, 1:length(allData))
		if ncol(df) > 5
			# If the dataset has more than 5 columns, it means it has
			# replicate Y values, so we weight the fit by their stddev.
			fits[i] = curve_fit(bindingModels[chosenModel],
						df.concentration,
						df.mean,
						df.std,
						initialValues)
		else
			# If the dataset has only 5 columns, it means it doesn't have
			# replicate Y values, so there are no stddev we can use as weights.
			fits[i] = curve_fit(bindingModels[chosenModel],
						df.concentration,
						df.mean,
						initialValues)
		end
	end
	fits
end

# ╔═╡ 264bf9ec-6af5-11eb-1ffd-79fb3466f596
begin
	dataPlot = initMainPlot()
	for (df, fit, title, initialVals) in zip(dfs, fits, datasetNames, initialParams)
		plotOneDataset!(dataPlot, df, fit, title, showInitialFit, initialVals)
	end
	dataPlot
end

# ╔═╡ a951b5dc-6af7-11eb-2401-5d11a14e3067
begin
	residualPlot = initResidualPlot()
	for (df, fit, title) in zip(dfs, fits, datasetNames)
		plotOneResiduals!(residualPlot, df, fit, title)
	end
	residualPlot
end

# ╔═╡ 7625d41a-dab1-4b10-947c-3667c03f85aa
begin
	residualHistogram = initResidualHistogram()
	for (df, fit, title) in zip(dfs, fits, datasetNames)
		plotOneResidualsHistogram!(residualHistogram, df, fit, title)
	end
	residualHistogram
end

# ╔═╡ 2109f516-6b99-11eb-05a0-99b9ecfd0f9d
PlutoUI.with_terminal() do
	println("Dataset\t\t\t\tSum of squared residuals")
	for (dataset, fit) in zip(datasetNames, fits)
		println(
			split(dataset, "/")[end],
			"\t\t",
			round(sum(fit.resid.^2), digits = 2)
		)
	end
end

# ╔═╡ 799680d0-6af1-11eb-321d-b7758a40f931
md"Degrees of freedom:"

# ╔═╡ 1f0384de-659b-11eb-043e-5b86fcdd36e6
dof.(fits)

# ╔═╡ 54501a10-6b9c-11eb-29de-77afc3772fb7
md"Best fit parameters:"

# ╔═╡ 5ed3ab64-6b9c-11eb-149e-43a1ef12ac7d
coef.(fits)

# ╔═╡ 8643b03c-6af1-11eb-0aa7-67acee28d2c0
md"Standard errors of best-fit parameters:"

# ╔═╡ a74998b4-659c-11eb-354d-09ff62710b87
paramsStdErrors = stderror.(fits)

# ╔═╡ 090347fc-6b8e-11eb-0e17-9d9d45749c0b
PlutoUI.with_terminal() do
	if length(initialParams[1]) == 3
		# No Hill coefficient to report.
		println("Dataset\t\t\t\tKd\t\t\t\tSmin\t\t\tSmax")
		for (dataset, fit, stderr) in zip(datasetNames, fits, paramsStdErrors)
			println(
				split(dataset, "/")[end],
				"\t\t",
				round(fit.param[3], digits = 1),
				" ± ",
				round(stderr[3], digits = 1),
				"\t\t",
				round(fit.param[1], digits = 1),
				" ± ",
				round(stderr[1], digits = 1),
				"\t\t",
				round(fit.param[2], digits = 1),
				" ± ",
				round(stderr[2], digits = 1)
			)
		end
	elseif length(initialParams[1]) == 4
		# There is a Hill coefficient to report.
		println("Dataset\t\t\t\tKd\t\t\t\tSmin\t\t\tSmax\t\t\th")
		for (dataset, fit, stderr) in zip(datasetNames, fits, paramsStdErrors)
			println(
				split(dataset, "/")[end],
				"\t\t",
				round(fit.param[3], digits = 1),
				" ± ",
				round(stderr[3], digits = 1),
				"\t\t",
				round(fit.param[1], digits = 1),
				" ± ",
				round(stderr[1], digits = 1),
				"\t\t",
				round(fit.param[2], digits = 1),
				" ± ",
				round(stderr[2], digits = 1),
				"\t\t",
				round(fit.param[4], digits = 1),
				" ± ",
				round(stderr[4], digits = 1),
			)
		end
	else
		# Other number of values in the parameters array make no sense.
		println("Error.")
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
CSV = "~0.10.4"
Chain = "~0.5.0"
DataFrames = "~1.4.1"
HTTP = "~1.4.1"
LsqFit = "~0.13.0"
Measurements = "~2.8.0"
Plots = "~1.35.3"
PlutoUI = "~0.7.44"
URIs = "~1.4.0"
Unitful = "~1.12.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "49e7e057576ab5c6a76f63f05027c88ded48a354"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e9f7992287edfc27b3cbe0046c544bace004ca5b"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.22"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "84259bb6172806304b9101094a7cc4bc6f56dbc6"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "558078b0b78278683a7445c626ee78c86b9bb000"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "992a23afdb109d0d2f8802a30cf5ae4b1fe7ea68"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "04db820ebcfc1e053bd8cbb8d8bccf0ff3ead3f7"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.76"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "5a2cff9b6b77b33b89f3d97a4d367747adce647e"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.15.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "187198a4ed8ccd7b5d99c41b69c679269ea2b2d4"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.32"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "00a9d4abadc05b9476e937a5557fcce476b9e547"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.69.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "bc9f7725571ddb4ab2c4bc74fa397c1c5ad08943"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.69.1+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "e8c58d5f03b9d9eb9ed7067a2f34c7c371ab130b"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.4.1"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "db619c421554e1e7e07491b85a8f4b96b3f04ca0"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "00f475f85c50584b12268675072663dfed5594b2"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.13.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "6872f9594ff273da6d13c7c1a1545d5a8c7d0c1c"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.6"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "12950d646ce04fb2e89ba5bd890205882c3592d7"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.8.0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "ebe81469e9d7b471d7ddb611d9e147ea16de0add"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.2.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "524d9ff1b2f4473fef59678c06f9f77160a204b1"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.35.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "6e33d318cf8843dade925e35162992145b4eb12f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.44"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "460d9e154365e058c4d886f6f7d6df5ffa1ea80e"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.1.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "3c009334f45dfd546a16a57960a821a1a023d241"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.5.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "612a4d76ad98e9722c8ba387614539155a59e30c"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.0"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "9b1c0c8e9188950e66fc28f40bfe0f8aac311fe0"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.7"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d57a4ed70b6f9ff1da6719f5f2713706d57e0d66"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─412eedac-658a-11eb-2326-93e5bf3d1a2c
# ╟─270cc0cc-660f-11eb-241e-b75746a39cc7
# ╠═a2c02fcf-9402-4938-bc3d-819b45a66afa
# ╠═14baf100-660f-11eb-1380-ebf4a860eed8
# ╟─2ce72e97-0133-4f15-bf1d-7fd04ccf3102
# ╟─77316a92-425a-4902-9828-52a7c4a74f27
# ╟─214acce6-6ae8-11eb-3abf-492e50140317
# ╟─d904fd76-6af1-11eb-2352-837e03072137
# ╟─0e8af3be-7ae7-4ec2-8d7a-670878cd52ee
# ╠═d5b0e2a1-865c-489c-9d0d-c4ae043828fb
# ╟─5c5e0392-658a-11eb-35be-3d940d4504cb
# ╟─3dd72c58-6b2c-11eb-210f-0b547bf38ebe
# ╟─3da83f72-6a11-11eb-1a74-49b66eb39c96
# ╟─d15cba72-6aeb-11eb-2c80-65702b48e859
# ╟─01b59d8a-6637-11eb-0da0-8d3e314e23af
# ╠═264bf9ec-6af5-11eb-1ffd-79fb3466f596
# ╟─4f4000b4-6b2c-11eb-015f-d76a0adda0a0
# ╟─c50cf18c-6b11-11eb-07d3-0b8e332ec5bc
# ╟─5a36fc3f-ce74-42c8-8284-19321e0d687f
# ╟─a951b5dc-6af7-11eb-2401-5d11a14e3067
# ╟─5392d99b-70f9-48cd-90b4-58cba5fc9681
# ╟─7625d41a-dab1-4b10-947c-3667c03f85aa
# ╟─be17b97e-663a-11eb-2158-a381c19ece3f
# ╟─090347fc-6b8e-11eb-0e17-9d9d45749c0b
# ╟─124c4f94-6b99-11eb-2921-d7c2cd00b893
# ╟─2109f516-6b99-11eb-05a0-99b9ecfd0f9d
# ╟─7e7a9dc4-6ae8-11eb-128d-83544f01b78b
# ╟─512e3028-6ae9-11eb-31b4-1bc9fc66b322
# ╠═393b2f5e-6556-11eb-2119-cf7309ee7392
# ╠═abc03f64-6a11-11eb-0319-ed7cea455cb5
# ╟─1884912a-6aeb-11eb-2b4a-d14d4a321dc5
# ╟─c94ef72d-bd12-4434-935e-01e94d5a4588
# ╟─54272681-dfeb-4034-b127-7e68c19fd576
# ╠═731492c6-95c7-449f-8c19-53e22ab438b8
# ╟─93ef0431-643b-4f6f-8c09-beabff59e0c6
# ╠═1c10231d-3bea-4468-a2d8-886c05c6474c
# ╠═730b693d-cca8-46de-8382-c151b5f63352
# ╠═24152a84-523b-4027-9b5a-7e3524b9c659
# ╟─97f73ee4-3db7-43ba-93eb-17025b485f4f
# ╠═8f2d959e-5e61-48a9-bbf0-538bc1d478a8
# ╟─c866d213-4680-4b08-8ecc-1faefc8661a4
# ╠═8e9cd30f-1722-4f2f-a26b-2f558805d4a1
# ╠═13417555-be95-4662-ad11-eca4dece81b5
# ╟─9c7e922e-c82f-41a4-9513-4462a0559c3f
# ╟─4f4b580d-507c-4ad0-b1d5-5967c8ed829e
# ╠═5eb607c7-172b-4a6c-a815-acbc195108f0
# ╟─992c86a2-6b13-11eb-1e00-95bdff2736d0
# ╠═1fe4d112-6a11-11eb-37a6-bf95fbe032b1
# ╟─5cc812af-4f8e-444b-9652-fb063cc6be06
# ╠═7c2af794-a0ed-4827-bc8e-1f15ad205eca
# ╟─a08d1b29-fa5e-4312-977c-39b050de4516
# ╠═0181a0b1-ec6c-4325-8b3a-851a3fe33846
# ╟─247c2416-6c67-11eb-01df-8dac01cdbf8f
# ╠═36ebe112-6c66-11eb-11f0-7fdb946865e4
# ╟─8e105fae-6aec-11eb-1471-83aebb776241
# ╟─97c5020c-6aec-11eb-024b-513b1e603d98
# ╠═caf4a4a2-6aec-11eb-2765-49d67afa47dd
# ╠═7c03fcbe-6a11-11eb-1b7b-cbad863156a6
# ╟─fc194672-6aed-11eb-0a06-2d967ec094b1
# ╠═14db987c-6aee-11eb-06cf-a11987b98f1e
# ╠═7f83b838-6a11-11eb-3652-bdff24f3473e
# ╟─9020fa5d-7408-4161-a52f-df37b3c2e6f5
# ╠═3184d209-1cc9-40ed-a9ec-9f094c5c94b5
# ╠═b38cd229-64e2-4ca4-a78b-8881ec166b09
# ╟─663a4cae-658a-11eb-382f-cf256c08c9d1
# ╟─594e7534-6aeb-11eb-1254-3b92b71877ed
# ╟─32dce844-6aee-11eb-3cf2-3ba420d311d3
# ╠═605f06ae-6a11-11eb-0b2f-eb81f6526829
# ╟─a1f56b0a-6aeb-11eb-0a44-556fad58f368
# ╠═5be2e5d2-6a11-11eb-1421-492f5af16f9c
# ╟─58617378-6aee-11eb-23e8-c13d89b4c57f
# ╠═6668c49a-6a11-11eb-2abf-5feecaee8972
# ╟─88d941e6-658a-11eb-08a2-0f021e5ae3a4
# ╠═e9fc3d44-6559-11eb-2da7-314e8fc76ee9
# ╟─9d1f24cc-6a0f-11eb-3b16-35f89aff5d4a
# ╠═78e664d0-6618-11eb-135b-5574bb05ddef
# ╟─b0b17206-6a0f-11eb-2f5e-5fc8fa06cd36
# ╠═5694f1da-6636-11eb-0fed-9fee5c48b114
# ╟─0f960a8a-6a0f-11eb-04e2-b543192f6354
# ╟─2eb890a4-658c-11eb-1fc4-af645d74109d
# ╠═ca2d2f12-6a1a-11eb-13ca-1f93df2b8e4a
# ╟─9be41e32-6af0-11eb-0904-d1cf3c288cab
# ╠═67b538f6-6a1b-11eb-3004-2d89c2f941e8
# ╟─213e8ffa-6a0f-11eb-357e-638146193c5d
# ╟─babcb896-6af0-11eb-194a-15922bc2df83
# ╠═47426056-6af2-11eb-17f8-6d27d35003ca
# ╟─799680d0-6af1-11eb-321d-b7758a40f931
# ╠═1f0384de-659b-11eb-043e-5b86fcdd36e6
# ╟─54501a10-6b9c-11eb-29de-77afc3772fb7
# ╠═5ed3ab64-6b9c-11eb-149e-43a1ef12ac7d
# ╟─8643b03c-6af1-11eb-0aa7-67acee28d2c0
# ╠═a74998b4-659c-11eb-354d-09ff62710b87
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
