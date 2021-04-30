### A Pluto.jl notebook ###
# v0.14.4

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

# ╔═╡ b6ff5052-04ec-11eb-043b-458c42b19edb
begin
	using Plots
	using PlutoUI
end

# ╔═╡ 87c1dc3b-d361-4f46-a68a-c28d7c7de0c3
md"""
# Equilibrium binding simulations

This notebook can either:

1. simulate the outcome of an equilibrium binding experiment given a hypothetical $K_D$ and desired experimental parameters (section [Simulate a binding curve](#d50721c7-34e4-411d-ae6b-dc6196b0f5a9)) or
2. calculate the fractional saturation of receptor given a $K_D$ and total concentrations of receptor and ligand (section [Predict saturation of receptor](#3954adfa-73cc-495a-93c1-bcde87a60ce8)).

The following reference explains very well the theory of equilibrium binding experiments, as well as many important practical considerations:

> Jarmoskaite I, AlSadhan I, Vaidyanathan PP & Herschlag D (2020) How to measure and evaluate binding affinities. *eLife* **9**: e57264 <https://doi.org/10.7554/eLife.57264>
"""

# ╔═╡ d50721c7-34e4-411d-ae6b-dc6196b0f5a9
md"""
## Simulate a binding curve

Enter values of $K_D$ and experimental parameters below to simulate the resulting data points (sampling the binding curve) as they would come out if you performed the experiment. The vertical line on the graph shows the value of $K_D$. The horizontal lines show zero and full saturation of the receptor.
"""

# ╔═╡ c4d66e7a-38e4-497d-b680-8c50c51fe146
md"""
Value of $K_D$: $(@bind kd PlutoUI.Slider(0.1:0.1:100.0, default = 10.0, show_value = true))

Concentration of ligand at highest titration point $L_{max}$: $(@bind lmax PlutoUI.Slider(kd:1.0:100000.0, default = 1000.0, show_value = true))

Dilution factor between titration points: $(@bind dilution PlutoUI.NumberField(1.1:0.1:10.0, default = 2.0))

Number of titration points: $(@bind points PlutoUI.NumberField(1:1:24, default = 15))

Concentration of receptor $R_{tot}$: $(@bind Rtot PlutoUI.Slider(0.001:1.0:50.0, default = 1.0, show_value = true))
"""

# ╔═╡ d49dee36-c12d-4206-9699-a181570caa80
md"""
Ratio $\frac{L_{max}}{K_D} =$ $(lmax / kd)

Ratio $\frac{K_D}{R_{tot}} =$ $(kd / Rtot)
"""

# ╔═╡ 3954adfa-73cc-495a-93c1-bcde87a60ce8
md"""
## Predict saturation of receptor

Value of $K_D$: $(@bind kd2 PlutoUI.NumberField(1.0:10.0:1000.0, default = 10.0))

Concentration of ligand $L_{tot}$: $(@bind ltot2 PlutoUI.NumberField(1.0:100.0:100000.0, default = 1000.0))

Concentration of receptor $R_{tot}$: $(@bind rtot2 PlutoUI.NumberField(0.001:0.01:100000.0, default = 1000.0))

Percent saturation of receptor $B =$
"""

# ╔═╡ ae1d6dac-92d7-4962-a641-63965f73a8f2
md"""
## Code

The quadratic model gives $B$, the fraction of bound receptor (or fractional saturation of receptor) for any values of receptor concentration $R_{tot}$ and ligand concentration $L_{tot}$, given a binding affinity determined by the value of $K_D$:

$B = \frac{(K_{D} + R_{tot} + L_{tot}) - \sqrt{(- K_{D} - R_{tot} - L_{tot})^2 - 4 \times R_{tot} \times L_{tot}}}{2 \times R_{tot}}$

The function below defines this binding model in code:
"""

# ╔═╡ 1eed6126-331b-4218-acc7-8b7c132ef355
function quadratic(Ltot)
	( (kd + Rtot + Ltot) - sqrt((- kd - Rtot - Ltot) ^ 2 - 4 * Rtot * Ltot) ) / (2 * Rtot)
end

# ╔═╡ a5ac735b-f748-4d37-9cf1-7e64e1ea4c53
function quadratic(kd, Rtot, Ltot)
	( (kd + Rtot + Ltot) - sqrt((kd + Rtot + Ltot) ^ 2 - 4 * Rtot * Ltot) ) / (2 * Rtot)
end

# ╔═╡ 4d816572-03a3-481c-a29c-522f0535a679
round(100 * quadratic(kd2, ltot2, rtot2), digits = 2)

# ╔═╡ 3794717b-1f61-4622-9a55-ad86ab947831
md"""
The code below calculates concentrations of the titration series ($X$ values for the plot) given the maximal ligand concentration $L_{max}$, the dilution factor and the number of titration points entered above:
"""

# ╔═╡ 7c43edeb-a5e6-4fb4-8494-e2072d0d56fd
x = [ lmax / dilution^i for i in 0:points-1 ]

# ╔═╡ 06f341e8-6d40-4eb7-b7e9-58c5d45f1f85
md"""
The code below uses the quadratic model to calculate fractional saturation of the receptor ($Y$ values for the plot) at each titration point for the values of $K_D$ and $R_{tot}$ entered above:
"""

# ╔═╡ 82af41ba-d9b9-4ce2-8736-6fe32c8e024e
y = quadratic.(x)

# ╔═╡ efc70ae1-abcc-472d-8df8-c55cbbd6529b
begin
	scatter(
		x,
		y,
		xscale = :log10,
		legend = :none,
		xlabel = "Ligand concentration",
		ylabel = "Predicted fractional saturation of receptor"
	)
	vline!([kd], color = :red)
	hline!([1.0], color = :black)
	hline!([0.0], color = :black)
end

# ╔═╡ b0340ca3-02d5-49a9-841b-7827e4aee2bb
PlutoUI.TableOfContents()

# ╔═╡ 9c5021f3-67ab-464d-8f22-25e269c36083
md"""
## Required packages
"""

# ╔═╡ Cell order:
# ╟─87c1dc3b-d361-4f46-a68a-c28d7c7de0c3
# ╟─d50721c7-34e4-411d-ae6b-dc6196b0f5a9
# ╟─c4d66e7a-38e4-497d-b680-8c50c51fe146
# ╟─d49dee36-c12d-4206-9699-a181570caa80
# ╟─efc70ae1-abcc-472d-8df8-c55cbbd6529b
# ╟─3954adfa-73cc-495a-93c1-bcde87a60ce8
# ╠═4d816572-03a3-481c-a29c-522f0535a679
# ╟─ae1d6dac-92d7-4962-a641-63965f73a8f2
# ╠═1eed6126-331b-4218-acc7-8b7c132ef355
# ╠═a5ac735b-f748-4d37-9cf1-7e64e1ea4c53
# ╟─3794717b-1f61-4622-9a55-ad86ab947831
# ╠═7c43edeb-a5e6-4fb4-8494-e2072d0d56fd
# ╟─06f341e8-6d40-4eb7-b7e9-58c5d45f1f85
# ╠═82af41ba-d9b9-4ce2-8736-6fe32c8e024e
# ╟─9c5021f3-67ab-464d-8f22-25e269c36083
# ╠═b6ff5052-04ec-11eb-043b-458c42b19edb
# ╠═b0340ca3-02d5-49a9-841b-7827e4aee2bb
