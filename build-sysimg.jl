using PackageCompiler
using Pkg

Pkg.activate(".")

create_sysimage(
    [
        :Chain,
        :CSV,
        :DataFrames,
	:Distributions,
        :HTTP,
        :LsqFit,
        :Measurements,
        :Plots,
        :Pluto,
        :PlutoUI,
	:Random,
        :Statistics,
        :Unitful,
        :UnitfulRecipes,
        :URIs
    ],
    sysimage_path = "sysimg.dylib"
)
