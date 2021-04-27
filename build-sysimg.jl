using PackageCompiler
using Pkg

Pkg.activate(".")

create_sysimage(
    [
        :Chain,
        :CSV,
        :DataFrames,
        :HTTP,
        :LsqFit,
        :Measurements,
        :Plots,
        :Pluto,
        :PlutoUI,
        :Statistics,
        :Unitful,
        :UnitfulRecipes,
        :URIs
    ],
    sysimage_path = "sysimg.dylib"
)
