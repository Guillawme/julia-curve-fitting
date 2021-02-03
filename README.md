# Non-linear curve fitting in Julia

This is an example [Pluto](https://github.com/fonsp/Pluto.jl) notebook to
perform non-linear curve fitting in [Julia](https://julialang.org).

## Datasets

Example datasets are from the following publication:

Gaullier G, Roberts G, Muthurajan UM, Bowerman S, Rudolph J, Mahadevan J, Jha A,
Rae PS & Luger K (2020) Bridging of nucleosome-proximal DNA double-strand breaks
by PARP2 enhances its interaction with HPF1. *PLOS ONE* **15**: e0240932
<https://doi.org/10.1371/journal.pone.0240932>

The CSV files in the repository's `datasets` directory have been reformatted to
be easy to read into a `DataFrame`. Original files can be downloaded at
<https://doi.org/10.5281/zenodo.3519435>

## How to use

To use this notebook, clone the repository, navigate into it, start Julia and
run the following commands (type `]` to get to the `pkg>` prompt, `backspace` at
an empty `pkg>` prompt to get back to the `julia>` prompt):

``` julialang
pkg> activate .
pkg> instantiate # this command is only needed the first time
julia> using Pluto
julia> Pluto.run()
```

