using ArgParse
using BifurcationKit
using DifferentialEquations
using ForwardDiff
using Measures
using Plots
using LaTeXStrings
using ColorSchemes

using PyCall

# Import our python package with PyCall
pushfirst!(pyimport("sys")."path", pwd())  # add current Julia working directory
model_generate_script = pyimport("generate_model")
asp_model = model_generate_script.generate_model()

s = ArgParseSettings()
@add_arg_table s begin
    "--output"
    default = "output/ATHB8_model"
    help = "Path to store output"
    "--plot_var"
    help = "Variable to plot on the y-axis of bifurcation diagrams"
    default = "athb8_tot"
end

parsed_args = parse_args(ARGS, s)



dir = pwd()                       # Or replace with explicit pat
sys = pyimport("sys")
push!(sys."path", dir)
push!(sys."path", "scripts")
println(pyimport("os").getcwd())

py"""
import sys

print(sys.path)
"""

plot_var = parsed_args["plot_var"]
output_dir_name = parsed_args["output"] * "_"  * "$(plot_var)"
output_dir = mkpath(output_dir_name)

# Print out default parameter labels
param_dict = asp_model.get_default_parameters(omit_unused=true)

param_labels = sort(collect(keys(param_dict)))
println(param_labels)

xlabel = "_ARF_1_"

state_labels = asp_model.get_state_variable_names(include_promoter_vars=false)

println("states are $(state_labels)")
# Create dictionary of parameters to change
new_params = Dict()

new_params["k__G2__R__X"] = 1.0e-2
new_params["k__G2__R__ATHB8"]= 0.0

new_params["k__G2__ARF_1__R__X"] = 0
new_params["k__G2__ARF_1__R__ATHB8"]= 1e-2

new_params["k__G2__ATHB8__R__ATHB8"] = 5.0e-2
new_params["k__G2__X__R__ATHB8"] = 0

new_params["k__G2__ATHB8__R__X"] = 0
new_params["k__G2__X__R__X"] = 5.0e-2

new_params["d__R"] = 1e-3

new_params["k__ATHB8"]  = 1.0e-1
new_params["k__X"]  = 1.0e-1

new_params["d__ATHB8"] = 0.1
new_params["d__X"] = 0.1

# Insert new parameters, ensuring that all are present in the model
# Only include parameters that already exist in the model
for key in keys(new_params)
    if key in keys(param_dict)
        param_dict[key] = new_params[key]
    else
        println("WARNING param $(key) is not in param dict")
    end
end

open(joinpath(output_dir, "parameters.txt"), "w") do io
    parameters = asp_model.get_parameters()
    dict = Dict(p.name=>p for p in parameters)

    _param_dict = sort(param_dict, by=(x->hash(dict[x].category)))

    pretty_par_dict = asp_model.generate_pretty_print_parameter_dict()

    pretty_state_dict = asp_model.generate_pretty_print_state_dict()

    pretty_state_dict = Dict(Regex("\\b$(escape_string(k))\\b")=>"\$ $(v) \$"
                             for (k, v) in pretty_state_dict)

    descs = Dict(k => replace(v.description, collect(pretty_state_dict)...)
                 for (k, v) in dict)

    s = join(["\$ $(pretty_par_dict[k])\$ & $v & $(descs[k]) \\\\ "
              for (k, v) in _param_dict],
             "\n")

    s = replace(s, collect(pretty_par_dict)...)
    s = replace(s, "ARF_"=>"A_", "iaa_"=>"I_", "__"=>"_")

    write(io, s)
end

# Generate derivative function using the Python package
f_deriv_func = asp_model.output_julia()

println(f_deriv_func)

f_deriv = eval(Meta.parse(f_deriv_func))

## Set up a dictionary of max values to use in plots
# Set all max_p values to 1.0 by default
max_p_dict = Dict()
for (k, v) in param_dict
    max_p_dict[k] = 1.0
end

max_p_dict["_ARF_1_"] = 1000.0
max_p_dict["k__arf__all"] = 1.0
max_p_dict["d__arf__all"] = 1.0

max_p_dict["d__iaa__all"] = 10.0
max_p_dict["k__D__ARF_1__ARF_1"] = 1e1
max_p_dict["k__D__ARF_1__iaa_1"] = 1e1
max_p_dict["k__D__iaa_1__iaa_1"] = 1e1
max_p_dict["k__RNA__IAA"] = 1e1
max_p_dict["d__R"] = 1e1
max_p_dict["d__D__ARF_1__iaa_1"] = 1e1
max_p_dict["auxin"] = 10.0

par_pp = Float64.(collect(values(sort(param_dict))))
push!(par_pp, 0.0)
println(par_pp)

# Solve system with high and low auxin
tspan = (0.0, 1e4)
z0 = Float64.(collect(values(sort(asp_model.get_default_initial_conditions()))))

println(z0)
println(param_dict)
println(par_pp)

print("Eval deriv func: ")
println(f_deriv(copy(z0), copy(z0), copy(par_pp)))
prob = ODEProblem(f_deriv, z0, tspan, par_pp)
sol = DifferentialEquations.solve(prob)
z0 = sol[:, size(sol, 2)]

legend_labels = reshape(state_labels, 1, length(state_labels))
s = plot(sol, label=legend_labels, legend=true)
savefig(s, "$(output_dir)/ode_sol1.pdf")

print("SS with Auxin conc. = 10.0, ")
println(Dict(zip(state_labels, z0)))

param_index = findfirst(==(xlabel), param_labels)

function athb8_total_func(x)
    include_indices = []
    for (i, key) in enumerate(state_labels)
        if occursin("ATHB8", key) && !occursin("R__", key) && !occursin("G__", key)
            for j in 1:length(collect(eachmatch(Regex("ATHB8"), key)))
                push!(include_indices, i)
            end
        end
    end
    return sum(x[include_indices])
end

const default_ic_dict = asp_model.get_default_initial_conditions()

plot_var = parsed_args["plot_var"]
if plot_var == "iaa_rna"
    plot_func = rna_iaa_tot
elseif plot_var == "arf1_rna"
    plot_func = (x) -> arf_rna_func(x, 1)
elseif plot_var == "athb8_tot"
    plot_func = athb8_total_func
else
    throw("plot var name not recognised")
end

plot_var_s = Symbol(plot_var)

recordFromSolution(x, p; k...) = (athb8_tot=athb8_total_func(x), )

z0 = copy(Float64.(collect(values(sort(asp_model.get_default_initial_conditions())))))

p_max = max_p_dict[xlabel]
p_init = p_max * 0.9
# par_pp[param_index] = p_init
tspan = (0.0, 1e4)

prob = ODEProblem(f_deriv, z0, tspan, par_pp)
sol = DifferentialEquations.solve(prob)
z0 = sol[:, size(sol, 2)]

print("Eval deriv func: ")
println(f_deriv(copy(z0), copy(z0), copy(par_pp)))

println(Dict(zip(state_labels, z0)))

par_pp[end] = 1.0
prob = ODEProblem(f_deriv, z0, tspan, par_pp)
sol = DifferentialEquations.solve(prob)

legend_labels = reshape(state_labels, 1, length(state_labels))
s = plot(sol, label=legend_labels, legend=true)
savefig(s, "$(output_dir)/ode_sol.pdf")

z0 = sol[:, size(sol, 2)]
println(Dict(zip(state_labels, z0)))
println(Dict(zip(state_labels, z0)))

dsmax = p_max * 1e-2
dsmin = p_init * 1e-9
ds = (dsmax + dsmin) / 2.0

# bifurcation problem
prob1 = BifurcationProblem(f_deriv, z0, par_pp,
                           (@optic _[param_index]),
                           record_from_solution=recordFromSolution,
                           )

opts_br = ContinuationPar(p_max = p_max, p_min=0.01,
                          dsmax=dsmax, dsmin=dsmin, ds = ds,
                          nev=10, max_steps=100000,
                          detect_bifurcation=3)

diagram1 = bifurcationdiagram(prob1, PALC(),
                              3,
                              opts_br,
                              bothside=true
                              )

par_pp[param_index] = 0.01
# bifurcation problem
prob1 = BifurcationProblem(f_deriv, z0, par_pp,
                           (@optic _[param_index]),
                           record_from_solution=recordFromSolution,
                           )

opts_br = ContinuationPar(p_max = p_max, p_min=1e-8,
                          dsmax=dsmax, dsmin=dsmin, ds = ds,
                          nev=10, max_steps=100000,
                          detect_bifurcation=3)

diagram2 = bifurcationdiagram(prob1, PALC(),
                              3,
                              opts_br,
                              bothside=true
                              )


scene = plot(diagram1; code=(), legend=true,
             vars=(:param, plot_var_s),
             xlabel=xlabel)
plot!(diagram2)
savefig(scene, "$(output_dir)/bifurcation_diagram.pdf")
