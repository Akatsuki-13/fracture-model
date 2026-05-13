# Pilot Test for Bone Repair Model 
#= This is to test the condition of the environment over time by testing:
    - overall tissue energy over time 
    - M1:M2 macrophage ratio over time
    - ROS count over time
    - osteoblast:osteoclast ratio over time
    - the number of each agent type over time 
This gives us an idea of how much time passes before the results plateau, and 
shows how results may (or may not, i.e., the null hypothesis) vary over time.
It may also provide insight into whether or not there are 'seasonal' pattenrs
over the course of the model run times. 
=#

include("C:/Users/Carrie/OneDrive - Nottingham Trent University/UG Research Project/Bone Repair Model.jl")
using Statistics
using Plots
using DataFrames
using HypothesisTests
using CSV

# Set save path
cd("C:/Users/Carrie/OneDrive - Nottingham Trent University/UG Research Project")

# Set the length of time the model runs here: 
t_end = 100.0

# Mitochondrial health per MDM
mito_health_MDM(a) = a isa MDM && !isempty(a.mitochondria) ?
    mean(h for (_, h) in a.mitochondria) : NaN

mito_count_MDM(a) = a isa MDM ? length(a.mitochondria) : 0 

# Agent-level data collectors 
adata = [
    # Molecular agent counts
    (a -> a isa Glucose,    count),
    (a -> a isa Oxygen,     count),
    (a -> a isa M1_cytokines, count),
    (a -> a isa M2_cytokines, count),
    (a -> a isa ROS,          count),
    (a -> a isa DAMPs,        count),
    (a -> a isa CellDebris,   count),

    # Cell counts 
    (a -> a isa MDM,          count),
    (a -> a isa Osteomac,     count),
    (a -> a isa Osteoblast,   count),
    (a -> a isa Osteoclast,   count),
    (a -> a isa Osteocyte,    count),

    # Overall tissue energy
    (:energy, mean, a -> a isa MDM || a isa Osteomac ||
                        a isa Osteoblast || a isa Osteoclast ||
                        a isa Osteocyte),
    
    # Energy by cell type 
    (:energy, mean, a -> a isa MDM),
    (:energy, mean, a -> a isa Osteomac),
    (:energy, mean, a -> a isa Osteoblast),
    (:energy, mean, a -> a isa Osteoclast),
    (:energy, mean, a -> a isa Osteocyte),

    # Internal ROS (oxidative stress)
    (:ROS, mean, a -> a isa MDM || a isa Osteomac ||
                      a isa Osteoblast || a isa Osteoclast ||
                      a isa Osteocyte),
    (mito_health_MDM, mean, a -> a isa MDM),
    (mito_count_MDM, mean, a -> a isa MDM),

    # Counting M1 macrophages
    (a -> (a isa MDM || a isa Osteomac) && a.plasticity < 0.5, count),

    # Counting M2 macrophages 
    (a -> (a isa MDM || a isa Osteomac) && a.plasticity >= 0.5, count)
]

# Initialize and run each condition 
low_low = initialize_bone_model(n_glucose = 25, n_oxygen = 150)
low_mid = initialize_bone_model(n_glucose = 25, n_oxygen = 650)
low_high = initialize_bone_model(n_glucose = 25, n_oxygen = 2600)
mid_low = initialize_bone_model(n_glucose = 100, n_oxygen = 150)
mid_mid = initialize_bone_model(n_glucose = 100, n_oxygen = 650)
mid_high = initialize_bone_model(n_glucose = 100, n_oxygen = 2600)
high_low = initialize_bone_model(n_glucose = 400, n_oxygen = 150)
high_mid = initialize_bone_model(n_glucose = 400, n_oxygen = 650)
high_high = initialize_bone_model(n_glucose = 400, n_oxygen = 2600)

println("Running hypoxia with hypoglycemic conditions...")
adf_hypo, _ = run!(low_low, t_end; adata)

println("Running normal oxygen with hypoglycemic conditions...")
adf_low_mid, _ = run!(low_mid, t_end; adata)

println("Running oxygen-rich and hypoglycemic conditions...")
adf_low_high, _ = run!(low_high, t_end; adata)

println("Running hypoxia with normal glucose conditions...")
adf_mid_low, _ = run!(mid_low, t_end; adata)

println("Running normal oxygen with normal glucose conditions...")
adf_mid, _ = run!(mid_mid, t_end; adata)

println("Running oxygen-rich environment with normal glucose conditions...")
adf_mid_high, _ = run!(mid_high, t_end; adata)

println("Running hypoxic and hyperglycemic conditions...")
adf_high_low, _ = run!(high_low, t_end; adata)

println("Running normal oxygen with hyperglycemic conditions...")
adf_high_mid, _ = run!(high_mid, t_end; adata)

println("Running oxygen-rich and hyperglycemic conditions...")
adf_hyper, _ = run!(high_high, t_end; adata)

# Tag each dataframe with condition
insertcols!(adf_hypo,     :condition => "depleted")
insertcols!(adf_low_mid,  :condition => "hypoglycemic/normal oxygen")
insertcols!(adf_low_high, :condition => "hypoglycemic/high oxygen")
insertcols!(adf_mid_low,  :condition => "normal glucose/hypoxia")
insertcols!(adf_mid,      :condition => "normal")
insertcols!(adf_mid_high, :condition => "normal glucose/high oxygen")
insertcols!(adf_high_low, :condition => "hyperglycemic/hypoxia")
insertcols!(adf_high_mid, :condition => "hyperglycemic/normal oxygen")
insertcols!(adf_hyper,    :condition => "hyperglycemic/high oxygen")

# Combine into one dataframe
all_data = vcat(adf_hypo, adf_low_mid, adf_low_high, adf_mid_low, adf_mid, adf_mid_high,
adf_high_low, adf_high_mid, adf_hyper)

# Round time to nearest integrate_average
all_data[!, :time_rounded] = round.(all_data.time, digits = 0)

# Rename columns here: 
let cols = names(all_data)
    # The anonymous lambda columns appear in adata order.
    # Columns 3 onward (after :time, :id) match adata entries in order.
    expected = [
        "count_Glucose", "count_Oxygen", "count_M1_cyto", "count_M2_cyto",
        "count_ROS", "count_DAMPs", "count_CellDebris",
        "count_MDM", "count_Osteomac", "count_Osteoblast",
        "count_Osteoclast", "count_Osteocyte",
        "mean_energy_all", "mean_energy_MDM", "mean_energy_Osteomac",
        "mean_energy_Osteoblast", "mean_energy_Osteoclast", "mean_energy_Osteocyte",
        "mean_ROS_cells", "mean_mito_health_MDM", "mean_mito_count_MDM",
        "count_M1_macro", "count_M2_macro"
    ]
    adata_cols = filter(c -> c != "time" && c != "id" && c != "condition", cols)
    rename!(all_data, Dict(zip(adata_cols, expected)))
end

get_condition(df, cond) = df[df.condition .== cond, :]

conditions = [
    "depleted",
    "hypoglycemic/normal oxygen",
    "hypoglycemic/high oxygen",
    "normal glucose/hypoxia",
    "normal",
    "normal glucose/high oxygen",
    "hyperglycemic/hypoxia",
    "hyperglycemic/normal oxygen",
    "hyperglycemic/high oxygen"
]

# Plot line colors for each condition 
condition_colors = Dict(
    "depleted"                    => :black,
    "hypoglycemic/normal oxygen"  => :steelblue1,
    "hypoglycemic/high oxygen"    => :royalblue4,
    "normal glucose/hypoxia"      => :darkorange,
    "normal"                      => :green4,
    "normal glucose/high oxygen"  => :purple2,
    "hyperglycemic/hypoxia"       => :violetred4,
    "hyperglycemic/normal oxygen" => :orangered4,
    "hyperglycemic/high oxygen"   => :mediumvioletred
)

# Counting the number or each type of cell over the course of the model for each environmental condition
cell_cols = [
    (:count_MDM,    "MDM Count Over Time"),
    (:count_Osteomac,   "Osteomac Count Over Time"),
    (:count_Osteoblast, "Osteoblast Count Over Time"),
    (:count_Osteoclast, "Osteoclast Count Over Time"),
    (:count_Osteocyte, "Osteocyte Count Over Time"),
]
for col in [:count_MDM, :count_Osteomac, :count_Osteoblast, :count_Osteoclast, :count_Osteocyte]
    p1 = Plots.plot(title = "Cells Over Time", xlabel = "Time", ylabel = "Count", legend = :outertopright, size = (800, 400))
    for cond in conditions 
        sub = get_condition(all_data, cond)
        Plots.plot!(p1, sub.time_rounded, sub[:, col], 
        label = cond, color = condition_colors[cond], linewidth = 3)
    end
    Plots.savefig(p1, "Cells Over Time.png")
end


# Tissue energy over time 
p2 = Plots.plot(title = "Mean Tissue Energy Over Time", xlabel = "Time", ylabel = "Mean Energy", legend = :outertopright, size = (800, 400))
for cond in conditions 
    sub = get_condition(all_data, cond)
    Plots.plot!(p2, sub.time_rounded, sub[:, :mean_energy_all], 
    label = cond, color = condition_colors[cond], linewidth = 3)
end
Plots.savefig(p2, "Mean Energy Over Time.png")


# Mean ROS over time 
p3 = Plots.plot(title = "Mean Intracellular ROS Over Time", xlabel = "Time", ylabel = "Mean ROS", legend = :outertopright, size = (800, 400))
for cond in conditions 
    sub = get_condition(all_data, cond)
    Plots.plot!(p3, sub.time_rounded, sub[:, :mean_ROS_cells], 
    label = cond, color = condition_colors[cond], linewidth = 3)
end
Plots.savefig(p3, "Mean Intracellular ROS Over Time.png")


# Ratio of M1:M2 macrophages over time 
p4 = Plots.plot(title = "M1:M2 Ratio", xlabel = "Time", ylabel = "Ratio", legend = :outertopright, size = (800, 400))
for cond in conditions
    sub = get_condition(all_data, cond)
    ratio = (sub[:, :count_M1_macro] .+ 1) ./ (sub[:, :count_M2_macro] .+ 1)
    Plots.plot!(p4, sub.time_rounded, ratio, label = cond, color = condition_colors[cond], linewidth = 3)
end
Plots.savefig(p4, "M1 to M2 Macrophage Ratio Over Time.png")

# Ratio of osteoblasts:osteoclasts over time 
p5 = Plots.plot(title = "Osteoblast:Osteoclast Ratio", xlabel = "Time", ylabel = "Ratio", legend = :outertopright, size = (800, 400))
for cond in conditions
    sub = get_condition(all_data, cond)
    ratio = (sub[:, :count_Osteoblast] .+ 1) ./ (sub[:, :count_Osteoclast] .+ 1)
    Plots.plot!(p5, sub.time_rounded, ratio, label = cond, color = condition_colors[cond], linewidth = 3)
end
Plots.savefig(p5, "Osteoblast to Osteoclast Ratio Over Time.png")

println("All plots saved as .png files - check your folder")