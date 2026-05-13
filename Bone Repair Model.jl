# Setup 

using Agents
using Random 

# Define cell types 

# For 'mitochondria' field:
# Int is the number assigned to the mitochondrion
# Float64 is the health value of the mitochondrion

@agent struct MDM(GridAgent{2})
    mitochondria::Vector{Tuple{Int, Float64}} 
    energy::Float64                     
    plasticity::Float64 
    mitophagy_propensity::Float64 
    mito_fusion_propensity::Float64 
    mito_fission_propensity::Float64 
    oxphos_propensity::Float64 
    glycolysis_propensity::Float64
    ROS::Int 
end

@agent struct Osteomac(GridAgent{2})
    mitochondria::Vector{Tuple{Int, Float64}}
    energy::Float64 
    plasticity::Float64 
    mitophagy_propensity::Float64 
    mito_fusion_propensity::Float64 
    mito_fission_propensity::Float64 
    oxphos_propensity::Float64 
    glycolysis_propensity::Float64
    ROS::Int 
end

@agent struct Osteoblast(GridAgent{2})
    mitochondria::Vector{Tuple{Int, Float64}} 
    energy::Float64 
    mitophagy_propensity::Float64 = 1.0
    mito_fusion_propensity::Float64 = 1.0
    mito_fission_propensity::Float64 = 1.0
    ROS::Int = 0 
end

@agent struct Osteoclast(GridAgent{2})
    mitochondria::Vector{Tuple{Int, Float64}}
    energy::Float64
    mitophagy_propensity::Float64 = 1.0
    mito_fusion_propensity::Float64 = 1.0
    mito_fission_propensity::Float64 = 1.0
    ROS::Int = 0 
end

@agent struct Osteocyte(GridAgent{2})
    mitochondria::Vector{Tuple{Int, Float64}}
    energy::Float64
    mitophagy_propensity::Float64 = 1.0
    mito_fusion_propensity::Float64 = 1.0
    mito_fission_propensity::Float64 = 1.0
    ROS::Int = 0
end

# Create molecular agents 

@agent struct Glucose(GridAgent{2})
end 

@agent struct Oxygen(GridAgent{2})
end

@agent struct M1_cytokines(GridAgent{2})
end

@agent struct M2_cytokines(GridAgent{2})
end

@agent struct DAMPs(GridAgent{2})
end

@agent struct CellDebris(GridAgent{2})
end

@agent struct ROS(GridAgent{2})
end

# Create fracture unit agent 
# During initialization, FUs will generate adjacent to one another, creating 
# a line across the model, representing a fracture line 

@agent struct FractureUnit(GridAgent{2})
end

agents = Union{MDM, Osteomac, Osteoblast, Osteoclast, Osteocyte, Glucose, Oxygen,
M1_cytokines, M2_cytokines, DAMPs, CellDebris, ROS, FractureUnit}

# Functions:
# Create a helper functions to compute oxphos_propensity + glycolysis_propensity

function update_oxphos!(agent)
    agent.oxphos_propensity = clamp(0.1*(100)^agent.plasticity, 0, 10)
end

function update_glycolysis!(agent)
    agent.glycolysis_propensity = clamp(1.7783*(0.003162)^agent.plasticity, 0, 10)
end

# Create an attack function for extracellular ROS

function attack!(agent, model)
    # Collect nearby agents
    contenders = collect(nearby_agents(agent, model, 1))

    # Filter for only cells
    nearby_cells = filter(n -> n isa MDM || n isa Osteomac
                          || n isa Osteoclast || n isa Osteoblast, contenders)
    isempty(nearby_cells) && return

    # Pick one cell to attack
    target = rand(nearby_cells)

    # Move plastiscity towards M1 side of spectrum
    if target isa MDM || target isa Osteomac
        target.plasticity -= 0.05
        update_oxphos!(target)
        update_glycolysis!(target)
    else 
        return
    end

    # Apply damage to random mitochondria
    isempty(target.mitochondria) && return
    mito_target = rand(target.mitochondria)

    target_health = last(mito_target) - 0.05

    # If target mitochondria's health = 0, remove instance from the tuple

    if target_health <= 0
        idx = findfirst(==(mito_target), target.mitochondria)
        deleteat!(target.mitochondria, idx)
    else
        # update mito health 
        idx = findfirst(==(mito_target), target.mitochondria)
        target.mitochondria[idx] = (mito_target[1], target_health)
    end

    # Reorganize tuples, where x-values are ordered 1:n_tuples
    for (i, (_, health)) in enumerate(target.mitochondria)
        target.mitochondria[i] = (i, health)
    end
end
       

# Create an internal attack function for intracellular ROS

function internal_attack!(agent, model)
    # If agent.ROS is more than 5, proceed
    agent.ROS <= 5 && return

    # Attack the health of random mitochondria 
   isempty(agent.mitochondria) && return
   mito_target = rand(agent.mitochondria)

    target_health = mito_target[2] - 0.05

    # -1 ROS after attack  
    agent.ROS -= 1

    # If target mitochondria's health = 0, remove instance from the tuple

    if target_health <= 0
        idx = findfirst(==(mito_target), agent.mitochondria)
        deleteat!(agent.mitochondria, idx)
    else 
        idx = findfirst(==(mito_target), agent.mitochondria)
        agent.mitochondria[idx] = (mito_target[1], target_health)
    end

    # Reorganize tuples, where x-values are ordered 1:n_tuples
    for (i, (_, health)) in enumerate(agent.mitochondria)
         agent.mitochondria[i] = (i, health)
    end

    # Release DAMPs 
    nbrs = collect(nearby_positions(agent.pos, model, 1))
    isempty(nbrs) && return 

    pos = rand(abmrng(model), nbrs)
    add_agent!(pos, DAMPs, model)
end



# Create a function where molecular agents that touch the fracture line 
# "disappear" into the system - represents hypoxia and nutrients deficiency 

function abyss!(agent, model)
    # Check if there is any FractureUnit at this position 
    fus = filter(a -> a isa FractureUnit,
                 collect(agents_in_position(agent.pos, model)))

    # If agent shares pos with fracture unit, delete agent 
    if !isempty(fus) 
        remove_agent!(agent, model)
    end
end


# Create an ROS removal function 

function removeROS!(agent, model)
    # Remove one ROS using one energy point 
    agent.ROS <= 0 && return 
    agent.energy <= 1.0 && return
    agent.ROS -= 1
    agent.energy -= 1
end


# Create active and transient move functions

function transient_move!(agent, model) 
    dx = rand(abmrng(model), -1:1)
    dy = rand(abmrng(model), -1:1)
    gs = spacesize(model)
    newpos = (mod1(agent.pos[1] + dx, gs[1]),
              mod1(agent.pos[2] + dy, gs[2]))
    move_agent!(agent, newpos, model)
end

function active_move!(agent, model)
    # Inhibit movement if energy < 3
    if hasfield(typeof(agent), :energy) && agent.energy < 3
        return 
    end
    dx = rand(abmrng(model), -1:1)
    dy = rand(abmrng(model), -1:1)
    gs = spacesize(model)
    newpos = (mod1(agent.pos[1] + dx, gs[1]),
              mod1(agent.pos[2] + dy, gs[2]))
    move_agent!(agent, newpos, model)
    agent.energy -= 0.1
end

# Create 'idle!' function

function idle!(agent, model)
    return 
end


# Create a damage function for cells affected by fracture line 

function damage!(fu::FractureUnit, model)
    occupants = collect(agents_in_position(fu.pos, model))
    
    # Filter for only cellular agents 
    cells = filter(a ->
    a isa MDM ||
    a isa Osteomac ||
    a isa Osteoblast ||
    a isa Osteoclast ||
    a isa Osteocyte, 
    occupants
    )

    isempty(cells) && return 

    # Reduce health by 3
    for cell in cells
        cell.energy -= 3
    end

    nbrs = collect(nearby_positions(fu.pos, model, 1))
    isempty(nbrs) && return

    pos = rand(abmrng(model), nbrs)
    add_agent!(pos, DAMPs, model)
end


# Create 'M1_polarization!' function

function M1_polarization!(agent, model)
    # Collects nearby agents
    contenders = collect(nearby_agents(agent, model, 1))

    # Filter for macrophages only 
    macrophages = filter(m -> m isa MDM || m isa Osteomac, contenders)
    isempty(macrophages) && return

    # Pick one macrophage to polarize
    target = rand(abmrng(model), macrophages)

    # Move macrophage plasticity towards M1 spectrum 
    target.plasticity -= 0.05
    update_oxphos!(target)
    update_glycolysis!(target)

    # Remove cytokine from the environment 
    remove_agent!(agent, model)
end


# Create an 'M2_polarization!' function 

function M2_polarization!(agent, model)
    # Collect nearby agents
    contenders = collect(nearby_agents(agent, model, 1))

    # Filter for macrophages only
    macrophages = filter(m -> m isa MDM || m isa Osteomac, contenders)
    isempty(macrophages) && return 

    # Pick one macrophage to polarize 
    target = rand(abmrng(model), macrophages)

    # Move macrophage plasticity towards M2 spectrum 
    target.plasticity += 0.05
    update_oxphos!(target)
    update_glycolysis!(target)

    # Remove cytokine from environment
    remove_agent!(agent, model)
end


# Create a function to release pro-inflammatory cytokines 

function releaseM1cyto!(agent, model)
    # Collect nearby agents and filter for M1 cytokines 
    contenders = collect(nearby_agents(agent, model, 1))
    M1cytokines = filter(c -> c isa M1_cytokines, contenders)
    isempty(M1cytokines) && return 

    # Pick one cytokine to interact with 
    target = rand(abmrng(model), M1cytokines)

    # Remove cytokine from model 
    remove_agent!(target, model)
    
    # Relase cytokine in response 
    nbrs = collect(nearby_positions(agent.pos, model, 1))
    isempty(nbrs) && return 

    pos = rand(abmrng(model), nbrs)
    add_agent!(pos, M1_cytokines, model)
end


# Create metabolism functions - glycolysis and OXPHOS

function glycolysis!(agent, model) 
    isempty(agent.mitochondria) && return
    # Collect surrounding reactants (glucose and oxygen)
    contenders = collect(nearby_agents(agent, model, 3))
    reactants = filter(m -> m isa Glucose || m isa Oxygen, contenders)
    isempty(reactants) && return 

    # Count how many Glucose and oxygen are around 
    n_Glucose = count(m -> m isa Glucose, reactants)
    n_oxygen = count(m -> m isa Oxygen, reactants)

    # Only proceed if 1 Glucose and 6 oxygen agents are present, 
    # then remove agents
    if n_Glucose >= 1 && n_oxygen >= 6
        Glucose_targets = filter(m -> m isa Glucose, reactants)[1:1]
        oxygen_targets = filter(m -> m isa Oxygen, reactants)[1:6]
        for r in vcat(Glucose_targets, oxygen_targets)
            remove_agent!(r, model)
        end
    end
    
    # Add 0.1 energy per mitochondria 
    agent.energy = agent.energy + (0.1 * length(agent.mitochondria))
    
    # Release ROS into agent
    agent.ROS += 1

    # Release M1 cytokine - only macrophages; surrounding agent 
    if agent isa MDM || agent isa Osteomac
        nbrs = collect(nearby_positions(agent.pos, model, 1))
        isempty(nbrs) && return
        pos = rand(abmrng(model), nbrs)
        add_agent!(pos, M1_cytokines, model)
    else return 
    end
end


# Create OXPHOS function 

function OXPHOS!(agent, model)
    isempty(agent.mitochondria) && return
    # Collect surrounding oxygen 
    contenders = collect(nearby_agents(agent, model, 3))
    oxygens = filter(o -> o isa Oxygen, contenders)
    isempty(oxygens) && return 

    # Count how many oxygen molecules are present 
    n_oxygen = count(o -> o isa Oxygen, oxygens)

    # Only proceed if 6 oxygen agents are present
    if n_oxygen >= 6
        oxygen_targets = oxygens[1:6]
        
        for r in vcat(oxygen_targets)
            remove_agent!(r, model)
        end
    end

    # Add sum of mito health values to 'energy' field
    agent.energy = agent.energy + sum(health for (_, health) in agent.mitochondria)

    # Release 1 ROS per mitochondria into ROS field 
    mito_number = length(agent.mitochondria)
    agent.ROS = agent.ROS + mito_number
    
    # Because defective mitochondria release more ROS, if mito_health <= 0.3,
    # release one more ROS per low health mitochondria
    defectives = filter(t -> t[2] <= 0.3, agent.mitochondria)
    extra_ROS = length(defectives)
    agent.ROS = agent.ROS + extra_ROS

    # If mito_health <= 0.1, relase another ROS 
    almost_dead = filter(t -> t[2] <= 0.1, agent.mitochondria)
    moreROS = length(almost_dead)
    agent.ROS = agent.ROS + moreROS
end


# Create a 'mitophagy!' function 
function mitophagy!(agent, model)
    isempty(agent.mitochondria) && return 
    
    # Select least healthy mitochondria 
    smallest_health = minimum(t -> t[2], agent.mitochondria)

    # Only perform if mito health is below 0.1
    smallest_health > 0.1 && return
    
    # Designate that mitochondria as the sick mitochondria 
    sick_list = filter(t -> t[2] == smallest_health, agent.mitochondria)
    
    # If sick_mito has 2 or more elements, choose only one - prevents a mass 
    # delete of mitochondria if all in one cell have health < 0.1
    sick_mito = rand(sick_list)

    # Delete least healthy mitochondria 
    # First find mito as index; 'sick_mito' is just a value 
    idx = findfirst(==(sick_mito), agent.mitochondria)
    deleteat!(agent.mitochondria, idx)

    # Rearrange mitochondria tuples so x's are in order 
    for (i, (_, health)) in enumerate(agent.mitochondria)
         agent.mitochondria[i] = (i, health)
    end
end


# Create a 'mito_fusion!' function 
function mito_fusion!(agent, model)
    # Find out the number of mitochondria
    n_mitochondria = length(agent.mitochondria)
    
    # Only perform if n_mitochondria >= 2
    if n_mitochondria < 2 && return
    
    # Select two mitochondria 
    first_mito = rand(agent.mitochondria)
    second_mito = rand(agent.mitochondria)
    
    # Create a loop where if second_mito = first_mito, finds another random 
    # mitochondria until they are not the same 
    while (second_mito) == (first_mito) 
        second_mito = rand(agent.mitochondria) 
    end

    # Take the two y-values (health score) and average them 
    average_health = (first_mito[2] + second_mito[2]) / 2

    # Remove the two mitochondria 
    # First need to find mitos by index since first_mito and second_mito are
    # just values 
    idxs = sort([findfirst(==(first_mito), agent.mitochondria),
                 findfirst(==(second_mito), agent.mitochondria)])
    deleteat!(agent.mitochondria, idxs)
    
    # Add a tuple into the mitochondria field where: 
        # x is n_mitochondria - 1 
        # y is the average of the two y-values 
    new_mito_count = length(agent.mitochondria)
    push!(agent.mitochondria, ((new_mito_count - 1), average_health))

    # Rearrange mitochondria so x's are in order  
    for (i, (_, health)) in enumerate(agent.mitochondria)
        agent.mitochondria[i] = (i, health)
    end
end
end
    

# Create a 'mito_fission!' function 
function mito_fission!(agent, model)
    isempty(agent.mitochondria) && return
    
    # Select one mitochondria 
    contender = rand(agent.mitochondria)

    # Find the number of mitochondria
    n_mitochondria = length(agent.mitochondria)

    # Add a tuple where;
        # x = highest x-value + 1
        # y = y-value of selected mitochondria 
    push!(agent.mitochondria, ((n_mitochondria + 1), contender[2]))
end


# Create an 'apoptosis!' function 

function apoptosis!(agent, model)
    # If energy <= 1.0, delete agent 
    agent.energy > 1.0 && return

    # Save position befoe removing the ceell
    pos = agent.pos
    
    # Remove dying cell
    remove_agent!(agent, model)

    # Add three cell debris agents at the same position 
    for _ in 1:3
        add_agent!(pos, CellDebris, model)
    end
end


# Create a 'phagocytosis!' function 

function phagocytosis!(agent, model)
    # Must have enough energy to perform phagocytosis 
    agent.energy < 2 && return 

    # Collect nearby agents 
    contenders = collect(nearby_agents(agent, model, 1))
    isempty(contenders) && return

    # Filter for cell debris 
    debris = filter(a -> a isa CellDebris, contenders)
    isempty(debris) && return 

    # Remove all debris found
    for d in debris
        remove_agent!(d, model)
    end

    # Pay energy cost 
    agent.energy -= 2

    # Move plasticity towards M1 polarization state 
    agent.plasticity -= 0.1
    update_oxphos!(agent)
    update_glycolysis!(agent)

    # Release cytokine and ROS in response
    nbrs = collect(nearby_positions(agent.pos, model, 1))
    isempty(nbrs) && return 

    pos = rand(abmrng(model), nbrs)
    add_agent!(pos, M1_cytokines, model)
    add_agent!(pos, ROS, model)
end


# Create a 'transfer!' function 

function transfer!(agent, model)
    isempty(agent.mitochondria) && return
    contenders = collect(nearby_agents(agent, model, 1))

    # Since mitochondria can be transferred from macrophages to any cell, 
    # including other macrophages, no need to filter 

    # Collect the energy levels of nearby cells 
    cells = filter(c -> hasfield(typeof(c), :energy), contenders)
    isempty(cells) && return
    energies = [c.energy for c in cells]

    # Select contender with lowest energy value 
    min_energy = minimum(energies)
    options = filter(c -> c.energy == min_energy, cells)

    # If two or more healths are the same, choose a random one
    sick_cell = rand(options)

    # Find average health of mitochondria within macrophage 
    avg_health = sum(t -> t[2], agent.mitochondria) / length(agent.mitochondria)

    # Select mitochondria with health equal to or closest to 'avg_health'
    _, closest = findmin(t -> abs(t[2] - avg_health), agent.mitochondria)
    chosen = agent.mitochondria[closest]

    # Add tuple = to selected mitochondria in cell if its energy value <= 1.0
    if sick_cell.energy <= 1.0
        push!(sick_cell.mitochondria, chosen)
        
        # Delete tuple from agent 
        deleteat!(agent.mitochondria, closest)

        # Rearrange mitochondria tuples in both agent and cell
        for (i, (_, health)) in enumerate(agent.mitochondria)
         agent.mitochondria[i] = (i, health)
        end

        for (i, (_, health)) in enumerate(sick_cell.mitochondria)
            sick_cell.mitochondria[i] = (i, health)
        end
    else end
end

    
# Create a sensing function

function sense!(agent, model)
    contenders = collect(nearby_agents(agent, model, 3))
    isempty(contenders) && return

    # Priority 1: Move towards cell debris - since apoptotic signals are relased 
    debris = filter(a -> a isa CellDebris, contenders)
    if !isempty(debris)
        # Choose closest debris particle
        target = argmin(a -> abs(a.pos[1] - agent.pos[1]) +
                             abs(a.pos[2] - agent.pos[2]), debris)
        dx = sign(target.pos[1] - agent.pos[1])
        dy = sign(target.pos[2] - agent.pos[2])
        gs = spacesize(model)
        new_pos = (mod1(agent.pos[1] + dx, gs[1]),
                   mod1(agent.pos[2] + dy, gs[2]))
        move_agent!(agent, new_pos, model)
        return
    end

    # Priority 2: Find nearby DAMPs released this step
    recent_damps = filter(a -> a isa DAMPs && hasfield(typeof(a), :released_at) && 
                                 a.released_at == model.tick, contenders)

    # If any recent DAMPs exist, move towards the closest one
    if !isempty(recent_damps)
        # Pick the closest DAMPs agent
        damp = argmin(a -> abs(a.pos[1] - agent.pos[1]) + abs(a.pos[2] - agent.pos[2]), recent_damps)
        dx = sign(damp.pos[1] - agent.pos[1])
        dy = sign(damp.pos[2] - agent.pos[2])

        gs = spacesize(model)
        move_agent!(agent, (mod1(agent.pos[1] + dx, gs[1]),
                            mod1(agent.pos[2] + dy, gs[2])), model)
        # Release cytokine in response
        nbrs = collect(nearby_positions(agent.pos, model, 1))
        isempty(nbrs) && return 
        pos = rand(abmrng(model), nbrs)
        add_agent!(pos, M1_cytokines, model)
        return
    end

    # Priority 3: Move towards ROS agents
    nearby_ROS = filter(r -> r isa ROS, contenders)
    if !isempty(nearby_ROS) 
        target = argmin(a -> abs(a.pos[1] - agent.pos[1]) +
                             abs(a.pos[2] - agent.pos[2]), nearby_ROS)
        dx = sign(target.pos[1] - agent.pos[1])
        dy = sign(target.pos[2] - agent.pos[2])
        gs = spacesize(model)
        new_pos = (mod1(agent.pos[1] + dx, gs[1]),
                   mod1(agent.pos[2] + dy, gs[2]))
        move_agent!(agent, new_pos, model)
        return
    end

    # Priority 4: Move towards lowest energy cell 
    nearby_cells = filter(a -> 
        a isa MDM || a isa Osteomac ||
        a isa Osteoblast || a isa Osteoclast ||
        a isa Osteocyte,
    contenders
    )

    isempty(nearby_cells) && return

    energies = [c.energy for c in nearby_cells]
    min_energy = minimum(energies)
    options = filter(c -> c.energy == min_energy, nearby_cells)

    # If two or more healths are the same, choose a random one
    sick_cell = rand(abmrng(model), options)

    # Move towards most damaged cell if agent isa MDM
    dx = sign(sick_cell.pos[1] - agent.pos[1])
    dy = sign(sick_cell.pos[2] - agent.pos[2])

    gs = spacesize(model)
    move_agent!(agent, (mod1(agent.pos[1] + dx, gs[1]),
                        mod1(agent.pos[2] + dy, gs[2])), model)
end



# Defnining events and their propensities 

attack_propensity = 1.0
internal_attack_propensity = 1.0
abyss_propensity = 1.0
removeROS_propensity = 1.0
transient_move_propenisty = 1.0
active_move_propensity = 1.0
idle_propensity = 1.0
damage_propensity = 1.0
polarize_to_M1_propensity = 1.0
polarize_to_M2_propensity = 1.0
releaseM1cyto_propensity = 1.0
glycolysis_propensity = 1.0 
oxphos_propensity = 1.0 
apoptosis_propensity = 1.0
phagocytosis_propensity = 1.0 
transfer_propensity = 1.0 
sensing_propensity = 1.0 

attack_event = AgentEvent(
    action! = attack!, propensity = attack_propensity,
    types = ROS
)
internal_attack_event = AgentEvent(
    action! = internal_attack!, propensity = internal_attack_propensity,
    types = Union{MDM, Osteomac, Osteoblast, Osteoclast, Osteocyte}
)
abyss_event = AgentEvent(
    action! = abyss!, propensity = abyss_propensity,
    types = Union{Glucose, Oxygen, ROS, DAMPs, M1_cytokines, M2_cytokines}
)
removeROS_event = AgentEvent(
    action! = removeROS!, propensity = removeROS_propensity,
    types = Union{MDM, Osteomac, Osteocyte, Osteoblast, Osteoclast}
)
active_move_event = AgentEvent(
    action! = active_move!, propensity = active_move_propensity, 
    types = Union{MDM, Osteoblast, Osteoclast, Osteocyte}
)
transient_move_event = AgentEvent(
    action! = transient_move!, propensity = transient_move_propenisty,
    types = Union{ROS, Glucose, Oxygen, CellDebris, M1_cytokines, M2_cytokines, DAMPs}
)
idle_event = AgentEvent(
    action! = idle!, propensity = idle_propensity,
    types = Union{Osteomac, FractureUnit}
)
damage_event = AgentEvent(
    action! = damage!, propensity = damage_propensity,
    types = FractureUnit
)
polarize_to_M1_event = AgentEvent(
    action! = M1_polarization!, propensity = polarize_to_M1_propensity,
    types = Union{M1_cytokines, DAMPs}
)
polarize_to_M2_event = AgentEvent(
    action! = M2_polarization!, propensity = polarize_to_M2_propensity,
    types = M2_cytokines
)
release_M1_event = AgentEvent(
    action! = releaseM1cyto!, propensity = releaseM1cyto_propensity,
    types = Union{MDM, Osteomac}
)
glycolysis_event = AgentEvent(
    action! = glycolysis!, propensity = glycolysis_propensity,
    types = Union{Osteoblast, Osteoclast, Osteocyte}
)
OXPHOS_event = AgentEvent(
    action! = OXPHOS!, propensity = oxphos_propensity, 
    types = Union{Osteoblast, Osteoclast, Osteocyte}
)
macro_OXPHOS_event = AgentEvent(
    action! = OXPHOS!, propensity = (o, model) -> o.oxphos_propensity,
    types = Union{MDM, Osteomac}
)
macro_glycolysis_event = AgentEvent(
    action! = glycolysis!, propensity = (o, model) -> o.glycolysis_propensity,
    types = Union{MDM, Osteomac}
)
mitophagy_event = AgentEvent(
    action! = mitophagy!, 
    propensity = (o, model) -> o.mitophagy_propensity,
    types = Union{MDM, Osteomac, Osteoclast, Osteoblast, Osteocyte}
)
mito_fusion_event = AgentEvent(
    action! = mito_fusion!,
    propensity = (o, model) -> o.mito_fusion_propensity,
    types = Union{MDM, Osteomac, Osteoclast, Osteoblast, Osteocyte}
)
mito_fission_event = AgentEvent(
    action! = mito_fission!, 
    propensity = (o, model) -> o.mito_fission_propensity,
    types = Union{MDM, Osteomac, Osteoclast, Osteoblast, Osteocyte}
)
apoptosis_event = AgentEvent(
    action! = apoptosis!, propensity = apoptosis_propensity,
    types = types = Union{MDM, Osteomac, Osteoclast, Osteoblast, Osteocyte}
)
phagocytosis_event = AgentEvent(
    action! = phagocytosis!, propensity = phagocytosis_propensity,
    types = Union{MDM, Osteomac}
)
transfer_event = AgentEvent(
    action! = transfer!, propensity = transfer_propensity,
    types = Union{MDM, Osteomac, Osteocyte}
)
sensing_event = AgentEvent(
    action! = sense!, propensity = sensing_propensity,
    types = MDM
)


events = (attack_event, internal_attack_event, abyss_event, removeROS_event, active_move_event, transient_move_event, idle_event, damage_event,
polarize_to_M1_event, polarize_to_M2_event, release_M1_event, glycolysis_event, macro_OXPHOS_event, macro_glycolysis_event, OXPHOS_event, mitophagy_event, 
mito_fusion_event, mito_fission_event, apoptosis_event, phagocytosis_event, transfer_event, sensing_event)


# Create EventQueueABM model 

space = GridSpace((50, 50))
rng = Xoshiro(13)

model = EventQueueABM(agents, events, space; rng, warn = false)

# Create a fracture line generator function to be embedded in the initializer 

function generate_fracture_line!(model)
    width, height = spacesize(model)

    # Random starting x-pos at the top of grid 
    x = rand(1:width)
    y = height  
    # Random initial width (1-5 cells thick)
    w = rand(1:5)

    # Helper to place a horizonral chunk of fracture units 
    function place_chunk!(x, y, w)
        half = div(w, 2)
        for dx in -half:half
            xx = clamp(x + dx, 1, width)
            pos = (xx, y)
            add_agent!(pos, FractureUnit, model)
        end
    end

    # Place first chunk
    place_chunk!(x, y, w)

    # Move downward until reaching the bottom
    while y > 1
        # Random horizontal drift 
        dx = rand((-1, 0, 1))
        x = clamp(x + dx, 1, width)

        # Move downward
        y -= 1

        # Randomly vary width (1-5)
        w = clamp(w + rand((-1, 0, 1)), 1, 5)

        # Place chunk for this row 
        place_chunk!(x, y, w)
    end
end


function initialize_bone_model(; x = 50, y = 50, seed = 13,
                                 n_glucose = 200,
                                 n_oxygen = 800)
    space = GridSpace((x, y))
    properties = Dict(
        :n_glucose => n_glucose,
        :n_oxygen => n_oxygen
    )
    model = EventQueueABM(agents, events, space; rng, properties = properties, warn = false)
    gs = spacesize(model)
    
    # Define osteon cluster centers 
    n_osteons = 4 
    osteon_centers = [(rand(abmrng(model), 4:gs[1]-4),
                       rand(abmrng(model), 4:gs[2]-4)) for _ in 1:n_osteons]
    osteon_radii = [rand(abmrng(model), 3:5) for _ in 1:n_osteons]

    # Helper: is a position inside a given osteon's interior 
    function in_interior(pos, center, radius)
        d = sqrt((pos[1]-center[1])^2 + (pos[2]-center[2])^2)
        return d < radius - 1 # strictly inside
    end

    # Helper: is a position on the ring of an osteon
    function on_ring(pos, center, radius)
        d = sqrt((pos[1]-center[1])^2 + (pos[2]-center[2])^2)
        return radius -1 <= d <= radius + 0.8 
    end

    # Helper: is a position in the interstitial space (outside all osteons)
    function is_interstitial(pos)
        for (c, r) in zip(osteon_centers, osteon_radii)
            if sqrt((pos[1]-c[1])^2 + (pos[2]-c[2])^2) <= r + 1
                return false 
            end
        end
        return true 
    end

    # Collect candidate positions for each role 
    all_positions = [(xi, yi) for xi in 1:gs[1], yi in 1:gs[2]]

    interior_positions = []
    ring_positions = []
    interstitial_positions = []

    for pos in all_positions
        for (c, r) in zip(osteon_centers, osteon_radii)
            if in_interior(pos, c, r)
                push!(interior_positions, pos)
                break
            elseif on_ring(pos, c, r)
                push!(ring_positions, pos)
                break
            end
        end
        if is_interstitial(pos)
            push!(interstitial_positions, pos)
        end
    end

    shuffle!(abmrng(model), interior_positions)
    shuffle!(abmrng(model), ring_positions)
    shuffle!(abmrng(model), interstitial_positions)

    # Place MDMs in interstitial space, branching between osteons
    n_MDM = min(20, length(interstitial_positions))
    for i in 1:n_MDM
        pos = interstitial_positions[i]
        mitos = [(i, rand(abmrng(model), 0.1:0.01:1.0)) for i in 1:5]
        energy0 = sum(h for (_, h) in mitos)
        
        add_agent!(pos, MDM, model;
            mitochondria = mitos,
            energy = energy0,
            plasticity = 0.5,
            mitophagy_propensity = 1.0,
            mito_fusion_propensity = 1.0,
            mito_fission_propensity = 1.0,
            oxphos_propensity = 0.1*(100)^0.5,
            glycolysis_propensity = 1.7783*0.003162^0.5,
            ROS = 0,
        )
    end

    # Place osteoblasts on osteon rings 
    n_osteoblasts = min(20, length(ring_positions))
    osteoblast_positions = ring_positions[1:n_osteoblasts]
    for pos in osteoblast_positions
        mitos = [(i, rand(abmrng(model))) for i in 1:5]
        energy0 = sum(h for (_, h) in mitos)

        add_agent!(pos, Osteoblast, model;
            mitochondria = mitos,
            energy = energy0
        )
    end

    # Place osteoclasts on osteo rings (remaining ring spots)
    n_osteoclasts = min(20, length(ring_positions) - n_osteoblasts)
    if n_osteoclasts > 0
        osteoclast_positions = ring_positions[n_osteoblasts + 1 : n_osteoblasts+n_osteoclasts]
        for pos in osteoclast_positions
            mitos = [(i, rand(abmrng(model))) for i in 1:5]
            energy0 = sum(h for (_, h) in mitos)

            add_agent!(pos, Osteoclast, model;
                mitochondria = mitos,
                energy = energy0
            )
        end
    end
    
    # Place osteomacs sharing/adjacent to osteoblast positions
    n_osteomacs = min(25, length(osteoblast_positions))
    for i in 1:n_osteomacs
        base_pos = osteoblast_positions[i]

        # Try to place on the same cell or an adjacent cell
        candidates = [base_pos,
                      (clamp(base_pos[1]+1, 1, gs[1]), base_pos[2]),
                      (clamp(base_pos[1]-1, 1, gs[1]), base_pos[2]),
                      (base_pos[1], clamp(base_pos[2]+1, 1, gs[2])),
                      (base_pos[1], clamp(base_pos[2]-1, 1, gs[2]))]
        pos = rand(abmrng(model), candidates)
        mitos = [(i, rand(abmrng(model))) for i in 1:5]
        energy0 = sum(h for (_, h) in mitos)
        
        add_agent!(pos, Osteomac, model;
            mitochondria = mitos,
            energy = energy0,
            plasticity  = 0.5,
            mitophagy_propensity = 1.0, 
            mito_fusion_propensity = 1.0, 
            mito_fission_propensity = 1.0,
            oxphos_propensity = 0.1*(100)^0.5,
            glycolysis_propensity = 1.7783*0.003162^0.5,
            ROS = 0
        )
    end

    # Place osteocytes in osteon interiors
    n_osteocytes = min(40, length(interior_positions))
    for i in 1:n_osteocytes
        pos = interior_positions[i]
        mitos = [(i, rand(abmrng(model))) for i in 1:5]
        energy0 = sum(h for (_, h) in mitos)

        add_agent!(pos, Osteocyte, model;
            mitochondria = mitos,
            energy = energy0,
        )
    end
    
    # Add molecular agents 
    # Glucose 
    for _ in 1:model.n_glucose
         pos = (rand(abmrng(model), 1:gs[1]),
               rand(abmrng(model), 1:gs[2]))
         add_agent!(pos, Glucose, model)
    end

    # Oxygen 
    for _ in 1:model.n_oxygen
         pos = (rand(abmrng(model), 1:gs[1]),
               rand(abmrng(model), 1:gs[2]))
         add_agent!(pos, Oxygen, model)
    end

    # Fracture units 
    generate_fracture_line!(model)
    
    # Don't add M1 cytokines, DAMPs, cell debris or ROS since cells produce these
    # M2 cytokines will be added by the user while model is running
    return model 
end


# Capture the returned model

bone_repair = initialize_bone_model()

# Test for bugs 
step!(bone_repair, 22.3)


# Data Visualization 
using InteractiveDynamics 
using GLMakie

# Color, shape, and size scheme 
agent_color(a::FractureUnit) = :red2
agent_color(a::MDM)          = :yellowgreen
agent_color(a::Osteomac)     = :chartreuse4
agent_color(a::Osteoblast)   = :royalblue
agent_color(a::Osteoclast)   = :hotpink2
agent_color(a::Osteocyte)    = :powderblue
agent_color(a::Glucose)      = :darkorange1
agent_color(a::Oxygen)       = :deepskyblue3
agent_color(a::ROS)          = :purple2
agent_color(a::M1_cytokines) = :crimson
agent_color(a::M2_cytokines) = :steelblue3
agent_color(a::DAMPs)        = :brown
agent_color(a::CellDebris)   = :gray

agent_marker(a::FractureUnit) = :rect
agent_marker(a::MDM)          = :circle
agent_marker(a::Osteomac)     = :circle
agent_marker(a::Osteoblast)   = :circle
agent_marker(a::Osteoclast)   = :circle
agent_marker(a::Osteocyte)    = :circle
agent_marker(a::Glucose)      = :hexagon
agent_marker(a::Oxygen)       = :star4
agent_marker(a::ROS)          = :star8
agent_marker(a::M1_cytokines) = :xcross
agent_marker(a::M2_cytokines) = :cross
agent_marker(a::DAMPs)        = :diamond
agent_marker(a::CellDebris)   = :star8

agent_size(a::FractureUnit) = 8
agent_size(a::MDM)          = 14
agent_size(a::Osteomac)     = 14
agent_size(a::Osteoblast)   = 13
agent_size(a::Osteoclast)   = 13
agent_size(a::Osteocyte)    = 12
agent_size(a::Glucose)      = 6
agent_size(a::Oxygen)       = 4
agent_size(a::ROS)          = 7
agent_size(a::M1_cytokines) = 5
agent_size(a::M2_cytokines) = 5
agent_size(a::DAMPs)        = 5
agent_size(a::CellDebris)   = 5

# Slider to change glucose and oxygen levels 
params = Dict(
    :n_glucose => 50:50:500,
    :n_oxygen => 200:100:1200
)

fig, ax, abmobs = abmplot(bone_repair;
    agent_color,
    agent_marker,
    agent_size,
    add_controls = true, 
    params = params,
    figure = (; size = (900, 900)),
    axis = (; title = "Bone Repair Model"),
)

# Slider/button that allows user to add M2 cytokines to illustrate potential therapuetic effects
# for resolving inflammation 
add_m2!(model) = add_agent!(
    rand(1:spacesize(model)[1]), rand(1:spacesize(model)[2]),
    M2_cytokines, model
)

# Add an M2 cytokine button to figure 
m2_button = Button(fig[3, 1]; label = "Add M2 Cytokines", tellwidth = false)
on(m2_button.clicks) do _
    for _ in 1:10 # adds burst of 10 M2 cytokines per click 
        add_m2!(abmobs.model[])
    end
end

display(fig)
