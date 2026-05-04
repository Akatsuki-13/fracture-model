# Fracture Model
This software models the bone tissue microenvironment following fracture in relation to nutrients availability and molecular signaling patterns, while demonstrating mitochondrial transfer dynamics and intracellular mitochodrial homeostasis. 

## How this model works
This model is written entirely in Julia, using its [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/) package. Cellular and molecular elements are created as agents, each designated a set of functions that are performed within each step. The list of agents used in this model are listed below with their respective representative icons.

<img width="234" height="320" alt="AgentsLegend" src="https://github.com/user-attachments/assets/214177d3-c1fb-47a1-92b0-4dd11a244ae3" />



Since agents are defined as structs, specific fields have been created within each cell type. The **mitochondria** field is a vector of tuples where the x-value is the id assigned to the mitochondrion and the y-value is the health value assigned to the mitochondrion. The **energy** field is 

The model operates in continuous time through the package's EventQueueABM. Events are defined by their propensities and the types of agents that perform that function are specified within the 'types' field; this can be found after the created functions and before initialization. The possible functions performed are described below:

-**Attack:** ROS agent reduces health of mitochondria within cell; if cell is a macrophage, plasticity moves towards M1 side of the spectrum 

-**Update OXPHOS:** updates the propensity for OXPHOS to occur within macrophages according to its plasticity; the closer to the M2 side of the spectrum, the more likely OXPHOS will be performed

-**Internal attack:** intracellular ROS attacks the health of a random mitochondria within the cell and 'DAMPs' agent is released

-**Abyss:** molecular agents that share a position with a fracture unit will 'disappear' from the model

-**Remove ROS:** intracellular ROS is removed from the cell using one energy point

-**Transient move:** allows molecular agents to move freely throughout the model 

-**Active move:** cellular agents move throughout model using 0.1 energy points 

-**Idle:** keeps 'Osteomac' and 'Fracture Unit' agents idle 

-**Damage:** cellular agents that share a position with a 'Fracture Unit' agent will continuously sustain damage and release a 'DAMPs' agent in response 

-**M1 polarization:** M1 cytokines and DAMPs agents move macrophage plasticity towards M1 side of spectrum

-**M2 polarization:** M2 cytokine agents move macrophage plasticity towards M2 side of spectrum 

-**Release M1 cytokines:** M1 cytokine agent is released by cellular agent in response to interaction with M1 cytokine agent in environment

-**Glycolysis:** cellular agents use one glucose agent to add 0.1 energy points per mitochondria and release one ROS agent and one M1 cytokine, the latter released only by macrophage agents 

-**OXPHOS:** oxidative phosphorylation is performed by cellular agents, where six oxygen agents are taken up to add the sum of all mitochondria health values to the cell's energy field; intracellular ROS are released according to the health of the mitochondria within the cell, the less healthy, the more ROS accumulate

-**Mitophagy:** removes mitochondria with health values below 0.1

-**Mitochondrial fusion:** 'fuses' two mitochondria within a cell; the resulting mitochondrion's health value is the average health value of the two mitochondria 

-**Mitochondrial fission:** creates a copy of a random mitochondrion within the cell

-**Apoptosis:** cell will destroy itself  if its energy level is below 1.0; three 'Cell Debris' agents are released as a result 

-**Phagocytosis:** macrophage agents uptake 'Cell Debris' agents, moving the plasticity towards the M1 side of the spectrum and releasing an M1 cytokine agent and ROS agent

-**Transfer:** cellular agents can transfer a mitochondrion with a health of average value compared to all mitochondria within the cell to another cellular agent

-**Sense:** 'MDM' agents move towards nearby cell debris, DAMPs, or cells of the lowest energy value 
