# Fracture Model
This software models the bone tissue microenvironment following fracture in relation to nutrients availability and molecular signaling patterns, while demonstrating mitochondrial transfer dynamics and intracellular mitochodrial homeostasis. 

## How this model works
This model is written entirely in Julia, using its [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/) package. Cellular and molecular elements are created as agents, each designated a set of functions that are performed within each step. The list of agents used in this model are listed below with their respective representative icons.

<img width="234" height="320" alt="AgentsLegend" src="https://github.com/user-attachments/assets/214177d3-c1fb-47a1-92b0-4dd11a244ae3" />

The model operates in continuous time through the package's EventQueueABM. Events are defined by their propensities and the types of agents that perform that function are specified within the 'types' field; this can be found after the created functions and before initialization. The possible functions performed are described below:

