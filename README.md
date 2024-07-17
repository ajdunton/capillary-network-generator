# Capillary Network Generator
Author: Aaron Dunton (ajdunton@gmail.com)

This program creates a representation of *capillary* infrastructure, which we define as the small-scale infrastructure that connect to every building. Specifically in the implemented program, we generate combined sewer networks. Modifications may be made to generate other types of capillary infrastructure.

# Inputs
The inputs, in the input directory, are listed below. Included in this repository are the inputs for a case study location in Seaside, Oregon. See Dunton and Gardoni (2023) for details on sources where these data can be collected for other locations. Note that all geospatial inputs need to be in the same projected coordinate reference system, with units in feet.
- area: shapefile that identifies the extent of the study area
- buildings: shapefile of building footprints 
- dem_ft: elevation raster, in feet
- per_imp: percent impervious area raster
- slope: ground slope raster
- streets: streets shapefile with every road segment that leaves the area removed, except the road segment that leads to the terminal node
- termnode: terminal node shapefile (i.e., where the capillary network is directed towards/from)
- pipe_sizes.csv: set of available pipe diameters with corresponding slopes and flow capacities

# Outputs
The script produces the following outputs that define the capillary network.
- outbldgs: shapefile of building footprints, including corresponding node for each building
- outedges: shapefile of network edges, including the diameter and slope for each pipe/edge
- outnodes: shapefile of network nodes, including attributes that define inflow at each node

<p align="center">
  <img src="https://github.com/ajdunton/capillary-network-generator/assets/147078788/74bebab9-113a-4f0d-9395-021f0c281d5d" width="400">
</p>
<p align="center">
  Output Network Nodes and Edges
</p>

# Procedure
The following is a brief overview of the 3-step procedure that is implemented in script.py. Note that params.py and funcs.py are modules that are imported into script.py. These modules contain parameters/constants and functions, respectively. 

## 1. Generate the Network Topology
To generate the topology, we first assign each building to a street segment and determine wastewater/demand nodes based on the assignment. For combined sewer networks, there are also nodes for stormwater inflow. We determine a tree topology for the combined sewer network as an optimal spanning tree of the wastewater and stormwater nodes. We use the elevation to determine weights for each edge. See Dunton and Gardoni (2023) where different formulations of this weight were tested. Comparison with the network topology for a case study location concluded that the Hybrid-Elevation formulation yields realistic results. This option is specified by ww_method = 'hybrid-elev' in ./params.py.

## 2. Estimate the Resource Demand at Each Node
We use either a building-based or population-based method to estimate the inflow at each wastewater node. The method to be used is specified by setting ww_method as 'buildings' or 'population' in ./params.py. The building-based method is based on the footprint area of each building and the population-based method uses a population attribute saved in the input buildings shapefile. That is, if using the population-based method, the building footprints shapefile must have a population column.
For stormwater inflow, we use Voronoi polygons to determine geographic areas associated with each stormwater node. Then, we estimate the inflow using the standard rational method.

## 3. Choose Components
We choose each pipe from a set of standard available pipes, each with a corresponding slope and flow capacity. We choose the smallest pipe with adequate flow capacity.

# Reference
See the following paper for further details about this procedure and validation for a case study.

Dunton, A., and Gardoni, P. (2024). Generating network representations of small-scale infrastructure using generally available data. *Computer-Aided Civil and Infrastructure Engineering*, 39(8):1143-1158. https://doi.org/10.1111/mice.13137
