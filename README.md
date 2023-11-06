# Capillary Network Generator

This program creates a representation of *capillary* infrastructure, which we define as the small-scale infrastructure that connect to every building. Specifically in the implemented program, we generate combined sewer networks. Modifications may be made to generate other types of capillary infrastructure.

# Inputs
The inputs, in the input directory, are listed below. Included in this repository are the inputs for a case study location in Seaside, Oregon. See [REF] for details on sources where these data can be collected for other locations. Note that all geospatial inputs need to be in the same projected coordinate reference system, with units in feet.
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
  <img src="https://github.com/ajdunton/capillary-network-generator/assets/147078788/5ff7e8f5-b79f-48ab-943b-6bf9c29b620b" width="300">
</p>
<p align="center">
  Output Network Nodes and Edges
</p>


# Procedure
The following is a brief overview of the 3-step procedure that is implemented in script.py. Note that params.py and funcs.py are modules that are imported into script.py. These modules contain parameters/constants and functions, respectively. 

## 1. Generate the Network Topology
To generate the topology, we first assign each building to a street segment and determine wastewater/demand nodes based on the assignment. For combined sewer networks, there are also nodes for stormwater inflow. We determine a tree topology for the combined sewer network as an optimal spanning tree of the wastewater and stormwater nodes.

## 2. Estimate the Resource Demand at Each Node
We use either a building-based or population-based method to estimate the inflow at each wastewater node. The building-based method is based on the footprint area of each building and the population-based method uses a population attribute saved in the input buildings shapefile. That is, if using the population-based method, the building footprints shapefile must have a population column.

For stormwater inflow, we use Voronoi polygons to determine geographic areas associated with each stormwater node. Then, we estimate the inflow using the rational method.

## 3. Choose Components
We choose each pipe from a set of standard available pipes, each with a corresponding slope and flow capacity. We choose the smallest pipe with adequate flow capacity.

# Reference
See the following paper for further details about this network generator and validation for a case study locaiton.

[Citation]
