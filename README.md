# Capillary Network Generator

This program creates a representation of *capillary* infrastructure, which we define as the small-scale infrastructure that connect to every building. Specifically in the implemented program, we generate combined sewer networks. Modifications may be made to generate other types of capillary infrastructure.

# Inputs
The inputs, in the input directory, are listed below. Included here are the inputs for a case study location in Seaside, Oregon. See [REF] for details on sources where these data can be collected for other locations. Note that for all geospatial inputs (shapefiles for vector; geotiff for raster) need to be in a consistent projected coordinate reference system, with units in feet.
- area: shapefile that identifies the extent of the study area
- buildings: shapefile of building footprints 
- dem_ft: elevation raster, in feet
- per_imp: percent impervious area raster
- slope: ground slope raster
- streets: streets shapefile with every road segment that leaves the area removed, except the road segment that leads to teh terminal node
- termnode: terminal node shapefile (i.e., where the capillary network is directed towards/from)
- pipe_sizes.csv: set of available pipe diameters with corresponding slopes and flow capacities

# Outputs
The script produces the following outputs that define the capillary network.
- outbldgs: building shapefile that identifies corresponding node for each building
- outedges: shapefile of edges of the network, including diameter, slope, edge index, and from/to node indices
- outnodes: shapefile of the nodes of the network, including node type, attributes that define the inflow at the node, and node index

[output_network.pdf](https://github.com/ajdunton/capillary-network-generator/files/13270543/output_network.pdf)

# Procedure
The following is a brief overview of the procedure that is implemented in script.py. Note that params.py and funcs.py are modules that are imported into script.py. These modules contain parameters/constants and functions, respectively. 

## Generate the Network Topology

## Estimate the Resource Demand at Each Node

Note that if using the population-based method, the necessary data needs to be included in the building footprints shapefile. 
## Choose Components

# Reference
See the following paper for further details about this network generator and validation for a case study locaiton.

[Citation]
