# Author: Aaron Dunton (ajdunton@gmail.com)

# %% Import Dependencies
import geopandas as gpd
import rasterio
import pandas as pd
from shapely.geometry import LineString, Point
import os

import funcs as fn
import params as pr

# %% Topology

# Open street shapefile and some cleaning
streets = fn.street_cleaning()

# Coordinate reference system
coref = streets.crs

# Open terminal node (i.e., outlet location) shapefile
termnode = gpd.read_file("./input/termnode/termnode.shp")

# Open building footprints
buildings = gpd.read_file("./input/buildings/buildings.shp")

# Assign each building to a street segment
buildings["street_index"] = fn.assign_to_street(buildings,streets)

# Nodes. ww_nodes has column 'attached_buildings'. sw_nodes includes both 
# corner and intermediary nodes
ww_nodes, sw_nodes = fn.nodes_and_split_streets(buildings,streets)

# Topology: shortest path spanning tree of the street network (i.e., no extra 
# nodes)
t, t_pos, tn_id = fn.tree_topology(termnode)

# Save the tree in ./inter
geoms = []
for i,j in t.edges():
    pt1 = Point(t_pos[i])
    pt2 = Point(t_pos[j])
    geoms.append(LineString((pt1,pt2)))
t_gdf = gpd.GeoDataFrame(geometry = geoms).set_crs(coref)
t_gdf.to_file('./inter/t_'+pr.weight_method)

# Full graph with wastewater nodes and design graph
ww_pos, ww_spst, design_graph = fn.graphs_from_t(t,t_pos)

# Add graph index to the nodes GeoDataFrame
ww_nodes['graph_node_ind'] = [fn.identify_ind(i,ww_pos) 
                              for i in ww_nodes['geometry']]

sw_nodes['graph_node_ind'] = [fn.identify_ind(i,ww_pos) 
                              for i in sw_nodes['geometry']]
sw_nodes.drop_duplicates(subset='graph_node_ind',inplace=True)
sw_nodes.reset_index(drop=True,inplace=True)

# %% Estimating Stormwater Inflow

# Open slope raster and extract array and affine transformation
with rasterio.open('./input/slope/slope.tif') as slope:
    slope_array = slope.read(1)
    slope_affine = slope.transform

# Open the percent impervious raster and extract array and affine
with rasterio.open('./input/per_imp/per_imp.tif') as per_imp:
    imp_array = per_imp.read(1)
    imp_affine = per_imp.transform

# Total area for sub-area
tot_area = gpd.read_file('./input/area/area.shp')

# Properties for each row of sw_nodes                                                    
sw_nodes['voroni_areas'] = fn.voronis(sw_nodes['geometry'], 
                                      tot_area['geometry'][tot_area.index[0]])
sw_nodes['area'] = sw_nodes.apply(fn.area_column,axis=1)
sw_nodes['width'] = sw_nodes.apply(fn.width_column,axis=1)
sw_nodes['slope'] = [fn.slope_column(geom,slope_array,slope_affine) 
                     for geom in sw_nodes['voroni_areas']]
sw_nodes['percent_impervious'] = [fn.imp_column(geom,imp_array,imp_affine) 
                                  for geom in sw_nodes['voroni_areas']]
sw_nodes['c10'] = sw_nodes.apply(fn.runoffc10, axis=1, HSG=pr.soil_HSG)
sw_nodes['nodal_inflow_ft3ps'] = sw_nodes.apply(fn.cia, axis=1)

# Save the nodal inflow to the graph
for node_num, flow in zip(sw_nodes['graph_node_ind'],
                          sw_nodes['nodal_inflow_ft3ps']):
    ww_spst.nodes[node_num]['sw_ft3ps'] = flow

# Accumulate nodal stormwater inflow over the graph
ww_spst = fn.accumulate(ww_spst, 'sw_ft3ps')

# %% Estimating Wastewater Inflow
# Building-level data
buildings['des_pop'] = buildings.apply(fn.ww_pop,axis=1)
buildings['ci_gpd'] = buildings.apply(fn.ww_ci,axis=1)

# Node-level data
ww_nodes['des_pop'] = [sum(buildings.iloc[bldgs]['des_pop']) 
                       for bldgs in ww_nodes['attached_buildings']]
ww_nodes['ci_gpd'] = [sum(buildings.iloc[bldgs]['ci_gpd']) 
                      for bldgs in ww_nodes['attached_buildings']]

# Save node-level data to graph
for node_num, des_pop, ci_gpd in zip(ww_nodes['graph_node_ind'],
                                     ww_nodes['des_pop'],
                                     ww_nodes['ci_gpd']):
    ww_spst.nodes[node_num]['des_pop'] = des_pop
    ww_spst.nodes[node_num]['ci_gpd'] = ci_gpd
    
# Accumulate values over tree
ww_spst = fn.accumulate(ww_spst,'des_pop')
ww_spst = fn.accumulate(ww_spst,'ci_gpd')

# Evaluate design wastewater flow, using peak factor
for i in ww_spst.nodes():
    pop = ww_spst.nodes[i]['des_pop_acc']
    q_gpd = pop*pr.ww_percap*fn.peak_f(pop)+ww_spst.nodes[i]['ci_gpd_acc']
    q_ft3ps = q_gpd/(7.48052*86400)
    ww_spst.nodes[i]['ww_ft3ps_acc'] = q_ft3ps
 

# %% Hydraulic Design
    
# Add stormwater and wastewater design flows to the design graph
for i,j in design_graph.edges():
    design_graph.edges[(i,j)]['sq'] = ww_spst.nodes[i]['sw_ft3ps_acc']

    ww_node = design_graph.edges[(i,j)]['path'][-2]
    design_graph.edges[(i,j)]['wwq'] = ww_spst.nodes[ww_node]['ww_ft3ps_acc']   

# Pipe inventory
pipes = pd.read_csv('./input/pipe_sizes.csv')

# Pipe Sizing
for i,j in design_graph.edges():
    des_q = max(design_graph.edges[(i,j)]['sq'], 
                design_graph.edges[(i,j)]['wwq'])
    
    # Identify subset of pipe sizes that meet capacity constraint
    s = pipes[pipes.cap_f3ps>des_q]
    
    # Choose the minimum diameter among the subset of sizes
    if len(s)!=0:
        diam = s.diam_in.min()
    else:
        diam = pipes.diam_in.max()
    
    # Add diameter and corresponding slope to the design_graph edge
    design_graph.edges[(i,j)]['diam_in'] = diam
    design_graph.edges[(i,j)]['slope'] = \
        float(pipes[pipes.diam_in==diam].slope_fpf)  

# %% Evaluate Pipe and Nodal Depths
# Add elevation data to the nodes in the design graph
with rasterio.open('./input/dem_ft/dem_ft.tif') as dem:  
    for i in design_graph.nodes():
        elevation = next(rasterio.sample.sample_gen(dem,[ww_pos[i]]))
        design_graph.nodes[i]['elev'] = elevation

# Add length data to the edges in the design graph
for i,j in design_graph.edges():
    design_graph.edges[(i,j)]['length'] = fn.coord_dist(ww_pos[i],ww_pos[j])
    # Here, assuming straight segments between nodes

design_graph = fn.depths(design_graph)

# %% Prepare Outputs
# Data from ww_nodes, add nodetype
ww_n_out = ww_nodes[['graph_node_ind','geometry']].copy()
ww_n_out['nodetype'] = 'waste'

# Copy data from storm_nodes, add nodetype
sw_n_out = sw_nodes[['graph_node_ind','area','width','slope',
                     'percent_impervious','c10','geometry']].copy()
sw_n_out['nodetype'] = 'storm'

# Combine storm and waste nodes into one gdf
outnodes = sw_n_out.append(ww_n_out).set_crs(coref)

# Rename columns for output 
new_names = {'graph_node_ind':'ind','percent_impervious':'perc_imp'}
outnodes.rename(new_names, axis=1, inplace=True)

# Add nodal depths for the nodes that are in design_graph
for i in design_graph.nodes():
    outnodes.loc[outnodes.ind==i,'depth'] = design_graph.nodes[i]['depth']
    
# Add nodal ground elevation
def sample_raster_at_gdf_row(row,rast):
    val = next(rasterio.sample.sample_gen(rast,
                                          [(row.geometry.x,row.geometry.y)]))
    return val[0]

with rasterio.open('./input/dem_ft/dem_ft.tif') as dem:
    outnodes['ground_el'] = outnodes.apply(sample_raster_at_gdf_row,
                                           axis=1,
                                           rast=dem)    

# Add indicator for outlet node
outnodes['outlet'] = 0
outnodes.loc[(outnodes.ind==tn_id),'outlet'] = 1

# Create outedges gdf
outedges = gpd.GeoDataFrame(columns=['edge','node_0','node_1','depth_0',
                                     'depth_1','diameter','slope','geometry'])
outedges.set_crs(coref,inplace=True)
outedges['edge'] = [x for x in ww_spst.edges()]
outedges['node_0'] = [i for i,j in ww_spst.edges()]
outedges['node_1'] = [j for i,j in ww_spst.edges()]

# Transfer data from design_graph to outedges
for i,j in design_graph.edges():
    path = design_graph.edges[(i,j)]['path']
    
    # If there are no intermediate nodes along (i,j), data transfer is trivial
    if len(path)==2:
        outedges.loc[outedges.edge==(i,j),'diameter'] = \
            design_graph.edges[(i,j)]['diam_in']
        outedges.loc[outedges.edge==(i,j),'slope'] = \
            design_graph.edges[(i,j)]['slope']
        outedges.loc[outedges.edge==(i,j),'depth_0'] = \
            design_graph.edges[(i,j)]['depth0']
        outedges.loc[outedges.edge==(i,j),'depth_1'] = \
            design_graph.edges[(i,j)]['depth1']
        outedges.loc[outedges.edge==(i,j),'geometry'] = \
            LineString([ww_pos[i],ww_pos[j]])
    
    # If there are intermediate nodes along (i,j)
    else:
        di = design_graph.edges[(i,j)]['depth0']
        dj = design_graph.edges[(i,j)]['depth1']
        lij = fn.coord_dist(ww_pos[i],ww_pos[j]) 
        
        # Loop over all edges from i to j
        cur_dn = di
        for n,m in zip(path[:-1],path[1:]):
            
            # Calculate depths of intermediate nodes by linear interpolation
            cur_dm = di+(dj-di)*(fn.coord_dist(ww_pos[i],ww_pos[m])/lij)
            
            # Store data to outedges
            outedges.loc[outedges.edge==(n,m),'diameter'] = \
                design_graph.edges[(i,j)]['diam_in']
            outedges.loc[outedges.edge==(n,m),'slope'] = \
                design_graph.edges[(i,j)]['slope']
            outedges.loc[outedges.edge==(n,m),'depth_0'] = cur_dn
            outedges.loc[outedges.edge==(n,m),'depth_1'] = cur_dm
            outedges.loc[outedges.edge==(n,m),'geometry'] = \
                LineString([ww_pos[n],ww_pos[m]])
            
            # Also store depth for intermediate nodes to outnodes
            if m!=j:
                outnodes.loc[outnodes.ind==m,'depth'] = cur_dm
                
            cur_dn = cur_dm

# Remove edge tuple from output GDF (cannot save to .shp attribute table)
outedges.drop('edge',axis=1,inplace=True)

# Buildings
try: 
    outbldgs = buildings[['population','type','nunits',
                          'area_ft2','geometry']].copy()
except KeyError:
    outbldgs = buildings[['population','geometry']].copy()
outbldgs['node_ind'] = ""
for ind, b_list in zip(ww_nodes.graph_node_ind, ww_nodes.attached_buildings):
    for b in b_list:
        outbldgs.iloc[b, outbldgs.columns.get_loc('node_ind')] = ind


# %% Save Outputs
if not os.path.isdir('output'):
    os.mkdir('output')
    
outnodes.to_file('output/outnodes')

types_edges = {'node_0':int,'node_1':int,'diameter':int,'slope':float,
               'depth_0':float,'depth_1':float}
outedges.astype(types_edges).to_file('output/outedges')

outbldgs.to_file('output/outbldgs')