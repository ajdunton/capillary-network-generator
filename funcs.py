# Author: Aaron Dunton (ajdunton@gmail.com)

import geopandas as gpd
import networkx as nx
import numpy as np
import math
import shapely
import rasterio
import os
from scipy.stats import norm
from shapely.geometry import Point, LineString, MultiPolygon
from shapely.ops import nearest_points, split
from geovoronoi import voronoi_regions_from_coords
from rasterstats import zonal_stats

import params as pr

# Functions used in the small-scale generator script

def graph_from_shp(path):
    """
    Creates an undirected graph and position dictionary from a shapefile, with
    integer node numbers.
    
    Parameters
    ----------
    path : string
        Path to the shapefile.

    Returns
    -------
    X : networkx graph 
        Undirected graph corresponding to the input shapefile.
    pos : dictionary
        Position dictionary for the graph {node #: (coord_x,coord_y)}.
    """
    G = nx.read_shp(path).to_undirected()
    pos = {k: v for k,v in enumerate(G.nodes())}
    X = nx.Graph()
    X.add_nodes_from(pos.keys())
    co_edges = [set(x) for x in G.edges()]
    int_edges = [tuple(x for x,y in pos.items() if y in edg) 
                 for edg in co_edges]
    X.add_edges_from(int_edges)
    return X, pos

# Function to open a .shp as a graph
def digraph_from_shp(path):
    """
    Creates a directed graph and position dictioary from a shapefile, with 
    integer node numbers.
    
    There are two directed edges for each line in the shapefile.

    Parameters
    ----------
    path : string
        Path to the shapefile.

    Returns
    -------
    X : networkx digraph
        Directed graph corresponding to the input shapefile.
    pos : dictionary
        Position dictionary for the graph {node #: (coord_x,coord_y)}.

    """
    G = nx.read_shp(path)
    pos = {k: v for k,v in enumerate(G.nodes())}
    X = nx.DiGraph() 
    X.add_nodes_from(pos.keys())
    l = list(G.edges())
    edg = [tuple(k for k,v in pos.items() if v in sl) for sl in l]
    edg_copy = edg.copy()
    for x,y in edg_copy: edg.append((y,x))
    X.add_edges_from(edg)
    
    return X, pos

def street_cleaning():
    """
    Opens the streets shapefile at ./input/streets; does some cleaning 
    (notably, if not connected, removes all but the largest component); saves 
    the cleaned GDF as ./inter/streets_clean; and returns the streets GDF.

    Returns
    -------
    streets : GeoDataFrame
        GDF of the cleaned streets.

    """
    streets = gpd.read_file('./input/streets/streets.shp')
    streets = streets.loc[streets['geometry'].is_valid, :]
    streets.drop_duplicates(inplace=True)
    streets.reset_index(drop=True, inplace=True)
    
    # Save cleaned streets to re-open with networkx
    if not os.path.isdir('inter'):
        os.mkdir('inter')
    streets.to_file('./inter/streets_clean')
    
    # Remove all but largest component if corresponding graph is not connected
    graph_check, pos_check = \
        graph_from_shp('./inter/streets_clean/streets_clean.shp')
    if not nx.is_connected(graph_check):
        nodes_to_delete = set(graph_check.nodes()) - \
            max(nx.connected_components(graph_check),key=len)
        s_inds = set()
        for i in nodes_to_delete:
            node_point = Point(pos_check[i]).buffer(1)
            for j in streets[streets.intersects(node_point)].index:
                s_inds.add(j)
            s_inds.add(streets.loc[streets.intersects(node_point)].index[0])
        streets.drop(labels=s_inds, axis=0, inplace=True)
        streets.reset_index(drop=True, inplace=True)
        streets.to_file('./inter/streets_clean')  

    return streets

def assign_to_street(bldgs,strts):
    """
    Assigns each building to a street using the rule-based approach.

    Parameters
    ----------
    bldgs : GeoDataFrame
        GeoDataFrame of buildings.
    strts : GeoDataFrame
        GeoDataFrame of streets.

    Returns
    -------
    street_if : list
        Street index (final) for each building, in same order as input bldgs.

    """
    if len(strts)>1:
        # Lists for saving first-order results, in order of buildings
        dists = [None]*len(bldgs)
        street_i1 = [None]*len(bldgs)
        # Call two_links for each building
        for i, bld in enumerate(bldgs['geometry']):
            dists[i],_,street_i1[i] = two_links(bld,strts)
        
        # Number of buildings connected to each street from nearest neighbor  
        n_con1 = dict.fromkeys(strts.index)
        for i in strts.index: 
            n_con1[i] = [x for x,_ in street_i1].count(i)
    
        # Set of streets that will have buildings attached to them
        street_set = set()
        # Follow rule to make decision between two candidates for each building
        for i, sts in enumerate(street_i1):    
            decision = 0
            if dists[i][1]<pr.maxratio*dists[i][0] or dists[i][1]<pr.maxdist:
                decision = np.argmax([n_con1[sts[0]], n_con1[sts[1]]])
            if decision == 0:
                street_set.add(sts[decision])
       
        # If the closest street is in the final street set, then the building is 
        # attached to the closest street. If not, attached to the second-closest
        street_if = [None]*len(bldgs)
        for i, sts in enumerate(street_i1):
            if sts[0] in street_set:
                street_if[i] = sts[0]
            else:
                street_if[i] = sts[1]
    
    else: street_if = [strts.index[0]]*len(bldgs)

    # Return list of street index for each building   
    return street_if

def nodes_and_split_streets(bldgs,strts):
    """
    Function that creates a wastewater node geodataframe, using "street_index" 
    in "bldgs." Also splits the street geodataframe at the wastewater nodes and 
    saves as a shapefile "splitgdf."

    Parameters
    ----------
    bldgs : GeoDataFrame
        Buildings, with column for associated street index for each building.
    strts : GeoDataFrame
        Streets, not yet modified.

    Returns
    -------
    waste_nodes : GeoDataFrame
        Wastewater nodes, with "attached_buildings" column for each node.
    storm_nodes : GeoDataFrame
        Stormwater nodes, geometry only, including both intermediary and corner 
        nodes.
    """
    # Nodes
    split_streets = []
    ww_nodes_geom = []
    ww_nodes_bi = []
    
    sw_nodes_geom = []
    
    for s,ind in zip(strts["geometry"],strts.index):
        # Buildings assigned to the street segment
        b_on_street = bldgs.loc[bldgs['street_index']==ind]
        
        # Add endpoints to list of stormwater node geometries
        sw_nodes_geom.append(s.boundary[0])
        sw_nodes_geom.append(s.boundary[1])
        
        # If there are any intermediary nodes along a street segment, add nodes 
        # to lists for returning and split the street segment and add it to the
        # geometry list for saving as splitgdf.
        if s.length>pr.max_drain_spacing or len(b_on_street)!=0:
            new_sw_nodes = []
            # Add point geometries to list of stormwater node geometries
            if s.length>pr.max_drain_spacing:
                n_new = int(s.length//pr.max_drain_spacing)
                for j in range(n_new):
                    new_pt = s.interpolate((2*(j+1)-1)/(2*n_new)+0.01/n_new,
                                           normalized=True)
                    sw_nodes_geom.append(new_pt)
                    new_sw_nodes.append(new_pt)
                    
            # Add point geometries to list of wastewater node geometries
            newnodes=gpd.GeoDataFrame(geometry=[])
            if len(b_on_street)!=0:
                # Get all new nodes
                newnodes = get_nodes(s, b_on_street, pr.max_node_spacing)
                
                # Drop nodes with no attached buildings
                drop_rows = []
                for ind,row in newnodes.iterrows():
                    if len(row["building_inds"])==0: drop_rows.append(ind)
                newnodes.drop(drop_rows,inplace=True)
                
                # Save remaining newnodes and building_inds attribute
                for i in newnodes["geometry"]: ww_nodes_geom.append(i)
                for i in newnodes["building_inds"]: ww_nodes_bi.append(i)
            
            # Split streetline at new nodes and append to list of lines
            for i in break_line(s,(list(newnodes["geometry"])+new_sw_nodes)):
                split_streets.append(i)
                
        # If there are no intermediary nodes, add the complete street segment
        # to the geometry list for saving as splitgdf.
        else:
            split_streets.append(s)
    
    # Create GDF for wastewater street segments and save as a shapefile
    split_gdf = gpd.GeoDataFrame(geometry=split_streets)
    split_gdf.set_crs(strts.crs,inplace=True)
    split_gdf.to_file("./inter/splitgdf")
    
    # Create wastewater nodes GDF
    waste_nodes = gpd.GeoDataFrame({"attached_buildings": ww_nodes_bi},
                             geometry=ww_nodes_geom)
    
    # Create stormwater nodes GDF
    storm_nodes = gpd.GeoDataFrame(geometry=sw_nodes_geom)

    return waste_nodes, storm_nodes

def coord_dist(p0,p1):
    """
    Calculates the distance between two points.

    Parameters
    ----------
    p0 : Tuple
        Coordinate tuple of point 0.
    p1 : TYPE
        Coordinate tuple of point 1.

    Returns
    -------
    float
        Distance between two points.
    """
    return math.sqrt((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2)

def get_spst(graph,target,keepweights=False,keepdists=False):
    """
    Determines the shortest path spanning tree of an input graph considering a 
    target, using Dijkstra's algorithm.

    Parameters
    ----------
    graph : networkx graph
        Input graph.
    target : int
        Node number corresponding to the target node.
    keepweights : boolean
        If true, store weights from the input graph in the output tree 
        (default False).
    keepdists : boolean
        If true, store the combined path weight from every node to the target 
        as node attribute 'w_to_target'.

    Returns
    -------
    t : networkx directed graph
        Shortest path spanning tree, where edge direction indicates the path 
        from each node to the target node.
    """
    dist, path = nx.single_source_dijkstra(graph,target)
    spst_edges = [(x[-1],x[-2]) for x in path.values() if len(x)>1]
    t = nx.DiGraph()
    t.add_nodes_from(graph.nodes())
    t.add_edges_from(spst_edges)
    if keepweights:
        for i,j in t.edges():
            t[i][j]['weight'] = graph[i][j]['weight']
    if keepdists:
        for i in t.nodes():
            t.nodes[i]['w_to_target'] = dist[i]
    return t

def graph_from_shp_2(path,cur_pos):
    """
    Creates a graph from a shapefile (specified with path), maintaining any 
    existing node numbering from a different graph (by specifying the position 
    dictionary of the existing graph).

    Parameters
    ----------
    path : string
        Location of shapefile to be opened as a graph.
    cur_pos : dictionary
        Existing position dictionary.

    Returns
    -------
    X : NetworkX Graph
        Resulting graph.
    return_pos : dictionary
        Position dictionary for the return graph.

    """
    G = nx.read_shp(path).to_undirected()
    new_nodes_co = set(G.nodes())
    for x in G.nodes():
        dists = [coord_dist(x,y) for y in set(cur_pos.values())]
        if min(dists)<1: new_nodes_co.remove(x)
    new_pos = {k:v for k,v in enumerate(new_nodes_co,max(cur_pos.keys())+1)}
    return_pos = cur_pos.copy()
    return_pos.update(new_pos)
    X = nx.Graph()
    X.add_nodes_from(return_pos.keys())
    co_edges = [set(x) for x in G.edges()]
    int_edges = [tuple(x for x,y in return_pos.items() if y in edg) 
                 for edg in co_edges]
    X.add_edges_from(int_edges)
    
    return X, return_pos

def two_links(poly,lines):
    """
    Function that takes in a polygon geometry and a line gdf and returns 
    information about the two closest points from the polygon to anywhere on 
    the lines.

    Parameters
    ----------
    poly : shapely polygon
        Polygon from which we are determining the nearest two points.
    lines : geopandas geodataframe
        Linestrings to which we are determining the nearest two points.

    Returns
    -------
    dist : list of floats
        List with two entries: the distance from the polygon to the nearest 
        line; the distance from the polygon to the second nearest line.
    link : list of tuples of tuples
        List with two entries: a tuple with the coordinates of the line 
        connecting the polygon to the nearest line; a tuple with the 
        coordinates of the line connecting the polygon to the second nearest 
        line.
    line_ind : list of ints
        List with two entries: the index of the closest line; the index of the 
        second closest line.
    """
    links = [None]*len(lines)
    for j, geom in enumerate(lines['geometry']):
        pt = nearest_points(poly,geom)
        links[j] = ((pt[0].x,pt[0].y),(pt[1].x,pt[1].y))
    zipped_list = zip([coord_dist(x,y) for x,y in links], 
                      links, list(lines.index))
    sorted_list = sorted(zipped_list)
    dist = [x for x,_,_ in sorted_list][0:2]
    link = [x for _,x,_ in sorted_list][0:2]
    line_ind = [x for _,_,x in sorted_list][0:2]
    return dist, link, line_ind 

def break_line(line,points):
    """
    Splits a line at specified points

    Parameters
    ----------
    line : TYPE
        DESCRIPTION.
    points : TYPE
        DESCRIPTION.

    Returns
    -------
    line_segs : TYPE
        DESCRIPTION.

    """
    dx = line.coords[0][0]-line.coords[-1][0]
    dy = line.coords[0][1]-line.coords[-1][1]
    pdx = 1/(math.sqrt(dx**2+dy**2))*dy
    pdy = -1/(math.sqrt(dx**2+dy**2))*dx
    
    lines = [None]*len(points)
    for i,j in enumerate(points):
        pt1 = Point(j.coords[0][0]+pdx*0.01*line.length, 
                    j.coords[0][1]+pdy*0.01*line.length)
        pt2 = Point(j.coords[0][0]-pdx*0.01*line.length,
                    j.coords[0][1]-pdy*0.01*line.length)
        lines[i] = LineString((pt1,pt2))
        
    break_lines = gpd.GeoSeries(lines)
    line_segs = [x for x in split(line,break_lines.unary_union)]
    
    return line_segs

def get_nodes(s,buildings_on_street,max_node_spacing):
    """
    Creates a GeoDataFrame of wastewater nodes for a street segment.

    Parameters
    ----------
    s : LineString
        Street segment.
    buildings_on_street : GeoDataFrame
        Buildings attached to the street segment.
    max_node_spacing : int
        Maximum spacing of nodes.

    Returns
    -------
    nodes : GeoDataFrame
        Wastewater nodes for the street segment, with column 'building_inds' to
        indicate which buildings are attached to each node.

    """
    # Nodes
    n_nodes = max(math.ceil(s.length/max_node_spacing)-1,1)
    nodes = gpd.GeoDataFrame({'geometry': 
                              [s.interpolate((2*(j+1)-1)/(2*n_nodes),
                                             normalized=True) 
                               for j in range(n_nodes)]
                              })

    # Buildings attached to each node
    buildings_at_node = [[] for k in range(n_nodes)]
    for geom,b_ind in zip(buildings_on_street["geometry"],
                          list(buildings_on_street.index)):
        _,to_node = nearest_points(geom,nodes["geometry"].unary_union)
        node_ind = nodes.loc[nodes["geometry"]==to_node].index[0]
        buildings_at_node[node_ind].append(b_ind)
    nodes['building_inds'] = buildings_at_node

    return nodes

def identify_ind(pt,pos):
    """
    Identifies the graph index corresponding to a point geometry.

    Parameters
    ----------
    pt : point
        Point geometry of the target node.
    pos : dicitonary
        Position dictionary of the graph.

    Returns
    -------
    ind : int
        Index of the target node.
    """
    dists = {i:coord_dist(pt.coords[0],xy) for i,xy in pos.items()}
    ind = [i for i,dist in dists.items() if dist==min(dists.values())][0]
    return ind

def sample_line(l,ras):
    """
    Funciton to get the profile of a raster along a line, with the number of 
    sample points based on cell size and line length.

    Parameters
    ----------
    l : LineString geometry (or GeoSeries)
        Line along which raster should be sampled.
    ras : rasterio.io.DatasetReader
        Raster to sample.

    Returns
    -------
    profile : List
        Values along the raster.

    """
    if isinstance(l,gpd.geoseries.GeoSeries):
        min_size = min((abs(ras.transform[0]),abs(ras.transform[4]))) 
        npts = max(math.ceil(2*l.length[l.index[0]]/min_size), 2) 
        samplepts = [l.interpolate(j/(npts-1), normalized=True) 
                     for j in range(npts)]
        profile = [[i for i in ras.sample([(j[l.index[0]].x, 
                                            j[l.index[0]].y)])][0][0] 
                   for j in samplepts]
   
    if isinstance(l,shapely.geometry.linestring.LineString):
        min_size = min((abs(ras.transform[0]),abs(ras.transform[4])))
        npts = max(math.ceil(2*l.length/min_size), 2)
        samplepts = [l.interpolate(j/(npts-1), normalized=True) 
                     for j in range(npts)]
        profile = [[i for i in ras.sample([(j.x, j.y)])][0][0] 
                   for j in samplepts]

    return profile

def max_min(val, max_val, min_val):
    """
    Returns value scaled using min-max scaling.

    Parameters
    ----------
    val : float
        Value to be scaled.
    max_val : float
        Maximum value in data (i.e., max_min(max_val) = 1).
    min_val : float
        Minimum value in data (i.e., max_min(min_val) = 0).

    Returns
    -------
    float
        Scaled value, between 0 and 1.

    """
    return (val-min_val)/(max_val-min_val)

def tree_topology(termnode):
    """
    Opens the cleaned streets at 'splitgdf/splitgdf.shp' as a graph, evaluates 
    weight (either using hybrid or distance method), and determines the 
    shortest path spanning tree).

    Parameters
    ----------
    termnode : GeoDataFrame
        GDF with a single point geometry: the location of the terminal node.

    Returns
    -------
    spst : networkx digraph
        Shortest path spanning tree.
    pos : dictionary
        Position dictionary for the output graph.
    tn_id : int
        Node id for the termnode.

    """
    if pr.weight_method == 'hybrid-slope':
        
        # Open splitgdf as an directed graph
        street_graph, pos = \
            digraph_from_shp('./inter/streets_clean/streets_clean.shp')
        
        # Evaluate weight of edges
        with rasterio.open('./input/dem_ft/dem_ft.tif') as dem:
            evaluated = set()
            for i,j in street_graph.edges():
                if (i,j) not in evaluated:
                    lij = coord_dist(pos[i],pos[j])
                    
                    hi = next(dem.sample([pos[i]]))[0]
                    hj = next(dem.sample([pos[j]]))[0]
                    
                    sij = (hi-hj)/lij
                    sji = (hj-hi)/lij
        
                    ftij = 1 - norm.cdf(sij, loc=pr.s_mu, scale=pr.s_sig)
                    ftji = 1 - norm.cdf(sji, loc=pr.s_mu, scale=pr.s_sig)
                    
                    street_graph[i][j]['weight'] = lij*ftij
                    street_graph[j][i]['weight'] = lij*ftji
                    evaluated.add((i,j))
                    evaluated.add((j,i))
            
        # Identify the index of the target node
        tn_id = identify_ind(termnode["geometry"][0],pos)
        
        # Shortest path spanning tree
        spst = get_spst(street_graph,tn_id,keepdists=True)            
    
    elif pr.weight_method == 'hybrid-elev':
        # Open splitgdf as a directed graph
        street_graph, pos = \
            digraph_from_shp('./inter/streets_clean/streets_clean.shp')
        
        # Sample elevation at every node
        elevs = dict.fromkeys(street_graph.nodes())
        with rasterio.open('./input/dem_ft/dem_ft.tif') as dem:
            for i in street_graph.nodes():
                elevs[i] = next(dem.sample([pos[i]]))[0]
                
        # Min-Max scaling of elevations
        elevs_scaled = dict.fromkeys(street_graph.nodes())
        max_elev = max(elevs.values())
        min_elev = min(elevs.values())
        for i in street_graph.nodes():
            elevs_scaled[i] = max_min(elevs[i], max_elev, min_elev)
            
        # Evaluate weight of each edge
        for i,j in street_graph.edges():
            lij = coord_dist(pos[i],pos[j])
            street_graph[i][j]['weight'] = lij*elevs_scaled[j]
            
        # Identify the index of the target node
        tn_id = identify_ind(termnode["geometry"][0],pos)
        
        # Shortest path spanning tree
        spst = get_spst(street_graph,tn_id,keepdists=True)   
            
    elif pr.weight_method == 'distance':
        # Open splitgdf as an undirected graph
        street_graph, pos = \
            graph_from_shp('./inter/streets_clean/streets_clean.shp')
        
        # Evaluate the weight as the distance between endpoints
        for i,j in street_graph.edges():
            street_graph[i][j]['weight'] = coord_dist(pos[i], pos[j])
            
        # Identify the index of the target node
        tn_id = identify_ind(termnode["geometry"][0],pos)
        
        # Shortest path spanning tree
        spst = get_spst(street_graph,tn_id,keepdists=True)
        
    return spst, pos, tn_id

def graphs_from_t(tree, tree_pos):
    """
    Function that takes the shortest path tree of the street network and 
    returns the equivalent tree for the split tree network, as well as the 
    associated design tree. 

    Parameters
    ----------
    tree : NetworkX DiGraph
        Shortest path tree of the street network.
    tree_pos : dictionary
        Position dictionary of tree.

    Returns
    -------
    split_pos : dictionary
        Position dictionary for the full_spt (note that this dictionary is 
        consistent with the input dictionary, just has extra key:value pairs 
        for the new nodes).
    full_spt : NetworkX DiGraph
        Full shortest path tree corresponding to the input tree (i.e., has all 
        of the nodes, including the new nodes).
    des_graph : NetworkX DiGraph
        Graph used for hydraulic design with no intermediate nodes, but 
        additional edges as compared to the input tree; includes 'path' 
        attribute for the equivalent path in full_spt.

    """
    # Open splitgdf as a graph, maintaining any existing node numbering
    split_graph, split_pos = \
        graph_from_shp_2('./inter/splitgdf/splitgdf.shp', tree_pos)

    # List of lists that fully defines the output graphs. Each list in the list 
    # represents an edge in the design graph. The list for each edge in the 
    # design_graph is the equivalent path in the full wastewater shortest path 
    # tree. 
    paths = []

    # For each edge in the tree, identify corresponding path in split_graph
    for i,j in tree.edges(): 
        for node_1 in split_graph.neighbors(i):
            path = [i,node_1]
            while path[-1] not in tree.nodes():
                next_node = [x for x in split_graph.neighbors(path[-1]) 
                             if x!=path[-2]][0]
                path.append(next_node)
            if path[-1]==j: 
                paths.append(path)
                break
    
    # Add paths for wastewater nodes that are not already in the tree
    nodes_included = set([x for p in paths for x in p])
    for i in split_graph.nodes():
        if i not in nodes_included:
            adj_nodes = [x for x in split_graph.neighbors(i)]
            path = [adj_nodes[0], i, adj_nodes[1]]
            while path[0] not in tree.nodes():
                next_node = [x for x in split_graph.neighbors(path[0]) 
                             if x!=path[1]][0]
                path.insert(0,next_node)
            while path[-1] not in tree.nodes():
                next_node = [x for x in split_graph.neighbors(path[-1]) 
                             if x!=path[-2]][0]
                path.append(next_node)
            w_i = tree.nodes[path[0]]['w_to_target']
            w_j = tree.nodes[path[-1]]['w_to_target']
            
            if w_i > w_j:
                path = path[1:]
            else:
                path.reverse()
                path = path[1:]
                
            paths.append(path)
            for x in path: nodes_included.add(x)
            
    # Prepare full shortest path tree for output   
    full_edges = []
    for path in paths:
        for i in range(len(path)-1):
            full_edges.append((path[i],path[i+1]))
    full_spt = nx.DiGraph()
    full_spt.add_nodes_from(split_graph.nodes())
    full_spt.add_edges_from(full_edges)
    
    # Prepare design graph for output
    des_edges = []
    des_nodes = set()
    for path in paths:
        des_edges.append((path[0],path[-1]))
        des_nodes.add(path[0])
        des_nodes.add(path[-1])
    des_graph = nx.DiGraph()
    des_graph.add_nodes_from(des_nodes)
    des_graph.add_edges_from(des_edges)
    
    for i,j in des_edges:
        path = [x for x in paths if x[0]==i and x[-1]==j][0]
        des_graph[i][j]['path'] = path

    return split_pos, full_spt, des_graph

def voronis(points_gs,area):
    """
    Function to get Voronoi areas for a set of points over an area

    Parameters
    ----------
    points_gs : GeoSeries
        Set of points, each to be associated with a Voronoi polygon.
    area : Polygon
        Total area covered by Voronoi polygons.

    Returns
    -------
    GeoSeries
        GeoSeries of Voronoi polygons, in same order as point_gs.

    """
    if len(points_gs)>2:
        polys, points = voronoi_regions_from_coords(points_gs,area)
        geom_list = [None]*len(points_gs)
        for x,y in polys.items():
            geom_list[points[x][0]] = y
            
    elif len(points_gs)==2:
        pt0 = points_gs[points_gs.index[0]]
        pt1 = points_gs[points_gs.index[1]]
        new_i = max(points_gs.index[0],points_gs.index[1])+1
        
        
        if pt0.distance(area.centroid)<pt1.distance(area.centroid):
            l = LineString((pt0,area.centroid))
            pt2 = l.interpolate(0.001,normalized = True)
            points_gs = points_gs.append(gpd.GeoSeries(pt2).set_axis([new_i]))
                
            polys, points = voronoi_regions_from_coords(points_gs,area)
            geom_list_3 = [None]*len(points_gs)
            for x,y in polys.items():
                geom_list_3[points[x][0]] = y
                
            geom_list = [None]*2
            geom_list[0] = geom_list_3[0].union(geom_list_3[2])
            geom_list[1] = geom_list_3[1]
            
        else:
            l = LineString((pt1,area.centroid))
            pt2 = l.interpolate(0.001,normalized = True)
            points_gs = points_gs.append(gpd.GeoSeries(pt2).set_axis([new_i]))
                
            polys, points = voronoi_regions_from_coords(points_gs,area)
            geom_list_3 = [None]*len(points_gs)
            for x,y in polys.items():
                geom_list_3[points[x][0]] = y
                
            geom_list = [None]*2
            geom_list[0] = geom_list_3[0]
            geom_list[1] = geom_list_3[1].union(geom_list_3[2])
    
    # Check that there are no multi-polygons, reassign any disconnected areas 
    for i,vor in enumerate(geom_list):
        if isinstance(vor, MultiPolygon):
            # List of areas to reassign
            reassign = []
            for vor_indiv in vor:
                if vor_indiv.intersects(points_gs.iloc[points_gs.index[i]]):
                    actual = vor_indiv
                else:
                    reassign.append(vor_indiv)
            
            # Change geometry to only the area that intersects the point
            try: geom_list[i] = actual
            # If there is some error in the voronoi polygons, assign a buffer 
            # instead
            except UnboundLocalError: 
                geom_list[i] = points_gs.iloc[points_gs.index[i]].buffer(100)
            
            # Merge areas to reassign with a random adjacent area
            for re in reassign:
                adj_inds = [j for j,vor2 in enumerate(geom_list) 
                            if re.intersects(vor2)]
                geom_list[adj_inds[0]] = re.union(geom_list[adj_inds[0]])
            
    return gpd.GeoSeries(geom_list) 

def area_column(x):
    """
    Calculates area from the Voronoi polygon associated with each storwater 
    node.

    Parameters
    ----------
    x : Row
        Row of the GeoDataFrame, from .apply(axis=1).

    Returns
    -------
    Float
        Area of the Voronoi polygon.

    """
    return x.voroni_areas.area

def width_column(x):
    """
    Calculates width of Voronoi area associated with each stormwater node.

    Parameters
    ----------
    x : Row
        Row of the GeoDataFrame, from .apply(axis=1).

    Returns
    -------
    Float
        Width of the Voronoi area.

    """
    ext_coords = zip(x.voroni_areas.exterior.xy[0],
                     x.voroni_areas.exterior.xy[1])
    max_dist = max([coord_dist(x.geometry.coords[0],(k,c)) 
                    for k,c in ext_coords])
    return x.area/max_dist

def slope_column(geom,ar,af):
    """
    Calculates the spatial average (median) of the slope raster over a polygon.
    Note: subsequent versions should put the units outside of this function.

    Parameters
    ----------
    geom : Polygon
        Polygon over which the raster is averaged.
    ar : Array
        Array of the slope raster. UNITS?
    af : Affine
        Affine transformation of the slope raster.

    Returns
    -------
    val : Float
        Spatial mediam of the slope raster over the area. UNITS?

    """
    slope_stat = zonal_stats(geom,ar,affine=af,stats=['median'],nodata=-9999)
    val = math.tan(slope_stat[0]['median']*math.pi/180)
    if val == 0: val=0.01 #Reasonable value to prevent ZeroDivisionError
    return val

def imp_column(geom,ar,af):
    """
    Calculates the spatial average (mean) of a raster over a polygon

    Parameters
    ----------
    geom : Polygon
        Polygon over which the raster is averaged.
    ar : Array
        Array of the raster.
    af : Affine
        Affine transformation of the raster.

    Returns
    -------
    Float
        Spatial mean of the raster over the area.

    """
    imp_stat = zonal_stats(geom, ar, affine=af, stats=['mean'], nodata=-999)
    # If no cells are sampled, re-sample with a buffer
    if imp_stat[0]['mean'] == None:
        imp_stat = zonal_stats(geom.buffer(50), ar, affine=af, stats=['mean'], 
                               nodata=-999)
    return imp_stat[0]['mean']

def runoffc10(row,HSG):
    """
    Caluclates the runoff coefficient for a 10-year storm.

    Parameters
    ----------
    row : Row
        Row of the GeoDataFrame, from .apply(axis=1).
    HSG : String
        Hydrologic soil group.

    Returns
    -------
    c5 : Float
        Runoff coefficient for the Voroni area.

    """
    if HSG=='A':
        c10 = 0.87*(row.percent_impervious/100)**1.232
    elif HSG=='B':
        c10 = 0.81*(row.percent_impervious/100)+0.057
    elif HSG=='C' or HSG=='D':
        c10 = 0.74*(row.percent_impervious/100)+0.132
    return c10

def cia(row):
    """
    Calculates nodal flow from rational method.

    Parameters
    ----------
    row : Row
        Row of the GeoDataFrame, from .apply(axis=1).

    Returns
    -------
    Float
        Nodal flow value from rational method (ft3/s).
        
    """
    c = row.c10
    i = pr.storm_intensity
    a = row.area
    return c*i*a/(12*3600)

def accumulate(t,data):
    """
    Accumulates a value over a tree, where data is stored in the nodes of the 
    tree, and accumulated data is stored as data+'_acc', also in the nodes.

    Parameters
    ----------
    t : networkx DiGraph
        Directed graph (must be a tree) over which the values stored in the 
        nodes will be accumulated.
    data : string
        String naming the data in the nodes of the graph.

    Returns
    -------
    networkx DiGraph
        Modified graph with the accumulated data stored in the nodes as 
        data+'_acc'.

    """
    # Check that input graph is a tree, print error if not
    if not nx.is_tree(t): 
        print('ERROR: graph passed to accumulate is not a tree')
        
    # Internal recursive function, continue until reach leaves
    def recurse_sum(node):
        if acc_val[node]==None:
            node_result = val[node]
            for j in t.predecessors(node):
                node_result = node_result+recurse_sum(j)
            acc_val[node] = node_result
            return node_result
        else:
            return acc_val[node]    
    
    # Dictionaries to store un-accumulated and accumulated values
    val = dict.fromkeys(t.nodes())
    acc_val = dict.fromkeys(t.nodes())
    
    for i in t.nodes():
        # Populate un-accumulated value for every node
        try: 
            val[i] = t.nodes[i][data] 
        # Zero if there is no data at a node
        except KeyError: 
            val[i] = 0
        # Transfer value to accumulated value for leaf nodes    
        if t.in_degree(i)==0: #
            acc_val[i] = val[i]
        # Identify root node
        if t.out_degree(i)==0:
            root = i
    
    # Recursive function to populate acc_val 
    _ = recurse_sum(root)
    
    # Store results in the tree and return
    for i in t.nodes():
        t.nodes[i][data+'_acc'] = acc_val[i]
    return t

def peak_f(pop):
    """
    Evaluates peak factor for wastewater inflow based on contributory 
    population.

    Parameters
    ----------
    pop : int/float
        Population.

    Returns
    -------
    f : float
        Peak factor.

    """
    f = (18+math.sqrt(pop/1000))/(4+math.sqrt(pop/1000))
    return f

def ww_pop(row):
    """
    Evaluates wastewater design population, using either building-based or 
    population-based method.

    Parameters
    ----------
    row : Row
        Row of the buildings GeoDataFrame, from .apply(axis=1).

    Returns
    -------
    pop : int/float
        Design population for each building.

    """
    if pr.ww_method == 'buildings':
        if row.type=='residential':
            if row.nunits==1:
                pop = pr.design_pop_single
            else:
                pop = pr.design_pop_multi*row.nunits
        else:
            pop = 0
    elif pr.ww_method == 'population':
        pop = row.population*pr.factor_designpop
    return pop

def ww_ci(row):
    """
    Evaluates and sums the commercial and industrial wastewater design flow for
    each building, using either the building-based or populaiton-based method.

    Parameters
    ----------
    row : Row
        Row of the buildings GeoDataFrame, from .apply(axis=1).

    Returns
    -------
    qci : float
        Commercial design flow + industrial design flow.

    """
    if pr.ww_method == 'buildings':
        if row.type == 'commercial':
            qci = row.area_ft2*pr.ww_comm
        elif row.type == 'industrial':
            qci = row.area_ft2*pr.ww_ind
        elif row.type == 'hotel':
            qci = row.nunits*pr.ww_hotel
        else:
            qci = 0
    elif pr.ww_method == 'population':
        pop_ci = (pr.factor_ci-1)*pr.factor_designpop*row.population
        qci = pop_ci*pr.ww_percap
    return qci

def remove_intermediate_nodes(nodes,t):
    """
    Creates a new directed tree that removes intermediate nodes by only 
    including the specified nodes, maintaining the connection through 
    intermediate nodes, and also returns the full path on the original graph as
    an attribute of the edges in the reduced graph.    

    Parameters
    ----------
    nodes : set
        Set of nodes to be included in the new graph.
    t : DiGraph
        Original graph with intermediate nodes.

    Returns
    -------
    G : DiGraph
        Graph without intermediate nodes, including the corresponding path in 
        the original graph as a list in the edge data.

    """
    # Create edge dictionary with keys = edge tuples and values = full paths
    edge_dict = dict()
    for i in nodes:
        if t.out_degree(i)==0:
            pass
        else:
            path = [i]
            j = [y for x,y in t.edges() if x==i][0]
            path.append(j)
            while j not in nodes:
                j = [y for x,y in t.edges() if x==j][0]
                path.append(j)
            edge_dict[(i,j)] = path
    
    # Create new graph
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edge_dict.keys())
    nx.set_edge_attributes(G,edge_dict,'path')
    return G

def depths(t):
    """
    DFS to determine edge and node depths.

    Parameters
    ----------
    t : DiGraph
        Directed tree with elevation attribute in nodes; and length and slope 
        attributes in edges.

    Returns
    -------
    DiGraph
        Returns the same graph, with nodal depths specified as node attribute 
        "depth", and elevations of two ends of each edge as edge attributes 
        "depth0" and "depth1".

    """
    if not nx.is_tree(t): print('ERROR: graph passed to depths is not a tree')
    
    # Internal recursive function, continue until reach leaves
    def recurse_depth(node):
        if node_depth[node]==None:
            node_depth[node] = max([recurse_depth(i) + depth_increase[(i,node)]
                                    for i in t.predecessors(node)])
            return node_depth[node]
        else:
            return node_depth[node]
    
    node_depth = dict.fromkeys(t.nodes())
    node_elev = dict.fromkeys(t.nodes())
    
    for i in t.nodes():
        node_elev[i] = t.nodes[i]['elev']
        if t.in_degree(i) == 0:
            node_depth[i] = pr.min_depth
        if t.out_degree(i) == 0:
            root = i
    
    depth_increase = dict.fromkeys(t.edges())
   
    for i,j in t.edges():
        l = t.edges[(i,j)]['length']
        s = t.edges[(i,j)]['slope']
        ground_slope = (node_elev[i]-node_elev[j])/l
        
        if ground_slope>s:
            s = ground_slope
            t.edges[(i,j)]['slope'] = ground_slope
  
        depth_increase[(i,j)] = (s-ground_slope)*l
            
    # Call to recursive function to get depth at root (and all other nodes)
    _ = recurse_depth(root)
    
    # Store the results in the tree and return the modified tree
    for i in t.nodes():
        t.nodes[i]['depth'] = node_depth[i]
    for i,j in t.edges():
        t.edges[(i,j)]['depth0'] = node_depth[i]
        t.edges[(i,j)]['depth1'] = node_depth[i]+depth_increase[(i,j)]
        
    return t