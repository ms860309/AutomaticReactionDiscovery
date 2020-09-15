import networkx as nx
import matplotlib.pyplot as plt
import os
import os.path as path
import sys
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'database'))
from connect import db

def extract_data():
    
    reaction_collection = db['reactions']
    query = {"unique":
                {"$in":
                    ['new one']
                }}
    targets = list(reaction_collection.find(query))
    reactant_smi = [target['reactant_smi'].split()[0] for target in targets]
    product_smi = [target['product_smi'].split()[0] for target in targets]
    barrier_energy = []
    for i in targets:
        barrier_energy.append(round(i['barrier_energy'], 2))
        
    generations = [target['generations'] for target in targets]
    zipped = zip(reactant_smi, product_smi, barrier_energy, generations)

    return zipped

def draw():
    G=nx.Graph() # create object
    
    zipped = extract_data()
    _dict = {}
    for i, j, k, l in list(zipped):
        if l == 1:
            G.add_edge(i,j,color='r')
            _dict[(i,j)] = k
        elif l == 2:
            G.add_edge(i,j,color='g')
            _dict[(i,j)] = k
        else:
            G.add_edge(i,j,color='b')
            _dict[(i,j)] = k
            
    plt.figure()
    
    colors = nx.get_edge_attributes(G,'color').values()
    weights = nx.get_edge_attributes(G,'weight').values()

    #pos = nx.circular_layout(G)
    pos = nx.spring_layout(G)  # positions for all nodes
    nx.draw_networkx_edge_labels(G, pos, edge_labels=_dict, font_size=8)
    nx.draw(G, pos, 
            edge_color=colors, 
            with_labels=True,
            node_color='lightgreen',
            font_size = 8)
    
    plt.savefig("simple_path_2.png") # save as png
    plt.show()


draw()