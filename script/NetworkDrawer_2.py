import networkx as nx
import matplotlib.pyplot as plt
import os
import os.path as path
import sys
sys.path.append(path.join(path.dirname( path.dirname( path.abspath(__file__))),'database'))
from connect import db
from openbabel import pybel

def extract_data():
    
    reaction_collection = db['reactions']
    qm_collection = db['qm_calculate_center']
    acceptable_condition = ['forward equal to reactant and reverse equal to product',
                            'reverse equal to reactant and forward equal to product',
                            'forward equal to reactant but reverse does not equal to product',
                            'reverse equal to reactant but forward does not equal to product']

    query = {'$and': 
                    [
                    {"unique":
                        {"$in":
                        ['new one']}},
                    {'irc_equal':
                        {'$in':acceptable_condition}}
                    ]
                }
    """
    q = {'irc_equal':{'$in':acceptable_condition}}
    a = list(reaction_collection.find(q))

    for i in a:
        try:
            print(i['reactant_inchi_key'].split()[0])
            pyMol_1 = xyz_to_pyMol(path.join(i['path'], 'irc_reactant.xyz'))
            print(pyMol_1.write('inchiKey').strip())
            print(i['reactant_smi'].split()[0])
            print(pyMol_1.write('can').split()[0])
            print(i['generations'])
            print(round(i['barrier_energy'], 2))
            print('-------')
        except:
            pass
    """
    targets = list(reaction_collection.find(query))
    reactant_smi = [target['reactant_smi'].split()[0] for target in targets]
    product_smi = [target['product_smi'].split()[0] for target in targets]
    barrier_energy = []
    for i in targets:
        ts = list(qm_collection.find({'path':i['path']}))[0]
        print((ts['product_xtb_hf'] - ts['reactant_xtb_hf'])*627.5095)
        #print(i['path'])
        #print(round(i['barrier_energy'], 2))
        barrier_energy.append(round(i['barrier_energy'], 2))
        
    generations = [target['generations'] for target in targets]
    zipped = zip(reactant_smi, product_smi, barrier_energy, generations)

    return zipped

def xyz_to_pyMol(xyz):
    mol = next(pybel.readfile('xyz', xyz))
    return mol

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
            with_labels=False,
            node_color='lightgreen',
            font_size = 8)
    
    plt.savefig("simple_path_2.png") # save as png
    plt.show()


draw()