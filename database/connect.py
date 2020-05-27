import pymongo
from pymongo import MongoClient

class Connector(object):

    def __init__(self):
        #self.host = host
        #self.port = port
        self.server = 'mongodb+srv://jianyi:aa123@cluster0-wo5fn.gcp.mongodb.net/test?retryWrites=true&w=majority'
        #self.server = 'mongodb://localhost:27017/'
        #self.mongo_db = mongo_db
        self.client = self.connect()
        self.db = self.client['network']

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=2000)
        return client

client = Connector()
db = getattr(Connector(), 'db')


"""
co = db['reactions']
reg_query = {"ts_status":
                {"$in": 
                    ["job_success"] 
                }
            }
targets = list(co.find(reg_query))
print(targets)

selected_targets = []
gens = []
for target in targets:
    dir_path = target['path']
    selected_targets.append(dir_path)
    gen = target['next_gen_num']
    gens.append(gen)
zipped = zip(selected_targets, gens)

for i in targets:
    
    print(i['ard_ssm_equal'])

# debug
product_pool = db['pool']
mol_obj_inchi_key = 'QEDGWKGQ--wGQREQGN'
reg_query = {"reactant_inchi_key":
                {"$in":
                    [mol_obj_inchi_key]
                }
            }
targets = list(product_pool.find(reg_query))
if not targets:
    print('Nothing in')
    product_pool.insert_one({'reactant_inchi_key':mol_obj_inchi_key})
else:
    for i in targets:
        update_field = {'reactant_inchi_key':'123123'}
        product_pool.update_one(i, {"$set": update_field}, True)
"""