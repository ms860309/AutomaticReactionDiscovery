import pymongo
from pymongo import MongoClient
import shutil
import os
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



# debug
"""
qm_collection = db['qm_calculate_center']
reg_query = {"irc_forward_status":
                {"$in":
                    ["job_launched"]
                }
            }
targets = list(qm_collection.find(reg_query))
for target in targets:
    dir_path = target['path']
    irc_path = os.path.join(dir_path, 'IRC')
    if os.path.exists(irc_path):
        shutil.rmtree(irc_path)
    update_field = {'irc_status':"unrun"}
    qm_collection.update_one({'path':dir_path}, {"$unset": {'irc_forward_jobid':"", 'irc_forward_status':"", 'irc_reverse_jobid':"", 'irc_reverse_status':""}, "$set": update_field}, True)



qm_collection = db['reactions']
query = {'$and':
                [{"for_debug":
                {"$in":
                    ['new one']
                }
            }, {'unique':
                {'$nin':
                    ['reactant equal to product']}
            }]}
query_2 = {'ssm':'job_unrun'}
targets = list(qm_collection.find(query_2))
for i in targets:
    print(i['reaction'][1])
"""
qm_collection = db['qm_calculate_center']

targets = list(qm_collection.aggregate(
    [
        {"$match": {'reactant_inchi_key':'AOPVWXNGMQWVNU-UHFFFAOYSA-N'}},
        {"$group":{"_id":{}, 'reactant_scf_energy':{'$min':"$reactant_scf_energy"}}}
    ]
))
print(targets[0]['reactant_scf_energy'])
