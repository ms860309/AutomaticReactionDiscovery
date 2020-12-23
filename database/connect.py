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
        self.db = self.client['xtb_irc']

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=2000)
        return client

client = Connector()
db = getattr(Connector(), 'db')



# debug

"""
reactions_collection = db['reactions']
reactions = list(reactions_collection.aggregate([{
                                        '$group':{
                                                '_id': "$reaction",
                                                'barrier': {'$min': "$barrier_energy"}}}
                                        ]))
for i in reactions:
    query = {'$and': 
                    [
                    { "reaction":i['_id']},
                    {'barrier_energy':i['barrier']}
                    ]
                }
    target = list(reactions_collection.find(query))[0]
    print(target['barrier_energy'])



qm_collection = db['qm_calculate_center']
query = {'low_opt_status':"job_success"}
a = list(qm_collection.find(query))
for i in a:
    print(i['low_energy'])

qm_collection = db['qm_calculate_center']
query = [{'$match':{'reactant_inchi_key':'OWCQMKVAAHGRRF-UHFFFAOYSA-N'}},
            {'$group':{'_id':'$reactant_inchi_key', 'reactant_mopac_hf':{'$min':'$reactant_mopac_hf'}}}]
a = list(qm_collection.aggregate(query))[0]['reactant_mopac_hf']
print(a)


qm_collection = db['qm_calculate_center']
query = {'barrier':'do not have reactant_energy'}
targets = list(qm_collection.find(query))
for target in targets:
    update_field = {'energy_status':"job_unrun"}
    qm_collection.update_one(target, {"$set": update_field}, True)
"""