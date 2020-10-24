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
        self.db = self.client['auto_bond_3']

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=2000)
        return client

client = Connector()
db = getattr(Connector(), 'db')



# debug
"""
qm_collection = db['qm_calculate_center']
query = {'reactant_inchi_key':'OWCQMKVAAHGRRF-UHFFFAOYSA-N'}
a = list(qm_collection.find(query))
for i in a:
    print(i['reactant_mopac_hf'])


qm_collection = db['qm_calculate_center']
query = [{'$match':{'reactant_inchi_key':'OWCQMKVAAHGRRF-UHFFFAOYSA-N'}},
            {'$group':{'_id':'$reactant_inchi_key', 'reactant_mopac_hf':{'$min':'$reactant_mopac_hf'}}}]
a = list(qm_collection.aggregate(query))[0]['reactant_mopac_hf']
print(a)


qm_collection = db['qm_calculate_center']
query = {'ssm_status':'job_fail'}
targets = list(qm_collection.find(query))
for target in targets:
    dir_path = target['path']
    ssm_path = os.path.join(dir_path, 'SSM')
    if os.path.exists(ssm_path):
        shutil.rmtree(ssm_path)
    update_field = {'ssm_status':"job_unrun"}
    qm_collection.update_one(target, {"$unset": {'ssm_jobid':""}, "$set": update_field}, True)
"""
