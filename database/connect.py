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
        self.db = self.client['new_reactant']

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=2000)
        return client

client = Connector()
db = getattr(Connector(), 'db')



# debug
"""
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
