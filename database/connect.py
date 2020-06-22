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
        self.db = self.client['network2']

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

"""
