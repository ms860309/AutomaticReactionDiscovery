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



# debug
"""
qm_collection = db['qm_calculate_center']
query = {'irc_forward_status':'job_launched'}
targets = list(qm_collection.find(query))
for target in targets:
    import os
    IRC_dir_path = os.path.join(target['path'], 'IRC/')
    import shutil
    shutil.rmtree(IRC_dir_path)
    qm_collection.update_one({"path":target['path']}, {"$unset": {'irc_forward_status':"", 'irc_forward_jobid':"", 'irc_reverse_jobid':"", 'irc_reverse_status':""}, "$set": {'irc_status':job_unrun}}, True)
"""