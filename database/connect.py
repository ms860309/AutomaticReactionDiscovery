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
targets = list(qm_collection.aggregate(
    [
        {"$match": {"$and":[
            {
                "ts_status": {"$in":['job_success']},
                "energy_status": {"$in":['job_success']}
            }
        ]}},
        {"$group":{"_id":"$reaction", 'barrier_energy':{'$min':"$barrier_energy"}}}
    ]
))
for i in targets:
    reaction = i['_id']
    barrier = i['barrier_energy']
    query = {
        '$and':[
            {
                'reaction':reaction    
            },
            {
                'barrier_energy':barrier
            }
        ]
    }

query = {'ard_status':'job_unrun'}
targets = list(qm_collection.find({},{'barrier_energy':1}))
for i in targets:
    print(i)
"""