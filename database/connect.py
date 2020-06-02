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
collect = db['qm_calculate_center']
query = {'$and': 
                [
                { "ts_status":
                    {"$in":
                    ['job_success']}},
                {'energy_status':
                    {'$in':
                        ['job_success']}}
                ]
            }
targets = list(collect.find(query))
num =[]
print(len(targets))
for idx, i in enumerate(targets):
    print(i)
    try:
        barrier_energy = i['barrier_energy']
    except:
        num.append(idx)
print(num)
print(len(targets))
"""