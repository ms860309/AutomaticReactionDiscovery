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
max_gen = qm_collection.find_one(sort=[("generations", -1)])
max_gen = max_gen['generations']
query = {'ard_status':'job_unrun', 'generations':max_gen}
targets = list(qm_collection.find(query))
print(len(targets))
print(targets)
raise
ssm_query = {'$and':
                [{"ssm_status":
                {"$in":
                    ['job_success', 'job_running']
                }
            },
            {
                'generations':1
            }]}
ssm_success_number = len(list(qm_collection.find(ssm_query)))
ts_query = {'$and':
            [{'$or': 
                [
                {'ts_status':
                    {"$in":
                        ['job_success', 'job_fail']}},
                {'energy_status':
                    {'$in':
                        ['job_success', 'job_fail']}}
                ]
            },
            {
                'generations':1
            }]}
ts_fail_and_success_number = len(list(qm_collection.find(ts_query)))
print(ssm_success_number)
print(ts_fail_and_success_number)
"""