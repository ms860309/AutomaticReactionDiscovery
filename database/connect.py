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
collect = db['reactions']
reg_query = {"ts_status":
                {"$in": 
                    ["job_success"] 
                }
            }
targets = list(collect.find(reg_query))

selected_targets = []
gens = []
equal = []
for target in targets:
    dir_path = target['path']
    selected_targets.append(dir_path)
    gen = target['next_gen_num']
    gens.append(gen)
    ard_ssm_equal = target['ard_ssm_equal']
    equal.append(ard_ssm_equal)
zipped = zip(selected_targets, gens, equal)
for target in list(zipped):
    dir_path, gen_num, ard_ssm_equal = target[0], target[1], target[2]
    print(dir_path)
"""