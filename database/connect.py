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
reaction_collection = db['reactions']
query = {'$and':
                [{"unique":
                {"$in":
                    ['new one']
                }
            }, {'for_debug':
                {'$nin':
                    ['from same']}
            }]}

targets = list(reaction_collection.find(query))
dic = {}
for idx, i in enumerate(targets):
    rxn = i['reaction']
    gen = i['generations']
    barrier = i['barrier_energy']
    name = 'reaction_{}'.format(idx)
    dic[name] = rxn

print(dic)
"""