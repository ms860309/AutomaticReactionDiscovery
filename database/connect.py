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
query = {'$and': 
                [
                {'product_inchi_key':
                    {'$in':
                    ['VHWYCFISAQVCCP-UHFFFAOYSA-N']}
                    },
                {'generations':
                    {'$in':
                        [1]
                    }}
                ]
            }

targets = list(qm_collection.find(query))
print(len(targets))
for i in targets:
    print(i['ssm_status'])
"""