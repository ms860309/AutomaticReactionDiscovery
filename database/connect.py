import pymongo

class Connector(object):

    def __init__(self, host, port, mongo_db, collection, data):
        self.host = host
        self.port = port
        self.mongo_db = mongo_db
        self.collection = collection
        self.data = data
        self.client = self.connect()

    def connect(self):
        client = pymongo.MongoClient(self.host, self.port)
        db = client[self.mongo_db]
        collection = db[self.collection]
        collection.insert_one(self.data)
        print(collection)


network = {'reaction2': ['00000', '00003'], 'reaction3': ['00000', '00004'], 'reaction0': ['00000', '00001'], 'reaction1': ['00000', '00002'], 'reaction6': ['00000', '00003', '00005', '00007'], 'reaction4': ['00000', '00003', '00005'], 'reaction5': ['00000', '00003', '00005', '00006']}
Connector(host = 'localhost', port = 27017, mongo_db = 'network', collection = 'reaction',data = network).connect()











