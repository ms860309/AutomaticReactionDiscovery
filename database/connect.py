import pymongo
from pymongo import MongoClient

class Connector(object):

    def __init__(self, server):
        #self.host = host
        #self.port = port
        self.server = server
        #self.mongo_db = mongo_db
        self.client = self.connect()

    def connect(self):
        client = MongoClient(self.server, serverSelectionTimeoutMS=2000)
        return client








