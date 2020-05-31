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
qm_cal_center = db['qm_calculate_center']
query = {'$or': 
                [
                {'ts_status':
                    {'$in':
                    ['job_fail', 'job_success']}
                    },
                {'ssm_status':
                    {'$in':
                        ['Exiting early', 'total dissociation', 'job_fail']
                    }}
                ]
            }
targets = list(qm_cal_center.find(query))
print(targets)
print(len(targets))
"""