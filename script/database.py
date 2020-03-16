import pymongo

client = pymongo.MongoClient(host = 'localhost', port = 27017)
# define database
db = client.network
# define collection