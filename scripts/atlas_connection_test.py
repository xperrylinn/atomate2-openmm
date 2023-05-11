from pymongo.mongo_client import MongoClient
import os


username, password = os.environ.get("ATLAS_USERNAME"), os.environ.get("ATLAS_PASSWORD")
if username is None or password is None:
    raise ValueError(f"Environment variables ATLAS_USERNAME and ATLAS_PASSWORD must be set.\n"
                     f"ATLAS_USERNAME - Atlas MongoDB username\n"
                     f"ATLAS_PASSWORD - database password\n")

# Setup store using Atlas MongoDB
uri = f"mongodb+srv://{username}:{password}@atomate2-openmm.vlzvqsg.mongodb.net/?retryWrites=true&w=majority"

# Create a new client and connect to the server
client = MongoClient(uri)

# Send a ping to confirm a successful connection
try:
    client.admin.command('ping')
    print("Pinged your deployment. You successfully connected to MongoDB!")
except Exception as e:
    print(e)