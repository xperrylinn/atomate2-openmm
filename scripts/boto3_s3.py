import boto3

session = boto3.Session(profile_name="atomate2-openmm-dev")
client = session.client("s3")

response = client.list_buckets()

buckets_by_name = {bucket["Name"]: bucket for bucket in response["Buckets"]}

assert "atomate2-openmm" in buckets_by_name
