import boto3


def clear_bucket(bucket_name):
    session = boto3.Session(profile_name="atomate2-openmm-dev")
    s3_client = session.client('s3')

    # List all objects in the bucket
    response = s3_client.list_objects_v2(Bucket=bucket_name)
    objects = [{'Key': obj['Key']} for obj in response.get('Contents', [])]

    # Delete objects in batches of 1000
    while objects:
        delete_keys = {'Objects': objects[:1000]}
        s3_client.delete_objects(Bucket=bucket_name, Delete=delete_keys)
        objects = objects[1000:]

# Example usage
bucket_name = "atomate2-openmm"
clear_bucket(bucket_name)
