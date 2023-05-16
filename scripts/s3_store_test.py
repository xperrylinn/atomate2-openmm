from maggma.stores.aws import S3Store
from maggma.stores import MemoryStore

index = MemoryStore(collection_name="index", key="blob_uuid")

index.connect()

s3_store = S3Store(
    index=index,
    bucket="atomate2-openmm",
    endpoint_url="https://s3.us-west-1.amazonaws.com",
    s3_profile="atomate2-openmm-dev",
    key="blob_uuid",
    s3_workers=1,
)

s3_store.connect()

s3_store.update(
    docs={"blob_uuid": "my_fancy_uid", "message": "hello world!"},
    key="blob_uuid",
)

print(list(s3_store.query()))
