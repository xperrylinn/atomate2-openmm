from maggma.stores import MemoryStore
from jobflow import JobStore, job, run_locally


additional_data = {"atomate2-openmm": "awesome", "value": 1}
data = {"hello world": "123abc"}


@job(additional_mem_store="my_additional_data")
def some_task():
    """
    The file_store argument in the @job decorator tells jobflow to store the value associated
    with the file_store in the additional stores provided to a JobStore. By default other outputs
    stored into the doc_store.
    """
    return {
        "my_additional_data": additional_data,
        "my_data": data
    }


some_task_job = some_task()

docs_store = MemoryStore()

additional_mem_store = MemoryStore()

store = JobStore(docs_store, additional_stores={"additional_mem_store": additional_mem_store})

output = run_locally(some_task_job, store=store, ensure_success=True)

# Assert document in doc store is linked with data in additional store
assert next(docs_store.query())["output"]["my_additional_data"]["blob_uuid"] == next(additional_mem_store.query())["blob_uuid"]

# Assert the data in the additional store matches correct what we put in
assert additional_data == next(additional_mem_store.query())["data"]

# Assert the data in the doc store matches what we put in
assert data == next(docs_store.query())["output"]["my_data"]
