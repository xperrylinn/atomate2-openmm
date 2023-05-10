from maggma.stores import MemoryStore
from maggma.stores.file_store import FileStore
from jobflow import JobStore, job, run_locally
from pathlib import Path

file_path = Path("data/hello_world.txt")


@job(file_store="file")
def some_task():
    """
    The file_store argument in the @job decorator tells jobflow to store the value associated
    with the file_store in the additional stores provided to a JobStore. By default other outputs
    stored into the doc_store.
    """
    return {
        "file": file_path,
        "doc_store": {"hello world": "123"}
    }


some_task_job = some_task()

docs_store = MemoryStore()

file_store = FileStore(path="./data/", read_only=False)

store = JobStore(docs_store, additional_stores={"file_store": file_store})

output = run_locally(some_task_job, store=store, ensure_success=True)

assert str(file_path) == next(file_store.query())["path"]
