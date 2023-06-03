from jobflow import job
from jobflow import Flow


@job
def add(a, b):
    return a + b

add_first = add(1, 5)

add_second = add(add_first.output, 3)

flow = Flow([add_first, add_second])

flow.draw_graph(figsize=(3, 3)).show()
