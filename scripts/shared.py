import time


def timeit(func, *args, repeats=10, **kwargs):
    total_time = 0
    for _ in range(repeats):
        start = time.time()
        func(*args, **kwargs)
        total_time += time.time() - start
    return total_time / repeats
