def fibonacci(start=1, **kwargs):
    prev_n = start
    n = start
    while True:
        yield n
        prev_n, n = n, n+prev_n

def linear(start=1, step=1, **kwargs):
    n = start
    while True:
        yield n
        n += step

def polynomial(start=1, step=1, power=2):
    n = start
    while True:
        yield n ** power
        n += step

def exponential(start=0, step=1, base=2):
    n = start
    while True:
        yield base ** n
        n += step
