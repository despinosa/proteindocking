def fibonacci(factor, start=1, **kwargs):
    prev_n = start
    n = start
    while True:
        yield int(factor * n)
        prev_n, n = n, n+prev_n

def linear(factor, start=1, step=1, **kwargs):
    n = start
    while True:
        yield int(factor * n)
        n += step

def polynomial(factor, start=1, step=1, power=1.5):
    n = start
    while True:
        yield int(factor * n ** power)
        n += step

def exponential(factor, start=0, step=1, base=2):
    n = start
    while True:
        yield int(factor * base ** n)
        n += step
