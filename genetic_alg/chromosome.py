import scipy as sp

sp_rand = sp.random.random
sp_randint = sp.random.randint

class Chromosome(sp.ndarray):
    @classmethod
    def setup(cls, fitness, encoding):
        cls.fitness = fitness
        cls.encoding = encoding
        cls.length = encoding.length

    def __new__(cls, birth=1, iterable=None):
        self.birth = birth
        if iterable is not None:
            arr = sp.concatenate([part[:] for part in iterable])
            if arr.length != cls.length:
                raise TypeError("longitud del cromosoma: %d".format(arr.length))
            return arr
        else:
            return sp_rand((cls.length,))

    def __lt__(self, other):
        return self.score < other.score
    def __le__(self, other):
        return self.score <= other.score
    def __gt__(self, other):
        return self.score > other.score
    def __ge__(self, other):
        return self.score >= other.score

    def __getattr__(self, name):
        if name == 'score':
            self.score = self.fitness()
            return self.score
        else:
            return super(Chromosome, self).__getattr__(self, name)

    def __array_finalize__(self, obj):
        if obj is None: return
        self.birth = getattr(obj, 'birth', None)

    def mutate(self):
        self[sp_randint(0, Chromosome.length)] = sp_rand()
        self.score = self.fitness()
