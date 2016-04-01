import scipy as sp

random = sp.random.random
randint = sp.radom.randint

class Chromosome(sp.ndarray):
    @classmethod
    def setup(cls, fitness):
        cls.fitness = fitness

    def __new__(cls, birth=1):
        self.birth = birth
        return random((Chromosome.length,))

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
        self[randint(0, Chromosome.length)] = random()
        self.score = fitness()
