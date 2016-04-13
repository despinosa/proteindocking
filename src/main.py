from datetime import datetime
from proteindocking.alpsimpl.dummyalps import ALPS

if __name__ == '__main__':
    with open('dump.txt', 'w') as f:
        for i in xrange(99):
            a = ALPS()
            start = datetime.now()
            a.run()
            f.write('time: {}\n'.format(datetime.now()-start))
            f.write('best: {}\n\n'.format(a.layers[-1].population[0]))
