from datetime import datetime
from proteindocking.alpsimpl.dummyalps import ALPS

if __name__ == '__main__':
    a = ALPS()
    start = datetime.now()
    a.run()
    print datetime.now() - start
