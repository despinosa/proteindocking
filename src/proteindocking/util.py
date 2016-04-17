class Singleton(type):
    """docstring for Singleton"""
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]        
