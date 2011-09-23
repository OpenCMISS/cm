"""Utility routines used by the OpenCMISS Python interface
"""

class CMISSError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def call_routine(routine,args=None):
    """Call a routine and check the return value, raise an
    exception if it is non-zero and return any other return
    values"""
    if args is None:
        r = routine()
    else:
        r = routine(*args)
    if isinstance(r,tuple):
        status = r[0]
        if len(r) == 1:
            return_val = None
        elif len(r) == 2:
            return_val = r[1]
        else:
            return_val = r[1:]
    else:
        status = r
        return_val = None
    if status != 0:
        #TODO: get error strings from OpenCMISS
        raise CMISSError, 'Non-zero return value'
    return return_val

