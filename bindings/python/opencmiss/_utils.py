"""Utility routines and classes used by OpenCMISS
"""

import _opencmiss_swig

class CMISSError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class CMISSType(object):
    """Base class for all OpenCMISS types
    """

    def __getattr__(self,attr_name):
        get_method = attr_name[0].upper()+attr_name[1:]+'Get'
        try:
            method = self.__getattribute__(get_method)
            return method()
        except AttributeError:
            raise AttributeError("CMISS type has no %s routine" % (get_method))

    def __setattr__(self,attr_name,value):
        set_method = attr_name[0].upper()+attr_name[1:]+'Set'
        try:
            method = self.__getattribute__(set_method)
            method(value)
        except AttributeError:
            self.__dict__[attr_name] = value


def wrap_cmiss_routine(routine, args=None):
    """Call a routine and check the return value, raise an
    exception if it is non-zero and return any other return values
    """

    if args is None:
        r = routine()
    else:
        #Replace wrapped cmiss types with the underlying type
        new_args = []
        for arg in args:
            if hasattr(arg,'cmiss_type'):
                new_args.append(arg.cmiss_type)
            else:
                try:
                    new_args.append([a.cmiss_type for a in arg])
                except (TypeError, AttributeError):
                    new_args.append(arg)
        r = routine(*new_args)
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
    if status != _opencmiss_swig.cvar.CMISSNoError:
        if status == _opencmiss_swig.cvar.CMISSPointerIsNULL:
            raise CMISSError, "CMISS type pointer is null"
        elif status == _opencmiss_swig.cvar.CMISSPointerNotNULL:
            raise CMISSError, "CMISS type pointer is not null"
        elif status == _opencmiss_swig.cvar.CMISSCouldNotAllocatePointer:
            raise CMISSError, "Could not allocate pointer"
        elif status == _opencmiss_swig.cvar.CMISSErrorConvertingPointer:
            raise CMISSError, "Error converting pointer"
        else:
            raise CMISSError, _opencmiss_swig.CMISSExtractErrorMessage()[1]
    return return_val

