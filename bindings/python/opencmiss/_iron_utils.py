"""Utility routines and classes used by OpenCMISS-Iron
"""

from . import _iron_swig


class CMISSError(Exception):
    """Base class for errors in the OpenCMISS library"""

    pass


class CMISSType(object):
    """Base class for all OpenCMISS types"""

    pass


class Enum(object):
    pass


def wrap_cmiss_routine(routine, args=None):
    """Wrap a call to the OpenCMISS SWIG module

    Call the routine and check the return value, and raise an
    exception if it is non-zero.

    Return any other remaining return values.

    """
    if args is None:
        r = routine()
    else:
        #Replace wrapped cmiss types with the underlying type
        new_args = []
        for arg in args:
            if hasattr(arg, 'cmiss_type'):
                new_args.append(arg.cmiss_type)
            else:
                try:
                    # Try to convert a list of CMISS types first.
                    # Check length to avoid empty strings being converted
                    # to an empty list
                    if len(arg) > 0:
                        new_args.append([a.cmiss_type for a in arg])
                    else:
                        new_args.append(arg)
                except (TypeError, AttributeError):
                    new_args.append(arg)
        r = routine(*new_args)
    # We will either have a list of multiple return values, or
    # a single status code as a return. Don't have to worry about
    # ever having a single return value as a list as there will always
    # be at least a return status.
    if isinstance(r, list):
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
    if status != _iron_swig.cvar.CMFE_NO_ERROR:
        if status == _iron_swig.cvar.CMFE_POINTER_IS_NULL:
            raise CMISSError("CMISS type pointer is null")
        elif status == _iron_swig.cvar.CMFE_POINTER_NOT_NULL:
            raise CMISSError("CMISS type pointer is not null")
        elif status == _iron_swig.cvar.CMFE_COULD_NOT_ALLOCATE_POINTER:
            raise CMISSError("Could not allocate pointer")
        elif status == _iron_swig.cvar.CMFE_ERROR_CONVERTING_POINTER:
            raise CMISSError("Error converting pointer")
        else:
            raise CMISSError(_iron_swig.cmfe_ExtractErrorMessage()[1])
    return return_val
