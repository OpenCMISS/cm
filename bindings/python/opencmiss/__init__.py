<<<<<<< HEAD

# OpenCMISS Python package initialisation file.
from pkgutil import extend_path
__path__ = extend_path(__path__, __name__)

=======
from pkgutil import extend_path


# Update the search path for this package to include other
# directories in the python path with the same name.
# This means that opencmiss.zinc can be in a different directory
# to opencmiss.iron and both can be imported.
__path__ = extend_path(__path__, __name__)
>>>>>>> e8f6663d0323aafdbd166678727f87b89fe300d3
