import parse


class Mock(object):
    def __init__(self, **kwargs):
        self.name = "test"
        self.comment = "Test comment"
        self.__dict__.update(kwargs)


# Scalars
input_integer = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='IN',
        pointer=False,
        array_dims=0,
        required_sizes=0,
        array_spec=[])

output_integer = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='OUT',
        pointer=False,
        array_dims=0,
        required_sizes=0,
        array_spec=[])

input_real = Mock(
        var_type=parse.Parameter.DOUBLE,
        intent='IN',
        pointer=False,
        array_dims=0,
        required_sizes=0,
        array_spec=[])

output_real = Mock(
        var_type=parse.Parameter.DOUBLE,
        intent='OUT',
        pointer=False,
        array_dims=0,
        required_sizes=0,
        array_spec=[])

# Arrays
input_array = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='IN',
        pointer=False,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

input_array_2d = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='IN',
        pointer=False,
        array_dims=2,
        required_sizes=2,
        array_spec=[':', ':'])

output_array_known_size = Mock(
        var_type=parse.Parameter.DOUBLE,
        intent='OUT',
        pointer=False,
        array_dims=1,
        required_sizes=0,
        array_spec=['2'])

output_array = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='OUT',
        pointer=False,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

input_array_pointer = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='IN',
        pointer=True,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

output_array_pointer = Mock(
        var_type=parse.Parameter.INTEGER,
        intent='OUT',
        pointer=True,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

# Strings
input_string = Mock(
        intent='IN',
        pointer=False,
        var_type=parse.Parameter.CHARACTER,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

output_string = Mock(
        intent='OUT',
        pointer=False,
        var_type=parse.Parameter.CHARACTER,
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])

input_string_array = Mock(
        intent='IN',
        pointer=False,
        var_type=parse.Parameter.CHARACTER,
        array_dims=2,
        required_sizes=2,
        array_spec=[':', ':'])

# CMFE Types
input_cmiss_type = Mock(
        intent='IN',
        pointer=False,
        var_type=parse.Parameter.CUSTOM_TYPE,
        type_name="cmfe_TestType",
        array_dims=0,
        required_sizes=0,
        array_spec=[])

output_cmiss_type = Mock(
        intent='OUT',
        pointer=False,
        var_type=parse.Parameter.CUSTOM_TYPE,
        type_name="cmfe_TestType",
        array_dims=0,
        required_sizes=0,
        array_spec=[])

input_cmiss_type_array = Mock(
        intent='IN',
        pointer=False,
        var_type=parse.Parameter.CUSTOM_TYPE,
        type_name="cmfe_TestType",
        array_dims=1,
        required_sizes=1,
        array_spec=[':'])
