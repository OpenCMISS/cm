# Redirection for field routines that have separate methods
# for different data types


def Field_ParameterSetDataGet(self, *args):
    variableType, fieldSetType = args
    routines = {
        FieldDataTypes.INTG: self.ParameterSetDataGetIntg,
        FieldDataTypes.SP: self.ParameterSetDataGetSP,
        FieldDataTypes.DP: self.ParameterSetDataGetDP,
        FieldDataTypes.L: self.ParameterSetDataGetL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetDataRestore(self, *args):
    variableType, fieldSetType, parameters = args
    routines = {
        FieldDataTypes.INTG: self.ParameterSetDataRestoreIntg,
        FieldDataTypes.SP: self.ParameterSetDataRestoreSP,
        FieldDataTypes.DP: self.ParameterSetDataRestoreDP,
        FieldDataTypes.L: self.ParameterSetDataRestoreL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ComponentValuesInitialise(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ComponentValuesInitialiseIntg,
        FieldDataTypes.SP: self.ComponentValuesInitialiseSP,
        FieldDataTypes.DP: self.ComponentValuesInitialiseDP,
        FieldDataTypes.L: self.ComponentValuesInitialiseL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetGetConstant(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetGetConstantIntg,
        FieldDataTypes.SP: self.ParameterSetGetConstantSP,
        FieldDataTypes.DP: self.ParameterSetGetConstantDP,
        FieldDataTypes.L: self.ParameterSetGetConstantL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetGetElement(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetGetElementIntg,
        FieldDataTypes.SP: self.ParameterSetGetElementSP,
        FieldDataTypes.DP: self.ParameterSetGetElementDP,
        FieldDataTypes.L: self.ParameterSetGetElementL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetGetGaussPoint(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetGetGaussPointIntg,
        FieldDataTypes.SP: self.ParameterSetGetGaussPointSP,
        FieldDataTypes.DP: self.ParameterSetGetGaussPointDP,
        FieldDataTypes.L: self.ParameterSetGetGaussPointL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetGetNode(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetGetNodeIntg,
        FieldDataTypes.SP: self.ParameterSetGetNodeSP,
        FieldDataTypes.DP: self.ParameterSetGetNodeDP,
        FieldDataTypes.L: self.ParameterSetGetNodeL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetUpdateConstant(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetUpdateConstantIntg,
        FieldDataTypes.SP: self.ParameterSetUpdateConstantSP,
        FieldDataTypes.DP: self.ParameterSetUpdateConstantDP,
        FieldDataTypes.L: self.ParameterSetUpdateConstantL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetUpdateElement(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetUpdateElementIntg,
        FieldDataTypes.SP: self.ParameterSetUpdateElementSP,
        FieldDataTypes.DP: self.ParameterSetUpdateElementDP,
        FieldDataTypes.L: self.ParameterSetUpdateElementL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetUpdateGaussPoint(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetUpdateGaussPointIntg,
        FieldDataTypes.SP: self.ParameterSetUpdateGaussPointSP,
        FieldDataTypes.DP: self.ParameterSetUpdateGaussPointDP,
        FieldDataTypes.L: self.ParameterSetUpdateGaussPointL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetUpdateNode(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetUpdateNodeIntg,
        FieldDataTypes.SP: self.ParameterSetUpdateNodeSP,
        FieldDataTypes.DP: self.ParameterSetUpdateNodeDP,
        FieldDataTypes.L: self.ParameterSetUpdateNodeL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetAddConstant(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetAddConstantIntg,
        FieldDataTypes.SP: self.ParameterSetAddConstantSP,
        FieldDataTypes.DP: self.ParameterSetAddConstantDP,
        FieldDataTypes.L: self.ParameterSetAddConstantL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetAddElement(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetAddElementIntg,
        FieldDataTypes.SP: self.ParameterSetAddElementSP,
        FieldDataTypes.DP: self.ParameterSetAddElementDP,
        FieldDataTypes.L: self.ParameterSetAddElementL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetAddNode(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.INTG: self.ParameterSetAddNodeIntg,
        FieldDataTypes.SP: self.ParameterSetAddNodeSP,
        FieldDataTypes.DP: self.ParameterSetAddNodeDP,
        FieldDataTypes.L: self.ParameterSetAddNodeL,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetGetGuassPoint(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.DP: self.ParameterSetGetGaussPointDP,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


def Field_ParameterSetInterpolateGauss(self, *args):
    variableType = args[0]
    routines = {
        FieldDataTypes.DP: self.ParameterSetInterpolateGaussDP,
    }
    data_type = self.DataTypeGet(variableType)
    return routines[data_type](*args)


Field.ParameterSetDataGet = Field_ParameterSetDataGet
Field.ParameterSetDataRestore = Field_ParameterSetDataRestore
Field.ComponentValuesInitialise = Field_ComponentValuesInitialise

Field.ParameterSetGetConstant = Field_ParameterSetGetConstant
Field.ParameterSetGetElement = Field_ParameterSetGetElement
Field.ParameterSetGetGaussPoint = Field_ParameterSetGetGaussPoint
Field.ParameterSetGetNode = Field_ParameterSetGetNode

Field.ParameterSetUpdateConstant = Field_ParameterSetUpdateConstant
Field.ParameterSetUpdateElement = Field_ParameterSetUpdateElement
Field.ParameterSetUpdateGaussPoint = Field_ParameterSetUpdateGaussPoint
Field.ParameterSetUpdateNode = Field_ParameterSetUpdateNode

Field.ParameterSetAddConstant = Field_ParameterSetAddConstant
Field.meterSetAddElement = Field_ParameterSetAddElement
Field.ParameterSetAddNode = Field_ParameterSetAddNode

Field.ParameterSetGetGuassPoint = Field_ParameterSetGetGuassPoint
Field.ParameterSetInterpolateGauss = Field_ParameterSetInterpolateGauss
