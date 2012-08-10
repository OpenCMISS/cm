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


def Matrix_DataGet(self, *args):
    routines = {
        MatrixVectorDataTypes.INTG: self.DataGetIntg,
        MatrixVectorDataTypes.SP: self.DataGetSP,
        MatrixVectorDataTypes.DP: self.DataGetDP,
        MatrixVectorDataTypes.L: self.DataGetL,
    }
    data_type = self.DataTypeGet()
    return routines[data_type](*args)


def Matrix_DataRestore(self, *args):
    routines = {
        MatrixVectorDataTypes.INTG: self.DataRestoreIntg,
        MatrixVectorDataTypes.SP: self.DataRestoreSP,
        MatrixVectorDataTypes.DP: self.DataRestoreDP,
        MatrixVectorDataTypes.L: self.DataRestoreL,
    }
    data_type = self.DataTypeGet()
    return routines[data_type](*args)


def Vector_DataGet(self, *args):
    routines = {
        MatrixVectorDataTypes.INTG: self.DataGetIntg,
        MatrixVectorDataTypes.SP: self.DataGetSP,
        MatrixVectorDataTypes.DP: self.DataGetDP,
        MatrixVectorDataTypes.L: self.DataGetL,
    }
    data_type = self.DataTypeGet()
    return routines[data_type](*args)


def Vector_DataRestore(self, *args):
    routines = {
        MatrixVectorDataTypes.INTG: self.DataRestoreIntg,
        MatrixVectorDataTypes.SP: self.DataRestoreSP,
        MatrixVectorDataTypes.DP: self.DataRestoreDP,
        MatrixVectorDataTypes.L: self.DataRestoreL,
    }
    data_type = self.DataTypeGet()
    return routines[data_type](*args)


Matrix.DataGet = Matrix_DataGet
Matrix.DataRestore = Matrix_DataRestore
Vector.DataGet = Vector_DataGet
Vector.DataRestore = Vector_DataRestore


def Matrix_ToSciPy(self):
    """Return a SciPy matrix representation of this matrix

    This works with sparse and full matrices and uses a view
    of the matrix data so there is no copying.
    Once finished with the matrix you should call the
    SciPyRestore method.
    """

    # Import scipy here as we don't want to require it unless
    # people are actually going to use it
    import numpy
    from scipy import sparse

    storageType = self.StorageTypeGet()
    dimensions = self.DimensionsGet()
    data = self.DataGet()

    if storageType == MatrixStorageTypes.BLOCK:
        # Not sparse, so just reshape the data
        matrix = data.reshape(dimensions, order='F')
    elif storageType == MatrixStorageTypes.DIAGONAL:
        offsets = numpy.array([0])
        matrix = sparse.dia_matrix((data, offsets), shape=dimensions)
    elif storageType == MatrixStorageTypes.COMPRESSED_ROW:
        # OpenCMISS has two types of distributed matrices, ones used internally
        # and ones used by PETSc. The PETSc ones use zero based arrays but the
        # internal ones use one based arrays. So for these to work with SciPy
        # we need to subtract one from the row and column index arrays.
        # This requires allocating a new array but it will take much less space
        # than the data array so this souldn't be an issue.
        # As the first row index is always zero, we can check if it is one instead
        rowIndices, columnIndices = self.StorageLocationsGet()
        if rowIndices[0] == 1:
            rowIndices = rowIndices - 1
            columnIndices = columnIndices - 1
        else:
            # Still copy these so that we're not hanging on do
            # data allocated within OpenCMISS, it's not that expensive
            rowIndices = rowIndices.copy()
            columnIndices = columnIndices.copy()
        matrix = sparse.csr_matrix(
                (data, columnIndices, rowIndices), shape=dimensions)
    elif storageType == MatrixStorageTypes.COMPRESSED_COLUMN:
        rowIndices, columnIndices = self.StorageLocationsGet()
        # As the first column index is always zero, we can check if it is one instead
        if columnIndices[0] == 1:
            rowIndices = rowIndices - 1
            columnIndices = columnIndices - 1
        else:
            rowIndices = rowIndices.copy()
            columnIndices = columnIndices.copy()
        matrix = sparse.csc_matrix(
                (data, rowIndices, columnIndices), shape=dimensions)
    else:
        self.DataRestore(data)
        raise ValueError("The storage type for this matrix is not "
            "supported by SciPy")
    return matrix


def Matrix_SciPyRestore(self, matrix):
    """Restores the data pointers used when creating a SciPy matrix

    Trying to use the SciPy matrix after this will not work
    """

    import numpy

    dimensions = self.DimensionsGet()
    if isinstance(matrix, numpy.ndarray):
        # For full matrices
        # This doesn't actually copy data if the
        # order='F' is used:
        matrix = numpy.reshape(matrix, -1, order='F')
        self.DataRestore(matrix)
    else:
        # For sparse matrices
        self.DataRestore(matrix.data)


Matrix.ToSciPy = Matrix_ToSciPy
Matrix.SciPyRestore = Matrix_SciPyRestore
