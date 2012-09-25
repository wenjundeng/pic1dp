# Copyright 2012 Wenjun Deng <wdeng@wdeng.info>
#
# This file is part of PIC1D-PETSc
#
# PIC1D-PETSc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PIC1D-PETSc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PIC1D-PETSc.  If not, see <http://www.gnu.org/licenses/>.


"""Provides an eXtended PetscBinaryIO class"""

from PetscBinaryIO import *

class XPetscBinaryIO(PetscBinaryIO):
    """eXtended PetscBinaryIO class"""

    def _update_dtypes(self):

        # define self._inttype and self._scalartype
        PetscBinaryIO._update_dtypes(self)

        if self.complexscalars:
            # fixing a bug in PetscBinaryIO,
            # complex number takes twice as many bytes than real number
            if self.precision == 'longlong':
                self._scalartype = '>c32'
            if self.precision == 'single':
                self._scalartype = '>c8'
            else:
                self._scalartype = '>c16'

            # define self._realtype
            if self.precision == 'longlong':
                self._realtype = '>f16'
            elif self.precision == 'single':
                self._realtype = '>f4'
            else:
                self._realtype = '>f8'
        else:
            # define self._realtype
            self._realtype = self._scalartype

    @decorate_with_conf
    def readInt(self, fh, n):
        """Reads n integers from binary file handle, returning an array of the data."""

        vals = np.fromfile(fh, dtype=self._inttype, count=n)
        return vals

    @decorate_with_conf
    def readReal(self, fh, n):
        """Reads n real numbers from binary file handle, returning an array of the data."""

        vals = np.fromfile(fh, dtype=self._realtype, count=n)
        return vals

    @decorate_with_conf
    def readScalar(self, fh, n):
        """Reads n scalars from binary file handle, returning an array of the data."""

        vals = np.fromfile(fh, dtype=self._scalartype, count=n)
        return vals

