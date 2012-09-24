import numpy as np
import matplotlib as mpl

class XScalarFormatter(mpl.ticker.ScalarFormatter):
    """eXtended ScalarFormatter"""

    def __init__(self, useOffset = True, useMathText = False, \
        useLocale = None, precision = None):

        mpl.ticker.ScalarFormatter.__init__(\
            self, useOffset, useMathText, useLocale)
        self._precision = precision

    def set_precision(self, precision):
        """use this function to set precision"""
        self._precision = precision

    def get_offset(self):
        """Return scientific notation, plus offset"""
        if len(self.locs)==0: return ''
        s = ''
        if self.orderOfMagnitude or self.offset:
            offsetStr = ''
            sciNotStr = ''
            if self.offset:
                offsetStr = self.format_data(self.offset)
                if self.offset > 0: offsetStr = '+' + offsetStr
            if self.orderOfMagnitude:
                if self._usetex or self._useMathText:
                    sciNotStr = self.format_data(10**self.orderOfMagnitude)
                else:
                    sciNotStr = '1e%d'% self.orderOfMagnitude
            if self._usetex or self._useMathText:
                if sciNotStr != '':
                    sciNotStr = r'\times%s' % sciNotStr
                s =  ''.join(('$',sciNotStr,offsetStr,'$'))
            else:
                s =  ''.join((sciNotStr,offsetStr))

        return self.fix_minus(s)

    def _set_format(self):
        # set the format string to format all the ticklabels
        # The floating point black magic (adding 1e-15 and formatting
        # to 8 digits) may warrant review and cleanup.
        locs = (np.asarray(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
        sigfigs = [len(str('%1.8f'% loc).split('.')[1].rstrip('0')) \
                   for loc in locs]
        sigfigs.sort()
        if self._precision is None:
            self.format = '%1.' + str(sigfigs[-1]) + 'f'
        else:
            self.format = '%1.' + str(self._precision) + 'f'

        if self._usetex or self._useMathText:
            self.format = '$%s$' % self.format

