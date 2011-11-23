class CheckLabels(object):
    def __init__(self, operators, operands):
        self._labels = set()
        for l, m, r in operators:
            for o in operands:
                if not hasattr(operands, '__iter__'):
                    self._labels.add('%s%s%s' % (l, o, r))
                elif len(operands) == 1:
                    self._labels.add('%s%s%s' % (l, o[0], r))
                elif len(operands) == 2:
                    self._labels.add('%s%s%s%s%s' % (l, o[1], m, o[2], r))

    def checkLabel(self, l):
        return l in self._labels

    def labels(self):
        return frozenset(self._labels)
