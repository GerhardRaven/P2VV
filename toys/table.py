from math import sqrt

def latex_table(name, values, gen_params):
    f = open('params_%s.tex' % name, 'w')

    string = '\\documentclass{article}\n'
    string += '\\begin{document}\n'

    def __s(n, fmt = '%.2e'):
        if type(n) == float:
            n = fmt % n
        if n.endswith('e+00'):
            return n[:-4]
        else:
            return n[:-4] + (' * 10^{%d}' % int(n[-3:]))

    offset = 0
    for name in values.iterkeys():
        name = name.replace('_', '\_')
        if len(name) > offset:
            offset = len(name)
    fmt = '{0:<' + str(offset) + '} & '

    string += '\\begin{tabular}{|c|c|c|c|}\n'
    string += '\\hline\n'
    string += 'parameter & bias & original value & bias/error \\\\ \n'
    string += '\\hline\n'
    string += '\\hline\n'
    for name, (mean, sigma) in values.iteritems():
        line = ''
        try:
            gen_param = gen_params.find(name[:-5])
            line += fmt.format(name[:-5].replace('_', '\_'))
            m = (mean.getVal(), mean.getError())
            s = (sigma.getVal(), sigma.getError())
            b = m[0] * s[0]
            e = sqrt((m[1]/m[0])**2 - (s[1] / s[0])**2) * b
            line += '$%s \pm %s$ & ' % (__s(b), __s(e))
            line += '$%s$ & ' % __s(gen_param.getVal())
            line +=  '$%s$ \\\\ \n' % __s(b / e)
        except ZeroDivisionError:
            pass
        else:
            string += line
    string += '\\hline\n'
    string += '\\end{tabular}\n'
    string += '\\end{document}\n'

    f.write(string)
    f.close()
