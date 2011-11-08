# http://www.math.univ-toulouse.fr/Archive-MIP/publis/files/06.34.pdf
# http://www.sciencedirect.com/science?_ob=MiamiImageURL&_cid=271382&_user=6891091&_pii=S016612809690531X&_check=y&_origin=&_coverDate=27-Sep-1996&view=c&wchp=dGLbVlB-zSkWb&md5=4d287100293db5e8c56ec03e48ef286e/1-s2.0-S016612809690531X-main.pdf
# http://www.sciencedirect.com/science/article/pii/S016612809690531X


# TODO: move into RooSpHarmonic as 'static' functions???
def gaunt( l1, l2, l3, m1, m2, m3 ) :
    import math
    _n = lambda x: 2*x+1
    n = sqrt( _n(l1)*_n(l2)*_n(l3) / ( 4 * math.pi ) )
    import ROOT
    _3j = ROOT.Math.wigner_3j
    w1 = _3j( l1, l2, l3,  0,  0,  0 )
    w2 = _3j( l1, l2, l3, m1, m2, m3 )
    return  n*w1*w2


def non_zero_coupling_coefficients( l1,m1,l2,m2 ) :
    l_max = l1+l2
    l_min = max( abs(l1-l2),min( abs(m1+m2),abs(m1-m2) ) )
    if ( l_min+l_max ) % 2  : l_min += 1
    return ( (l,m) for l in range( l_min, l_max+1, 2 ) for m in range( -l, l+1 ) )
            

