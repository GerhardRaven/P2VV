# http://www.math.univ-toulouse.fr/Archive-MIP/publis/files/06.34.pdf
# http://www.sciencedirect.com/science?_ob=MiamiImageURL&_cid=271382&_user=6891091&_pii=S016612809690531X&_check=y&_origin=&_coverDate=27-Sep-1996&view=c&wchp=dGLbVlB-zSkWb&md5=4d287100293db5e8c56ec03e48ef286e/1-s2.0-S016612809690531X-main.pdf
# http://www.sciencedirect.com/science/article/pii/S016612809690531X

# adapted from H.H.H.Homeier & E.O. Steinborn, "Some properties of the coupling coefficiencs
# of real spherical harmonics and their relation to Gaunt coefficients", Journal of Molecular Structure (Theochem) 368 (1996) 31-37

class lm :
    def __init__( self, l, m) :
        self.l = l
        self.m = m


# transformation between spherical harmonics and real spherical harmonics
#  Y_{l0} = Y_l^0   (m=0)
#  Y_{lm} = \frac{1}{\sqrt{2}} \left(        Y_l^{+m} + (-1)^m Y_l^{-m} \right) (m>0)
#  Y_{lm} = \frac{i}{\sqrt{2}} \left( (-1)^m Y_l^{+m} -        Y_l^{-m} \right) (m<0)
#
#  write as:
# 
# Y_{lm} = sum_mu _u(l,m,mu)Y_l^mu
def _u(l,m,mu) :
    if abs(m)!=abs(mu) : return complex(0,0)
    if m==0   :          return complex(1,0)
    from math import sqrt
    is2 = 1.0/sqrt(2)
    sgn = lambda x : +1 if x%2 == 0 else -1
    if m>0 : return complex(is2, 0) if mu>0 else complex(sgn(mu)*is2,0)
    else   : return complex(0,-is2) if mu>0 else complex(0,sgn(mu)*is2)


# TODO: move into RooSpHarmonic as 'static' functions???
def gaunt( l1, l2, l3, m1, m2, m3 ) :
    import math
    _n = lambda x: 2*x+1
    n = math.sqrt( _n(l1)*_n(l2)*_n(l3) / ( 4 * math.pi ) )
    if m1%2 == 1 : n = -n
    import ROOT
    MathMore = ROOT.MathMore
    _3j = ROOT.Math.wigner_3j  # WARNING: wigner_3j uses half-integer units!!!!
    w1 = _3j( 2*l1, 2*l2, 2*l3,  0,  0,  0 )  
    w2 = _3j( 2*l1, 2*l2, 2*l3, -2*m1, 2*m2, 2*m3 )
    #print l1,l2,l3,m1,m2,m3,' -> ',n,w1,w2
    return  n*w1*w2

def real_gaunt( l1, l2, l3, m1, m2, m3 ) :
     c = sorted( [ lm(l1,m1), lm(l2,m2), lm(l3,m3) ], key = lambda x: abs(x.m) )
     if c[2].m == 0 : # m1=m2=m3=0
         return gaunt( c[0].l, c[1].l, c[2].l, c[0].m, c[1].m, c[2].m ) if c[0].m == 0 else 0 
     elif c[1].m == 0 : # m1=m2=0, m3!=0
         g = gaunt( c[0].l, c[1].l, c[2].l, c[0].m, c[1].m, c[2].m )
         u = ( _u(c[0].l,c[0].m,c[1].m ).conjugate()*_u(c[1].l,c[1].m,c[1].m ) ).real
         return 2*g*u
     else :
         u1 = (_u(c[0].l,c[0].m,c[1].m+c[2].m).conjugate()*_u(c[1].l,c[1].m,c[1].m)*_u(c[2].l,c[2].m, c[2].m)).real
         g1 = gaunt( c[0].l, c[1].l, c[2].l, c[1].m+c[2].m, c[1].m,  c[2].m )
         u2 = (_u(c[0].l,c[0].m,c[1].m-c[2].m).conjugate()*_u(c[1].l,c[1].m,c[1].m)*_u(c[2].l,c[2].m,-c[2].m)).real
         g2 = gaunt( c[0].l, c[1].l, c[2].l, c[1].m-c[2].m, c[1].m, -c[2].m )
         return 2*(g1*u1+g2*u2)

def non_zero_coupling_coefficients( l1,m1,l2,m2 ) :
    l_max = l1+l2
    l_min = max(abs(l1-l2),abs(m1+m2))
    if (l_min+l_max)%2 : l_min += 1
    #print 'selected lm for ',l1,m1,l2,m2,' : ', [ (l,m1+m2) for l in xrange( l_min,l_max+1,2) ]
    return( (l,m1+m2) for l in xrange( l_min,l_max+1,2) )


def non_zero_real_coupling_coefficients( l1,m1,l2,m2 ) :
    # (37)
    l_max = l1+l2
    # (38)
    l_min = max( abs(l1-l2),min( abs(m1+m2),abs(m1-m2) ) )
    if ( l_min+l_max ) % 2  : l_min += 1
    # note: odd number of negative (m,m1,m2) give also zero...
    #_range = lambda l : xrange(0,l+1) if m1*m2<0 else xrange(-l,1)
    _range = lambda l : xrange(-l,1) if m1*m2<0 else xrange(0,l+1)
    return ( (l,m) for l in xrange( l_min, l_max+1, 2 ) for m in _range( l ) )
            

# CLAIM: for real spherical harmonics Y_l1,m1 * Y_l2,m2 = sum c * Y_l,m for (c,l,m) in gaunt_expansion(l1,m1, l2,m2)
#        where c = int Y_l,m Y_l1,m1 Y_l2,m2
def real_gaunt_expansion(l1,m1, l2,m2) :
    return filter( lambda x : x[-1]!=0 
                 , ( (l,m,real_gaunt(l,l1,l2,m,m1,m2)) for (l,m) in non_zero_real_coupling_coefficients(l1,m1,l2,m2) ) 
                 )


def gaunt_expansion(l1,m1, l2,m2) :
    return filter( lambda x : x[-1]!=0
                 , ( (l,m,gaunt(l,l1,l2,m,m1,m2)) for (l,m) in non_zero_coupling_coefficients(l1,m1,l2,m2) ) 
                 )
