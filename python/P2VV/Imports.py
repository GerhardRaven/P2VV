###########################################################################################################################################
## P2VVImports                                                                                                                           ##
##                                                                                                                                       ##
## authors':                                                                                                                             ##
##   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                                                                            ##
##                                                                                                                                       ##
###########################################################################################################################################


# parameter names dictionary:  'P2VV internal name' : ( 'text name', 'LaTeX name' )
parNames = {  'phiCP'                     : ( 'phi_s',                       '$\\phi_\\text{s}$'                          )
            , '__phiCP__'                 : ( 'phi_s (b)',                   '$\\phi_\\text{s}$ (b)'                      )
            , 'lambdaCP'                  : ( '|lambda_s|',                  '$|\\lambda_\\text{s}|$'                     )
            , '__lambdaCP__'              : ( '|lambda_s| (b)',              '$|\\lambda_\\text{s}|$ (b)'                 )
            , 'Gamma'                     : ( 'Gamma_s',                     '$\\Gamma_\\text{s}$'                        )
            , 'Gamma_p2011'               : ( 'Gamma_s - beta_2011',         '$\\Gamma_\\text{s} - \\beta_\\{2011}$'      )
            , 'Gamma_p2012'               : ( 'Gamma_s - beta_2012',         '$\\Gamma_\\text{s} - \\beta_\\{2012}$'      )
            , 'dGamma'                    : ( 'Delta Gamma_s',               '$\\Delta\\Gamma_\\text{s}$'                 )
            , '__dGamma__'                : ( 'Delta Gamma_s (b)',           '$\\Delta\\Gamma_\\text{s}$ (b)'             )
            , 'dM'                        : ( 'Delta m_s',                   '$\\Delta m_\\text{s}$'                      )
            , 'A0Mag2'                    : ( '|A_0|^2',                     '$|A_0|^2$'                                  )
            , 'AperpMag2'                 : ( '|A_perp|^2',                  '$|A_\\perp|^2$'                             )
            , 'f_S'                       : ( 'F_S',                         '$F_\\text{S}$'                              )
            , 'f_S_bin0'                  : ( 'F_S_0',                       '${F_\\text{S}}_0$'                          )
            , 'f_S_bin1'                  : ( 'F_S_1',                       '${F_\\text{S}}_1$'                          )
            , 'f_S_bin2'                  : ( 'F_S_2',                       '${F_\\text{S}}_2$'                          )
            , 'f_S_bin3'                  : ( 'F_S_3',                       '${F_\\text{S}}_3$'                          )
            , 'f_S_bin4'                  : ( 'F_S_4',                       '${F_\\text{S}}_4$'                          )
            , 'f_S_bin5'                  : ( 'F_S_5',                       '${F_\\text{S}}_5$'                          )
            , 'AparPhase'                 : ( 'delta_para - delta_0',        '$\\delta_\\parallel-\\delta_0$'             )
            , 'AperpPhase'                : ( 'delta_perp - delta_0',        '$\\delta_\\perp-\\delta_0$'                 )
            , 'ASOddPhase'                : ( 'delta_S - delta_perp',        '$\\delta_\\text{S}-\\delta_\\perp$'         )
            , 'ASOddPhase_bin0'           : ( 'delta_S_0 - delta_perp',      '${\\delta_\\text{S}}_0-\\delta_\\perp$'     )
            , 'ASOddPhase_bin1'           : ( 'delta_S_1 - delta_perp',      '${\\delta_\\text{S}}_1-\\delta_\\perp$'     )
            , 'ASOddPhase_bin2'           : ( 'delta_S_2 - delta_perp',      '${\\delta_\\text{S}}_2-\\delta_\\perp$'     )
            , 'ASOddPhase_bin3'           : ( 'delta_S_3 - delta_perp',      '${\\delta_\\text{S}}_3-\\delta_\\perp$'     )
            , 'ASOddPhase_bin4'           : ( 'delta_S_4 - delta_perp',      '${\\delta_\\text{S}}_4-\\delta_\\perp$'     )
            , 'ASOddPhase_bin5'           : ( 'delta_S_5 - delta_perp',      '${\\delta_\\text{S}}_5-\\delta_\\perp$'     )
            , 'timeResSigmaSF'            : ( 'S_sigma_t',                   '$S_{\\sigma_{t}}$'                          )
            , 'betaTimeEff'               : ( 'beta',                        '$\\beta$'                                   )
            , 'betaTimeEff_p2011'         : ( 'beta_2011',                   '$\\beta_{2011}$'                            )
            , 'betaTimeEff_p2012'         : ( 'beta_2012',                   '$\\beta_{2012}$'                            )
            , 'wTagP0OS'                  : ( 'p_0 OS',                      '$p_0$ OS'                                   )
            , 'wTagP1OS'                  : ( 'p_1 OS',                      '$p_1$ OS'                                   )
            , 'wTagDelP0OS'               : ( 'Delta p_0 OS',                '$\\Delta p_0$ OS'                           )
            , 'wTagDelP1OS'               : ( 'Delta p_1 OS',                '$\\Delta p_1$ OS'                           )
            , 'wTagP0SS'                  : ( 'p_0 SS',                      '$p_0$ SS'                                   )
            , 'wTagP1SS'                  : ( 'p_1 SS',                      '$p_1$ SS'                                   )
            , 'wTagDelP0SS'               : ( 'Delta p_0 SS',                '$\\Delta p_0$ SS'                           )
            , 'wTagDelP1SS'               : ( 'Delta p_1 SS',                '$\\Delta p_1$ SS'                           )
            , 'N_sigMass'                 : ( 'N_sig',                       '$N_\\text{sig}$'                            )
            , 'N_sigMass_bin0'            : ( 'N_sig_0',                     '${N_\\text{sig}}_0$'                        )
            , 'N_sigMass_bin1'            : ( 'N_sig_1',                     '${N_\\text{sig}}_1$'                        )
            , 'N_sigMass_bin2'            : ( 'N_sig_2',                     '${N_\\text{sig}}_2$'                        )
            , 'N_sigMass_bin3'            : ( 'N_sig_3',                     '${N_\\text{sig}}_3$'                        )
            , 'N_sigMass_bin4'            : ( 'N_sig_4',                     '${N_\\text{sig}}_4$'                        )
            , 'N_sigMass_bin5'            : ( 'N_sig_5',                     '${N_\\text{sig}}_5$'                        )
            , 'N_sigMass_{notExclB;bin0}' : ( 'N_sig_U0',                    '${N_\\text{sig}}_{U0}$'                     )
            , 'N_sigMass_{notExclB;bin1}' : ( 'N_sig_U1',                    '${N_\\text{sig}}_{U1}$'                     )
            , 'N_sigMass_{notExclB;bin2}' : ( 'N_sig_U2',                    '${N_\\text{sig}}_{U2}$'                     )
            , 'N_sigMass_{notExclB;bin3}' : ( 'N_sig_U3',                    '${N_\\text{sig}}_{U3}$'                     )
            , 'N_sigMass_{notExclB;bin4}' : ( 'N_sig_U4',                    '${N_\\text{sig}}_{U4}$'                     )
            , 'N_sigMass_{notExclB;bin5}' : ( 'N_sig_U5',                    '${N_\\text{sig}}_{U5}$'                     )
            , 'N_sigMass_{exclB;bin0}'    : ( 'N_sig_B0',                    '${N_\\text{sig}}_{B0}$'                     )
            , 'N_sigMass_{exclB;bin1}'    : ( 'N_sig_B1',                    '${N_\\text{sig}}_{B1}$'                     )
            , 'N_sigMass_{exclB;bin2}'    : ( 'N_sig_B2',                    '${N_\\text{sig}}_{B2}$'                     )
            , 'N_sigMass_{exclB;bin3}'    : ( 'N_sig_B3',                    '${N_\\text{sig}}_{B3}$'                     )
            , 'N_sigMass_{exclB;bin4}'    : ( 'N_sig_B4',                    '${N_\\text{sig}}_{B4}$'                     )
            , 'N_sigMass_{exclB;bin5}'    : ( 'N_sig_B5',                    '${N_\\text{sig}}_{B5}$'                     )
            , 'N_bkgMass'                 : ( 'N_bkg',                       '$N_\\text{bkg}$'                            )
            , 'N_bkgMass_bin0'            : ( 'N_bkg_0',                     '${N_\\text{bkg}}_0$'                        )
            , 'N_bkgMass_bin1'            : ( 'N_bkg_1',                     '${N_\\text{bkg}}_1$'                        )
            , 'N_bkgMass_bin2'            : ( 'N_bkg_2',                     '${N_\\text{bkg}}_2$'                        )
            , 'N_bkgMass_bin3'            : ( 'N_bkg_3',                     '${N_\\text{bkg}}_3$'                        )
            , 'N_bkgMass_bin4'            : ( 'N_bkg_4',                     '${N_\\text{bkg}}_4$'                        )
            , 'N_bkgMass_bin5'            : ( 'N_bkg_5',                     '${N_\\text{bkg}}_5$'                        )
            , 'N_bkgMass_{notExclB;bin0}' : ( 'N_bkg_U0',                    '${N_\\text{bkg}}_{U0}$'                     )
            , 'N_bkgMass_{notExclB;bin1}' : ( 'N_bkg_U1',                    '${N_\\text{bkg}}_{U1}$'                     )
            , 'N_bkgMass_{notExclB;bin2}' : ( 'N_bkg_U2',                    '${N_\\text{bkg}}_{U2}$'                     )
            , 'N_bkgMass_{notExclB;bin3}' : ( 'N_bkg_U3',                    '${N_\\text{bkg}}_{U3}$'                     )
            , 'N_bkgMass_{notExclB;bin4}' : ( 'N_bkg_U4',                    '${N_\\text{bkg}}_{U4}$'                     )
            , 'N_bkgMass_{notExclB;bin5}' : ( 'N_bkg_U5',                    '${N_\\text{bkg}}_{U5}$'                     )
            , 'N_bkgMass_{exclB;bin0}'    : ( 'N_bkg_B0',                    '${N_\\text{bkg}}_{B0}$'                     )
            , 'N_bkgMass_{exclB;bin1}'    : ( 'N_bkg_B1',                    '${N_\\text{bkg}}_{B1}$'                     )
            , 'N_bkgMass_{exclB;bin2}'    : ( 'N_bkg_B2',                    '${N_\\text{bkg}}_{B2}$'                     )
            , 'N_bkgMass_{exclB;bin3}'    : ( 'N_bkg_B3',                    '${N_\\text{bkg}}_{B3}$'                     )
            , 'N_bkgMass_{exclB;bin4}'    : ( 'N_bkg_B4',                    '${N_\\text{bkg}}_{B4}$'                     )
            , 'N_bkgMass_{exclB;bin5}'    : ( 'N_bkg_B5',                    '${N_\\text{bkg}}_{B5}$'                     )
            , 'm_sig_frac'                : ( 'f_mass_1',                    '${f_1}_\\text{mass}$'                       )
            , 'm_sig_mean'                : ( 'mu_mass',                     '$\\mu_\\text{mass}$'                        )
            , 'm_sig_sigma_sf'            : ( 'sigma_mass 2:1 SF',           '$\\text{SF}_{2:1}\\ \\sigma_\\text{mass}$'  )
            , 'm_sig_sigma_1'             : ( 'sigma_mass_1',                '${\\sigma_1}_\\text{mass}$'                 )
            , 'm_sig_sigma_2'             : ( 'sigma_mass_2',                '${\\sigma_2}_\\text{mass}$'                 )
            , 'm_sig_widthPar0'           : ( 'width_par_mass_0',            '${\\text{width}_0}_\\text{mass}$'           )
            , 'm_sig_widthPar1'           : ( 'width_par_mass_1',            '${\\text{width}_1}_\\text{mass}$'           )
            , 'm_sig_widthPar2'           : ( 'width_par_mass_2',            '${\\text{width}_2}_\\text{mass}$'           )
            , 'm_bkg_exp'                 : ( 'alpha_mass',                  '$\\alpha_\\text{mass}$'                     )
            , 'm_bkg_exp_bin0'            : ( 'alpha_mass_0',                '${\\alpha_\\text{mass}}_0$'                 )
            , 'm_bkg_exp_bin1'            : ( 'alpha_mass_1',                '${\\alpha_\\text{mass}}_1$'                 )
            , 'm_bkg_exp_bin2'            : ( 'alpha_mass_2',                '${\\alpha_\\text{mass}}_2$'                 )
            , 'm_bkg_exp_bin3'            : ( 'alpha_mass_3',                '${\\alpha_\\text{mass}}_3$'                 )
            , 'm_bkg_exp_bin4'            : ( 'alpha_mass_4',                '${\\alpha_\\text{mass}}_4$'                 )
            , 'm_bkg_exp_bin5'            : ( 'alpha_mass_5',                '${\\alpha_\\text{mass}}_5$'                 )
            , 'm_bkg_arg'                 : ( 'alpha_mass',                  '$\\alpha_\\text{mass}$'                     )
            , 'm_bkg_arg_bin0'            : ( 'alpha_mass_0',                '${\\alpha_\\text{mass}}_0$'                 )
            , 'm_bkg_arg_bin1'            : ( 'alpha_mass_1',                '${\\alpha_\\text{mass}}_1$'                 )
            , 'm_bkg_arg_bin2'            : ( 'alpha_mass_2',                '${\\alpha_\\text{mass}}_2$'                 )
            , 'm_bkg_arg_bin3'            : ( 'alpha_mass_3',                '${\\alpha_\\text{mass}}_3$'                 )
            , 'm_bkg_arg_bin4'            : ( 'alpha_mass_4',                '${\\alpha_\\text{mass}}_4$'                 )
            , 'm_bkg_arg_bin5'            : ( 'alpha_mass_5',                '${\\alpha_\\text{mass}}_5$'                 )
           }

# nominal values of physics parameters
parValues = {  'A0Mag2'           : (  5.2077e-01,  3.44e-03, -1. )
             , 'ASOddPhase_bin0'  : (  8.0132e-01,  1.85e-01, -1. )
             , 'ASOddPhase_bin1'  : (  2.3158e+00,  1.96e-01, -1. )
             , 'ASOddPhase_bin2'  : (  4.5875e-01,  2.24e-01, -1. )
             , 'ASOddPhase_bin3'  : ( -3.5722e-01,  1.82e-01, -1. )
             , 'ASOddPhase_bin4'  : ( -6.7780e-01,  1.84e-01, -1. )
             , 'ASOddPhase_bin5'  : ( -8.2175e-01,  1.23e-01, -1. )
             , 'AparPhase'        : (  3.2583e+00,  1.24e-01, -1. )
             , 'AperpMag2'        : (  2.5360e-01,  4.92e-03, -1. )
             , 'AperpPhase'       : (  3.1311e+00,  1.14e-01, -1. )
             , 'Gamma_p2011'      : (  6.6900e-01,  4.45e-03, -1. )
             , 'Gamma_p2012'      : (  6.7442e-01,  3.20e-03, -1. )
             , '__dGamma__'       : (  6.7108e-02,  9.13e-03, -1. )
             , '__phiCP__'        : ( -1.2021e-01,  5.05e-02, -1. )
             , 'dM'               : (  1.7765e+01,  2.21e-02, -1. )
             , 'f_S_bin0'         : (  4.4388e-01,  5.49e-02, -1. )
             , 'f_S_bin1'         : (  6.3626e-02,  1.77e-02, -1. )
             , 'f_S_bin2'         : (  8.6575e-03,  6.46e-03, -1. )
             , 'f_S_bin3'         : (  8.9579e-03,  5.81e-03, -1. )
             , 'f_S_bin4'         : (  4.4347e-02,  1.56e-02, -1. )
             , 'f_S_bin5'         : (  2.0901e-01,  2.61e-02, -1. )
             , 'lambdaCP'         : (  9.6445e-01,  1.79e-02, -1. )
             , 'timeResSigmaSF'   : (  1.4670e+00,  5.83e-02, -1. )
             , 'wTagDelP0OS'      : (  1.0018e-02,  1.00e-03, -1. )
             , 'wTagDelP0SS'      : ( -1.5953e-02,  2.00e-03, -1. )
             , 'wTagDelP1OS'      : (  6.9817e-02,  1.00e-02, -1. )
             , 'wTagDelP1SS'      : (  1.5012e-02,  1.90e-02, -1. )
             , 'wTagP0OS'         : (  3.9051e-01,  9.16e-03, -1. )
             , 'wTagP0SS'         : (  4.4062e-01,  6.68e-03, -1. )
             , 'wTagP1OS'         : (  1.0310e+00,  5.89e-02, -1. )
             , 'wTagP1SS'         : (  9.4109e-01,  1.03e-01, -1. )
            }


# trigger selection strings
triggerSelStrings = dict(  noSelection    = ''
                         , HLT1Unbiased   = 'hlt1_unbiased_dec==1 && hlt2_biased==1'
                         , HLT1ExclBiased = 'hlt1_excl_biased_dec==1 && hlt2_biased==1'
                         , paper2012      = '(hlt1_excl_biased_dec==1 || hlt1_unbiased_dec==1) && hlt2_biased==1'
                         , timeEffFit     = '(hlt1_excl_biased_dec==1 || hlt1_unbiased_dec==1) && (hlt2_biased==1 || hlt2_unbiased==1)'
                         , unbiased       = 'hlt1_unbiased==1 && hlt2_unbiased==1'
                        )

# cut selection strings
cutSelStrings = dict(  noSelection = ''
                     , nominal2011 = 'sel == 1 && sel_cleantail == 1'\
                                     ' && time>0.3 && time<14. && sigmat<0.12'\
                                     ' && mass>5200. && mass<5550. && abs(mdau1-3090.)<60. && abs(mdau2-1020.)<30.'\
                                     ' && muplus_track_chi2ndof < 4. && muminus_track_chi2ndof < 4.'\
                                     ' && Kplus_track_chi2ndof < 4. && Kminus_track_chi2ndof < 4.'
                    )

# external constraint values dictionary
class ExtConstrValsDict( dict ) :
    def __init__( self, valsDict ) :
        for key, val in valsDict.iteritems() : self[key] = val

    def valsDict(self)           : return self.copy()
    def getVal( self, key )      : return self[key]
    def setVal( self, key, val ) : self[key] = val

    def getSetVal( self, key, val = None ) :
        if not key in self :
            assert val != None, 'ExtConstrValsDict.getSetVal(): no value provided'
            self[key] = val
        return self[key]

constrVals = dict(  DM      = (  17.768, 0.024  )
                  , beta    = ( -8.3e-3, 4.0e-3 )
                  , P0OS    = (  0.392,  0.008, 0.392 )
                  , P0SS    = (  0.350,  0.017, 0.350 )
                  , P1OS    = (  1.000,  0.023  )
                  , P1SS    = (  1.00,   0.16   )
                  , DelP0OS = (  0.0110, 0.0034 )
                  , DelP0SS = ( -0.019,  0.005  )
                  , DelP1OS = (  0.000,  0.001  )
                  , DelP1SS = (  0.00,   0.01   )
                 )
global extConstraintValues
extConstraintValues = ExtConstrValsDict( constrVals )
