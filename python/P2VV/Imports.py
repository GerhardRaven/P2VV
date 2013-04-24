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
            , 'timeResSigmaSF'            : ( 'sigma_t',                     '$\\sigma_t$ s.f.'                           )
            , 'wTagP0OS'                  : ( 'p_0 OS',                      '$p_0$ OS'                                   )
            , 'wTagP1OS'                  : ( 'p_1 OS',                      '$p_1$ OS'                                   )
            , 'wTagDelP0OS'               : ( 'Delta p_0 OS',                '$\\Delta p_0$ OS'                           )
            , 'wTagP0SS'                  : ( 'p_0 SS',                      '$p_0$ SS'                                   )
            , 'wTagP1SS'                  : ( 'p_1 SS',                      '$p_1$ SS'                                   )
            , 'wTagDelP0SS'               : ( 'Delta p_1 SS',                '$\\Delta p_0$ SS'                           )
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
parValues = {  'phiCP'           : (  9.6919e-02,  8.93e-02, -1.       )
             , 'lambdaCP'        : (  9.4355e-01,  3.13e-02, -1.       )
             , 'Gamma'           : (  6.7132e-01,  4.80e-03, -1.       )
             , 'dGamma'          : (  1.0122e-01,  1.62e-02, -1.       )
             , 'dM'              : (  1.7765e+01,  2.34e-02, -1.       )
             , 'A0Mag2'          : (  5.2127e-01,  6.07e-03, -1.       )
             , 'AperpMag2'       : (  2.4828e-01,  8.67e-03, -1.       )
             , 'f_S_bin0'        : (  2.2996e-01,  7.78e-02, -1.       )
             , 'f_S_bin1'        : (  6.8150e-02,  2.86e-02, -1.       )
             , 'f_S_bin2'        : (  8.1262e-03,  1.16e-02, -1.       )
             , 'f_S_bin3'        : (  1.5272e-02,  1.06e-02, -1.       )
             , 'f_S_bin4'        : (  5.6655e-02,  2.63e-02, -1.       )
             , 'f_S_bin5'        : (  1.6512e-01,  4.27e-02, -1.       )
             , 'AparPhase'       : (  3.3127e+00,  1.52e-01, -1.       )
             , 'AperpPhase'      : (  3.2382e+00,  1.86e-01, -1.       )
             , 'ASOddPhase_bin0' : (  1.2465e+00,  5.86e-01, -1.       )
             , 'ASOddPhase_bin1' : (  7.6076e-01,  2.71e-01, -1.       )
             , 'ASOddPhase_bin2' : (  5.0161e-01,  4.60e-01, -1.       )
             , 'ASOddPhase_bin3' : ( -5.3139e-01,  2.74e-01, -1.       )
             , 'ASOddPhase_bin4' : ( -4.5022e-01,  2.02e-01, -1.       )
             , 'ASOddPhase_bin5' : ( -6.5691e-01,  2.02e-01, -1.       )
             , 'timeResSigmaSF'  : (  1.4497,      5.48e-02, -1.       )
             , 'wTagP0OS'        : (  3.9098e-01,  7.69e-03, -1.       )
             , 'wTagP1OS'        : (  1.0003e+00,  2.29e-02, -1.       )
             , 'wTagP0SS'        : (  3.6110e-01,  1.59e-02, -1.       )
             , 'wTagP1SS'        : (  1.0204e+00,  1.57e-01, -1.       )
             , 'wTagDelP0OS'     : (  1.1008e-02,  3.40e-03, -1.       )
             , 'wTagDelP0SS'     : ( -1.9036e-02,  5.00e-03, -1.       )
             , 'm_sig_mean'      : (  5.36822e+03, 4.84e-02, -1.       )
             , 'm_sig_frac'      : (  7.5966e-01,  3.49e-02, -1.       )
             , 'm_sig_sigma_1'   : (  6.0797e+00,  1.33e-01, -1.       )
             , 'm_sig_sigma_sf'  : (  2.0652e+00,  8.87e-02, -1.       )
             , 'm_sig_widthPar0' : (  1.7415e+01,  1.00,     -1.       )
             , 'm_sig_widthPar1' : ( -5.9840e+01,  1.00,     -1.       )
             , 'm_sig_widthPar2' : (  4.2391e+01,  1.00,     -1.       )
             , 'm_bkg_exp_bin0'  : ( -1.7313e-03,  1.47e-04, -1.       )
             , 'm_bkg_exp_bin1'  : ( -1.7486e-03,  1.77e-04, -1.       )
             , 'm_bkg_exp_bin2'  : ( -1.7360e-03,  1.98e-04, -1.       )
             , 'm_bkg_exp_bin3'  : ( -1.2316e-03,  1.98e-04, -1.       )
             , 'm_bkg_exp_bin4'  : ( -1.5366e-03,  1.55e-04, -1.       )
             , 'm_bkg_exp_bin5'  : ( -1.5859e-03,  9.95e-05, -1.       )
            }

parValues2011 = {  'phiCP'           : (  6.9616e-02,  9.09e-02, -1.       )
                 , 'lambdaCP'        : (  9.4368e-01,  3.04e-02, -1.       )
                 , 'Gamma'           : (  6.7132e-01,  4.80e-03, -1.       )
                 , 'dGamma'          : (  1.0053e-01,  1.62e-02, -1.       )
                 , 'dM'              : (  1.7667e+01,  7.68e-02, -1.       )
                 , 'A0Mag2'          : (  5.2122e-01,  6.08e-03, -1.       )
                 , 'AperpMag2'       : (  2.4848e-01,  8.69e-03, -1.       )
                 , 'f_S_bin0'        : (  2.2751e-01,  7.74e-02, -1.       )
                 , 'f_S_bin1'        : (  6.7228e-02,  2.86e-02, -1.       )
                 , 'f_S_bin2'        : (  7.9881e-03,  1.16e-02, -1.       )
                 , 'f_S_bin3'        : (  1.6367e-02,  1.08e-02, -1.       )
                 , 'f_S_bin4'        : (  5.4997e-02,  2.64e-02, -1.       )
                 , 'f_S_bin5'        : (  1.6757e-01,  4.24e-02, -1.       )
                 , 'AparPhase'       : (  3.3097e+00,  1.54e-01, -1.       )
                 , 'AperpPhase'      : (  3.0759e+00,  2.17e-01, -1.       )
                 , 'ASOddPhase_bin0' : (  1.3140e+00,  6.81e-01, -1.       )
                 , 'ASOddPhase_bin1' : (  7.6351e-01,  2.76e-01, -1.       )
                 , 'ASOddPhase_bin2' : (  4.9424e-01,  4.60e-01, -1.       )
                 , 'ASOddPhase_bin3' : ( -5.0928e-01,  2.55e-01, -1.       )
                 , 'ASOddPhase_bin4' : ( -4.5728e-01,  2.08e-01, -1.       )
                 , 'ASOddPhase_bin5' : ( -6.5000e-01,  1.98e-01, -1.       )
                 , 'wTagP0OS'        : (  3.9065e-01,  7.69e-03, -1.       )
                 , 'wTagP0SS'        : (  3.6119e-01,  1.59e-02, -1.       )
                 , 'wTagP1OS'        : (  1.0001e+00,  2.29e-02, -1.       )
                 , 'wTagP1SS'        : (  1.0184e+00,  1.57e-01, -1.       )
                 , 'wTagDelP0OS'     : (  1.1002e-02,  3.40e-03, -1.       )
                 , 'wTagDelP0SS'     : ( -1.9038e-02,  5.00e-03, -1.       )
                 , 'm_sig_frac'      : (  7.5966e-01,  3.49e-02, -1.       )
                 , 'm_sig_mean'      : (  5.3682e+03,  4.84e-02, -1.       )
                 , 'm_sig_sigma_1'   : (  6.0797e+00,  1.33e-01, -1.       )
                 , 'm_sig_sigma_sf'  : (  2.0652e+00,  8.87e-02, -1.       )
                 , 'm_bkg_exp_bin0'  : ( -1.7313e-03,  1.47e-04, -1.       )
                 , 'm_bkg_exp_bin1'  : ( -1.7486e-03,  1.77e-04, -1.       )
                 , 'm_bkg_exp_bin2'  : ( -1.7360e-03,  1.98e-04, -1.       )
                 , 'm_bkg_exp_bin3'  : ( -1.2316e-03,  1.98e-04, -1.       )
                 , 'm_bkg_exp_bin4'  : ( -1.5366e-03,  1.55e-04, -1.       )
                 , 'm_bkg_exp_bin5'  : ( -1.5859e-03,  9.95e-05, -1.       )
                }

# values of physics parameters with HLT1 unbiased data
parValuesUnbiased = {  'phiCP'           : (  8.5502e-02,  9.19e-02, -1.       )
                     , 'lambdaCP'        : (  9.5703e-01,  3.27e-02, -1.       )
                     , 'Gamma'           : (  6.7176e-01,  5.18e-03, -1.       )
                     , 'dGamma'          : (  1.1044e-01,  1.74e-02, -1.       )
                     , 'dM'              : (  1.7765e+01,  2.37e-02, -1.       )
                     , 'A0Mag2'          : (  5.2182e-01,  6.53e-03, -1.       )
                     , 'AperpMag2'       : (  2.4959e-01,  9.28e-03, -1.       )
                     , 'f_S_bin0'        : (  1.9923e-01,  7.47e-02, -1.       )
                     , 'f_S_bin1'        : (  7.7526e-02,  3.18e-02, -1.       )
                     , 'f_S_bin2'        : (  1.5327e-02,  1.52e-02, -1.       )
                     , 'f_S_bin3'        : (  1.2846e-02,  1.05e-02, -1.       )
                     , 'f_S_bin4'        : (  5.0058e-02,  2.79e-02, -1.       )
                     , 'f_S_bin5'        : (  1.5209e-01,  4.51e-02, -1.       )
                     , 'AparPhase'       : (  3.3111e+00,  1.62e-01, -1.       )
                     , 'AperpPhase'      : (  3.2840e+00,  1.94e-01, -1.       )
                     , 'ASOddPhase_bin0' : (  1.3590e+00,  5.19e-01, -1.       )
                     , 'ASOddPhase_bin1' : (  7.3779e-01,  2.66e-01, -1.       )
                     , 'ASOddPhase_bin2' : (  3.2396e-01,  2.58e-01, -1.       )
                     , 'ASOddPhase_bin3' : ( -5.2282e-01,  3.14e-01, -1.       )
                     , 'ASOddPhase_bin4' : ( -4.9460e-01,  2.42e-01, -1.       )
                     , 'ASOddPhase_bin5' : ( -6.8709e-01,  2.28e-01, -1.       )
                     , 'timeResSigmaSF'  : (  1.4497,      5.48e-02, -1.       )
                     , 'wTagP0OS'        : (  3.9062e-01,  7.72e-03, -1.       )
                     , 'wTagP0SS'        : (  3.6251e-01,  1.60e-02, -1.       )
                     , 'wTagP1OS'        : (  1.0012e+00,  2.30e-02, -1.       )
                     , 'wTagP1SS'        : (  1.0142e+00,  1.57e-01, -1.       )
                     , 'wTagDelP0OS'     : (  1.0992e-02,  3.40e-03, -1.       )
                     , 'wTagDelP0SS'     : ( -1.9022e-02,  5.00e-03, -1.       )
                     , 'm_sig_mean'      : (  5.368223e+03,  5.26e-02, -1.       )
                     , 'm_sig_frac'      : (  7.7116e-01,  3.73e-02, -1.       )
                     , 'm_sig_sigma_1'   : (  6.1485e+00,  1.43e-01, -1.       )
                     , 'm_sig_sigma_sf'  : (  2.0689e+00,  1.03e-01, -1.       )
                     , 'm_bkg_exp_bin0'  : ( -1.7821e-03,  1.58e-04, -1.       )
                     , 'm_bkg_exp_bin1'  : ( -1.7280e-03,  1.90e-04, -1.       )
                     , 'm_bkg_exp_bin2'  : ( -1.8612e-03,  2.12e-04, -1.       )
                     , 'm_bkg_exp_bin3'  : ( -1.2156e-03,  2.11e-04, -1.       )
                     , 'm_bkg_exp_bin4'  : ( -1.5947e-03,  1.66e-04, -1.       )
                     , 'm_bkg_exp_bin5'  : ( -1.6132e-03,  1.06e-04, -1.       )
                    }

# values of physics parameters with HLT1 unbiased magnet down data
parValuesUnbMagDown = {  'phiCP'           : (  7.3384e-02,  1.23e-01, -1.       )
                       , 'lambdaCP'        : (  9.3791e-01,  8.80e-02, -1.       )
                       , 'Gamma'           : (  6.6573e-01,  6.81e-03, -1.       )
                       , 'dGamma'          : (  1.1965e-01,  2.32e-02, -1.       )
                       , 'dM'              : (  1.7766e+01,  2.40e-02, -1.       )
                       , 'A0Mag2'          : (  5.2588e-01,  8.57e-03, -1.       )
                       , 'AperpMag2'       : (  2.4172e-01,  1.21e-02, -1.       )
                       , 'f_S_bin0'        : (  1.6508e-01,  8.46e-02, -1.       )
                       , 'f_S_bin1'        : (  8.6053e-02,  4.83e-02, -1.       )
                       , 'f_S_bin2'        : (  9.2474e-03,  1.45e-02, -1.       )
                       , 'f_S_bin3'        : (  1.5459e-02,  1.77e-02, -1.       )
                       , 'f_S_bin4'        : (  3.5602e-02,  3.43e-02, -1.       )
                       , 'f_S_bin5'        : (  1.3670e-01,  7.61e-02, -1.       )
                       , 'AparPhase'       : (  3.3419e+00,  1.65e-01, -1.       )
                       , 'AperpPhase'      : (  3.3440e+00,  3.02e-01, -1.       )
                       , 'ASOddPhase_bin0' : (  1.5833e+00,  5.86e-01, -1.       )
                       , 'ASOddPhase_bin1' : (  7.3696e-01,  4.06e-01, -1.       )
                       , 'ASOddPhase_bin2' : (  1.0790e+00,  1.46e+00, -1.       )
                       , 'ASOddPhase_bin3' : ( -5.2537e-01,  4.38e-01, -1.       )
                       , 'ASOddPhase_bin4' : ( -5.6514e-01,  4.19e-01, -1.       )
                       , 'ASOddPhase_bin5' : ( -8.0370e-01,  4.88e-01, -1.       )
                       , 'wTagP0OS'        : (  3.9249e-01,  7.86e-03, -1.       )
                       , 'wTagP0SS'        : (  3.5899e-01,  1.64e-02, -1.       )
                       , 'wTagP1OS'        : (  1.0002e+00,  2.30e-02, -1.       )
                       , 'wTagP1SS'        : (  1.0168e+00,  1.58e-01, -1.       )
                       , 'wTagDelP0OS'     : (  1.0996e-02,  3.40e-03, -1.       )
                       , 'wTagDelP0SS'     : ( -1.9046e-02,  5.00e-03, -1.       )
                       , 'm_sig_mean'      : (  5.368161e+03,  7.00e-02, -1.       )
                       , 'm_sig_frac'      : (  8.3960e-01,  3.52e-02, -1.       )
                       , 'm_sig_sigma_1'   : (  6.4803e+00,  1.57e-01, -1.       )
                       , 'm_sig_sigma_sf'  : (  2.2898e+00,  2.15e-01, -1.       )
                       , 'm_bkg_exp_bin0'  : ( -1.8396e-03,  2.04e-04, -1.       )
                       , 'm_bkg_exp_bin1'  : ( -1.5986e-03,  2.46e-04, -1.       )
                       , 'm_bkg_exp_bin2'  : ( -2.0421e-03,  2.77e-04, -1.       )
                       , 'm_bkg_exp_bin3'  : ( -1.1750e-03,  2.77e-04, -1.       )
                       , 'm_bkg_exp_bin4'  : ( -1.4716e-03,  2.16e-04, -1.       )
                       , 'm_bkg_exp_bin5'  : ( -1.7165e-03,  1.38e-04, -1.       )
                      }

# values of physics parameters with phase space MC angular acceptance
parValuesPHSPAcc = {  'phiCP'           : (  7.0691e-02,  9.08e-02, -1.       )
                    , 'lambdaCP'        : (  9.4518e-01,  2.97e-02, -1.       )
                    , 'Gamma'           : (  6.7113e-01,  4.83e-03, -1.       )
                    , 'dGamma'          : (  1.0038e-01,  1.62e-02, -1.       )
                    , 'dM'              : (  1.7664e+01,  7.70e-02, -1.       )
                    , 'A0Mag2'          : (  5.2531e-01,  6.07e-03, -1.       )
                    , 'AperpMag2'       : (  2.4478e-01,  8.80e-03, -1.       )
                    , 'f_S_bin0'        : (  2.3063e-01, -0.073649,  0.080953 )
                    , 'f_S_bin1'        : (  6.9315e-02, -0.027325,  0.029634 )
                    , 'f_S_bin2'        : (  8.7124e-03, -0.007704,  0.014229 )
                    , 'f_S_bin3'        : (  1.6198e-02, -0.009186,  0.012207 )
                    , 'f_S_bin4'        : (  5.6242e-02, -0.025579,  0.027261 )
                    , 'f_S_bin5'        : (  1.6906e-01, -0.042141,  0.042779 )
                    , 'AparPhase'       : (  3.4047e+00, -0.153221,  0.104434 )
                    , 'AperpPhase'      : (  3.1016e+00,  2.15e-01, -1.       )
                    , 'ASOddPhase_bin0' : (  1.2889e+00, -0.471392,  0.758919 )
                    , 'ASOddPhase_bin1' : (  7.5367e-01, -0.225195,  0.36704  )
                    , 'ASOddPhase_bin2' : (  4.7920e-01, -0.288453,  1.36977  )
                    , 'ASOddPhase_bin3' : ( -5.2024e-01, -0.36217,   0.218655 )
                    , 'ASOddPhase_bin4' : ( -4.5437e-01, -0.25763,   0.183822 )
                    , 'ASOddPhase_bin5' : ( -6.5442e-01, -0.225561,  0.184025 )
                    , 'timeResSigmaSF'  : (  1.4458e+00,  4.53e-02, -1.       )
                    , 'wTagP0OS'        : (  3.9056e-01,  7.69e-03, -1.       )
                    , 'wTagP1OS'        : (  1.0001e+00,  2.29e-02, -1.       )
                    , 'wTagDelP0OS'     : (  1.1008e-02,  3.40e-03, -1.       )
                    , 'wTagP0SS'        : (  3.6112e-01,  1.59e-02, -1.       )
                    , 'wTagP1SS'        : (  1.0189e+00,  1.57e-01, -1.       )
                    , 'wTagDelP0SS'     : ( -1.9023e-02,  5.00e-03, -1.       )
                   }

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