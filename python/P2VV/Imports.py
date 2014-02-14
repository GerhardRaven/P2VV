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
            , 'phiCPAv'                   : ( 'phi_s^av',                    '$\\phi_\\text{s}^\\text{av}$'               )
            , '__phiCPAv__'               : ( 'phi_s^av (b)',                '$\\phi_\\text{s}^\\text{av}$ (b)'           )
            , 'phiCP_A0'                  : ( 'phi_s^0',                     '$\\phi_\\text{s}^0$'                        )
            , '__phiCP_A0__'              : ( 'phi_s^0 (b)',                 '$\\phi_\\text{s}^0$ (b)'                    )
            , 'phiCP_Apar'                : ( 'phi_s^para',                  '$\\phi_\\text{s}^\\parallel$'               )
            , '__phiCP_Apar__'            : ( 'phi_s^para (b)',              '$\\phi_\\text{s}^\\parallel$ (b)'           )
            , 'phiCP_Aperp'               : ( 'phi_s^perp',                  '$\\phi_\\text{s}^\\perp$'                   )
            , '__phiCP_Aperp__'           : ( 'phi_s^perp (b)',              '$\\phi_\\text{s}^\\perp$ (b)'               )
            , 'phiCP_AS'                  : ( 'phi_s^S',                     '$\\phi_\\text{s}^\\text{S}$'                )
            , '__phiCP_AS__'              : ( 'phi_s^S (b)',                 '$\\phi_\\text{s}^\\text{S}$ (b)'            )
            , 'phiCP_m'                   : ( 'phi_s^0',                     '$\\phi_\\text{s}^0$'                        )
            , '__phiCP_m__'               : ( 'phi_s^0 (b)',                 '$\\phi_\\text{s}^0$ (b)'                    )
            , 'phiCPRel_Apar'             : ( 'Delta phi_s^para',            '$\\Delta \\phi_\\text{s}^\\parallel$'       )
            , '__phiCPRel_Apar__'         : ( 'Delta phi_s^para (b)',        '$\\Delta \\phi_\\text{s}^\\parallel$ (b)'   )
            , 'phiCPRel_Aperp'            : ( 'Delta phi_s^perp',            '$\\Delta \\phi_\\text{s}^\\perp$'           )
            , '__phiCPRel_Aperp__'        : ( 'Delta phi_s^perp (b)',        '$\\Delta \\phi_\\text{s}^\\perp$ (b)'       )
            , 'phiCPRel_AperpApar'        : ( 'Delta phi_s^perp\'',          '$\\Delta \\phi_\\text{s}^\\perp\'$'         )
            , '__phiCPRel_AperpApar__'    : ( 'Delta phi_s^perp\' (b)',      '$\\Delta \\phi_\\text{s}^\\perp\'$ (b)'     )
            , 'phiCPRel_AS'               : ( 'Delta phi_s^S',               '$\\Delta \\phi_\\text{s}^\\text{S}$'        )
            , '__phiCPRel_AS__'           : ( 'Delta phi_s^S (b)',           '$\\Delta \\phi_\\text{s}^\\text{S}$ (b)'    )
            , 'lambdaCP'                  : ( '|lambda_s|',                  '$|\\lambda_\\text{s}|$'                     )
            , '__lambdaCP__'              : ( '|lambda_s| (b)',              '$|\\lambda_\\text{s}|$ (b)'                 )
            , 'rhoCP_A0'                  : ( '|lambda_s^0|',                '$|\\lambda_\\text{s}^0|$'                   )
            , 'rhoCP_Apar'                : ( '|lambda_s^para|',             '$|\\lambda_\\text{s}^\\parallel|$'          )
            , 'rhoCP_Aperp'               : ( '|lambda_s^perp|',             '$|\\lambda_\\text{s}^\\perp|$'              )
            , 'rhoCP_AS'                  : ( '|lambda_s^S|',                '$|\\lambda_\\text{s}^\\text{S}|$'           )
            , 'lambdaCP_A0'               : ( '|lambda_s^0|',                '$|\\lambda_\\text{s}^0|$'                   )
            , 'lambdaCP_Apar'             : ( '|lambda_s^para|',             '$|\\lambda_\\text{s}^\\parallel|$'          )
            , 'lambdaCP_Aperp'            : ( '|lambda_s^perp|',             '$|\\lambda_\\text{s}^\\perp|$'              )
            , 'lambdaCP_AS'               : ( '|lambda_s^S|',                '$|\\lambda_\\text{s}^\\text{S}|$'           )
            , 'CCP'                       : ( 'C_s',                         '$C_\\text{s}$'                              )
            , 'CCPAv'                     : ( 'C_s^av',                      '$C_\\text{s}^\\text{av}$'                   )
            , 'CCPAv_AS'                  : ( 'C_s^avS',                     '$C_\\text{s}^\\text{avS}$'                  )
            , 'CCPRel_Apar'               : ( 'Delta C_s^para',              '$\\Delta C_\\text{s}^\\parallel$'           )
            , 'CCPRel_Aperp'              : ( 'Delta C_s^perp',              '$\\Delta C_\\text{s}^\\perp$'               )
            , 'Gamma'                     : ( 'Gamma_s',                     '$\\Gamma_\\text{s}$'                        )
            , 'Gamma_p2011'               : ( 'Gamma_s - beta_2011',         '$\\Gamma_\\text{s} - \\beta_\\{2011}$'      )
            , 'Gamma_p2012'               : ( 'Gamma_s - beta_2012',         '$\\Gamma_\\text{s} - \\beta_\\{2012}$'      )
            , 'dGamma'                    : ( 'Delta Gamma_s',               '$\\Delta\\Gamma_\\text{s}$'                 )
            , '__dGamma__'                : ( 'Delta Gamma_s (b)',           '$\\Delta\\Gamma_\\text{s}$ (b)'             )
            , 'dM'                        : ( 'Delta m_s',                   '$\\Delta m_\\text{s}$'                      )
            , 'A0Mag2'                    : ( '|A_0|^2',                     '$|A_0|^2$'                                  )
            , 'AperpMag2'                 : ( '|A_perp|^2',                  '$|A_\\perp|^2$'                             )
            , 'avA02'                     : ( 'A^av_0^2',                    '${A^\\text{av}_0}^2$'                       )
            , 'avAperp2'                  : ( 'A^av_perp^2',                 '${A^\\text{av}_\\perp}^2$'                  )
            , 'delA02'                    : ( 'Delta A_0^2',                 '$\\Delta A_0^2$'                            )
            , 'delApar2'                  : ( 'Delta A_para^2',              '$\\Delta A_\\parallel^2$'                   )
            , 'delAperp2'                 : ( 'Delta A_perp^2',              '$\\Delta A_\\perp^2$'                       )
            , 'f_S'                       : ( 'F_S',                         '$F_\\text{S}$'                              )
            , 'f_S_bin0'                  : ( 'F_S_0',                       '${F_\\text{S}}_0$'                          )
            , 'f_S_bin1'                  : ( 'F_S_1',                       '${F_\\text{S}}_1$'                          )
            , 'f_S_bin2'                  : ( 'F_S_2',                       '${F_\\text{S}}_2$'                          )
            , 'f_S_bin3'                  : ( 'F_S_3',                       '${F_\\text{S}}_3$'                          )
            , 'f_S_bin4'                  : ( 'F_S_4',                       '${F_\\text{S}}_4$'                          )
            , 'f_S_bin5'                  : ( 'F_S_5',                       '${F_\\text{S}}_5$'                          )
            , 'avf_S'                     : ( 'F_S^av',                      '$F_\\text{S}$'                              )
            , 'avf_S_bin0'                : ( 'F_S^av_0',                    '${F_\\text{S}}_0$'                          )
            , 'avf_S_bin1'                : ( 'F_S^av_1',                    '${F_\\text{S}}_1$'                          )
            , 'avf_S_bin2'                : ( 'F_S^av_2',                    '${F_\\text{S}}_2$'                          )
            , 'avf_S_bin3'                : ( 'F_S^av_3',                    '${F_\\text{S}}_3$'                          )
            , 'avf_S_bin4'                : ( 'F_S^av_4',                    '${F_\\text{S}}_4$'                          )
            , 'avf_S_bin5'                : ( 'F_S^av_5',                    '${F_\\text{S}}_5$'                          )
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

# common parameter names dictionary:  'P2VV internal name' : 'common name'
commonParNames = {  'phiCP'                  : 'phis'
                  , '__phiCP__'              : 'phis'
                  , 'phiCPAv'                : 'avphis'
                  , '__phiCPAv__'            : 'avphis'
                  , 'phiCP_A0'               : 'phiszero'
                  , '__phiCP_A0__'           : 'phiszero'
                  , 'phiCP_Apar'             : 'phispara'
                  , '__phiCP_Apar__'         : 'phispara'
                  , 'phiCP_Aperp'            : 'phisperp'
                  , '__phiCP_Aperp__'        : 'phisperp'
                  , 'phiCP_AS'               : 'phisS'
                  , '__phiCP_AS__'           : 'phisS'
                  , 'phiCP_m'                : 'phiszero'
                  , '__phiCP_m__'            : 'phiszero'
                  , 'phiCPRel_Apar'          : 'delphispara'
                  , '__phiCPRel_Apar__'      : 'delphispara'
                  , 'phiCPRel_Aperp'         : 'delphisperp'
                  , '__phiCPRel_Aperp__'     : 'delphisperp'
                  , 'phiCPRel_AperpApar'     : 'delphisperpprime'
                  , '__phiCPRel_AperpApar__' : 'delphisperpprime'
                  , 'phiCPRel_AS'            : 'delphisS'
                  , '__phiCPRel_AS__'        : 'delphisS'
                  , 'lambdaCP'               : 'lambda'
                  , '__lambdaCP__'           : 'lambda'
                  , 'rhoCP_A0'               : 'lambdazero'
                  , 'rhoCP_Apar'             : 'lambdapara'
                  , 'rhoCP_Aperp'            : 'lambdaperp'
                  , 'rhoCP_AS'               : 'lambdaS'
                  , 'rhoCP_m'                : 'lambdazero'
                  , 'lambdaCP_A0'            : 'lambdazero'
                  , 'lambdaCP_Apar'          : 'lambdapara'
                  , 'lambdaCP_Aperp'         : 'lambdaperp'
                  , 'lambdaCP_AS'            : 'lambdaS'
                  , 'CCPAv'                  : 'avC'
                  , 'CCPAv_AS'               : 'avCS'
                  , 'CCPRel_Apar'            : 'delCpara'
                  , 'CCPRel_Aperp'           : 'delCperp'
                  , 'Gamma'                  : 'Gamma'
                  , 'Gamma_p2011'            : 'Gamma_2011'
                  , 'Gamma_p2012'            : 'Gamma_2012'
                  , 'dGamma'                 : 'DelGam'
                  , '__dGamma__'             : 'DelGam'
                  , 'dM'                     : 'Delm'
                  , 'A0Mag2'                 : 'AzeroSq'
                  , 'AperpMag2'              : 'AperpSq'
                  , 'avA02'                  : 'avAzeroSq'
                  , 'avAperp2'               : 'avAperpSq'
                  , 'f_S'                    : 'FS'
                  , 'f_S_bin0'               : 'FS1'
                  , 'f_S_bin1'               : 'FS2'
                  , 'f_S_bin2'               : 'FS3'
                  , 'f_S_bin3'               : 'FS4'
                  , 'f_S_bin4'               : 'FS5'
                  , 'f_S_bin5'               : 'FS6'
                  , 'avf_S'                  : 'avFS'
                  , 'avf_S_bin0'             : 'avFS1'
                  , 'avf_S_bin1'             : 'avFS2'
                  , 'avf_S_bin2'             : 'avFS3'
                  , 'avf_S_bin3'             : 'avFS4'
                  , 'avf_S_bin4'             : 'avFS5'
                  , 'avf_S_bin5'             : 'avFS6'
                  , 'AparPhase'              : 'delpara'
                  , 'AperpPhase'             : 'delperp'
                  , 'ASOddPhase'             : 'delS'
                  , 'ASOddPhase_bin0'        : 'delS1'
                  , 'ASOddPhase_bin1'        : 'delS2'
                  , 'ASOddPhase_bin2'        : 'delS3'
                  , 'ASOddPhase_bin3'        : 'delS4'
                  , 'ASOddPhase_bin4'        : 'delS5'
                  , 'ASOddPhase_bin5'        : 'delS6'
                  , 'betaTimeEff'            : 'beta'
                  , 'betaTimeEff_p2011'      : 'beta_2011'
                  , 'betaTimeEff_p2012'      : 'beta_2012'
                 }

# nominal values of physics parameters
parValues = {  'A0Mag2'            : (  5.2101e-01, 3.44e-03, -1 )
             , 'ASOddPhase_bin0'   : (  8.4476e-01, 1.99e-01, -1 )
             , 'ASOddPhase_bin1'   : (  2.1642e+00, 2.66e-01, -1 )
             , 'ASOddPhase_bin2'   : (  4.7383e-01, 2.12e-01, -1 )
             , 'ASOddPhase_bin3'   : ( -3.6993e-01, 1.90e-01, -1 )
             , 'ASOddPhase_bin4'   : ( -6.1136e-01, 1.58e-01, -1 )
             , 'ASOddPhase_bin5'   : ( -8.9845e-01, 1.40e-01, -1 )
             , 'AparPhase'         : (  3.2707e+00, 1.18e-01, -1 )
             , 'AperpMag2'         : (  2.5322e-01, 4.95e-03, -1 )
             , 'AperpPhase'        : (  3.1817e+00, 1.12e-01, -1 )
             , 'Gamma'             : (  6.6126e-01, 5.98e-03, -1 )
             , '__dGamma__'        : (  8.7297e-02, 9.14e-03, -1 )
             , '__phiCP__'         : (  7.9247e-02, 5.00e-02, -1 )
             , 'betaTimeEff_p2011' : ( -8.2753e-03, 4.00e-03, -1 )
             , 'betaTimeEff_p2012' : ( -1.3459e-02, 6.53e-03, -1 )
             , 'dM'                : (  1.77621e+01, 2.20e-02, -1 )
             , 'f_S_bin0'          : (  4.4288e-01, 5.48e-02, -1 )
             , 'f_S_bin1'          : (  6.0081e-02, 1.75e-02, -1 )
             , 'f_S_bin2'          : (  1.0347e-02, 6.82e-03, -1 )
             , 'f_S_bin3'          : (  8.9317e-03, 5.78e-03, -1 )
             , 'f_S_bin4'          : (  4.9146e-02, 1.56e-02, -1 )
             , 'f_S_bin5'          : (  1.8988e-01, 2.56e-02, -1 )
             , 'lambdaCP'          : (  9.6233e-01, 1.82e-02, -1 )
             , 'wTagP0OS'          : (  3.8769e-01, 4.25e-03, -1 )
             , 'wTagP0SS'          : (  4.5075e-01, 4.64e-03, -1 )
             , 'wTagP1OS'          : (  1.0113e+00, 3.44e-02, -1 )
             , 'wTagP1SS'          : (  9.5740e-01, 8.23e-02, -1 )
            }

parValuesCPVDecay = {  'ASOddPhase_bin0'        : (  8.6911e-01, 2.01e-01, -1 )
                     , 'ASOddPhase_bin1'        : (  2.15913e+00, 2.90e-01, -1 )
                     , 'ASOddPhase_bin2'        : (  4.9734e-01, 2.14e-01, -1 )
                     , 'ASOddPhase_bin3'        : ( -3.5862e-01, 1.95e-01, -1 )
                     , 'ASOddPhase_bin4'        : ( -6.0399e-01, 1.64e-01, -1 )
                     , 'ASOddPhase_bin5'        : ( -8.8556e-01, 1.47e-01, -1 )
                     , 'AparPhase'              : (  3.26381e+00, 1.15e-01, -1 )
                     , 'AperpPhase'             : (  3.16213e+00, 1.17e-01, -1 )
                     , 'CCPAv'                  : (  1.13358e-02, 3.81e-02, -1 )
                     , 'CCPAv_AS'               : (  4.66513e-02, 2.90e-02, -1 )
                     , 'CCPRel_Apar'            : ( -8.7215e-03, 1.20e-01, -1 )
                     , 'CCPRel_Aperp'           : (  6.0304e-02, 1.63e-01, -1 )
                     , 'Gamma'                  : (  6.61237e-01, 5.99e-03, -1 )
                     , '__dGamma__'             : (  8.7100e-02, 9.18e-03, -1 )
                     , '__phiCPAv__'            : (  1.34655e-01, 5.03e-02, -1 )
                     , '__phiCPRel_AS__'        : ( -2.05711e-01, 6.22e-02, -1 )
                     , '__phiCPRel_Apar__'      : ( -1.35727e-01, 4.14e-02, -1 )
                     , '__phiCPRel_AperpApar__' : (  9.9645e-02, 2.48e-02, -1 )
                     , 'avA02'                  : (  5.2095e-01, 3.45e-03, -1 )
                     , 'avAperp2'               : (  2.53186e-01, 4.98e-03, -1 )
                     , 'avf_S_bin0'             : (  4.41733e-01, 5.48e-02, -1 )
                     , 'avf_S_bin1'             : (  5.8722e-02, 1.78e-02, -1 )
                     , 'avf_S_bin2'             : (  1.03081e-02, 7.03e-03, -1 )
                     , 'avf_S_bin3'             : (  8.7243e-03, 5.68e-03, -1 )
                     , 'avf_S_bin4'             : (  4.8540e-02, 1.57e-02, -1 )
                     , 'avf_S_bin5'             : (  1.89062e-01, 2.57e-02, -1 )
                     , 'betaTimeEff_p2011'      : ( -8.2595e-03, 4.00e-03, -1 )
                     , 'betaTimeEff_p2012'      : ( -1.34396e-02, 6.53e-03, -1 )
                     , 'dM'                     : (  1.77589e+01, 2.23e-02, -1 )
                     , 'wTagP0OS'               : (  3.87495e-01, 4.25e-03, -1 )
                     , 'wTagP0SS'               : (  4.50871e-01, 4.65e-03, -1 )
                     , 'wTagP1OS'               : (  1.01095e+00, 3.44e-02, -1 )
                     , 'wTagP1SS'               : (  9.61067e-01, 8.26e-02, -1 )
                    }

parValues20131203 = {  'A0Mag2'            : (  5.2080e-01, 3.45e-03, -1 )
                     , 'ASOddPhase_bin0'   : (  8.2353e-01, 1.85e-01, -1 )
                     , 'ASOddPhase_bin1'   : (  2.2694e+00, 2.22e-01, -1 )
                     , 'ASOddPhase_bin2'   : (  4.2763e-01, 1.90e-01, -1 )
                     , 'ASOddPhase_bin3'   : ( -3.6125e-01, 1.86e-01, -1 )
                     , 'ASOddPhase_bin4'   : ( -6.2751e-01, 1.65e-01, -1 )
                     , 'ASOddPhase_bin5'   : ( -8.9047e-01, 1.38e-01, -1 )
                     , 'AparPhase'         : (  3.2492e+00, 1.29e-01, -1 )
                     , 'AperpMag2'         : (  2.5361e-01, 4.92e-03, -1 )
                     , 'AperpPhase'        : (  3.1662e+00, 1.17e-01, -1 )
                     , 'Gamma'             : (  6.6112e-01, 6.00e-03, -1 )
                     , '__dGamma__'        : (  8.6749e-02, 9.16e-03, -1 )
                     , '__phiCP__'         : (  7.0521e-02, 5.14e-02, -1 )
                     , 'betaTimeEff_p2011' : ( -8.2819e-03, 4.00e-03, -1 )
                     , 'betaTimeEff_p2012' : ( -1.3524e-02, 6.54e-03, -1 )
                     , 'dM'                : (  1.7762e+01, 2.21e-02, -1 )
                     , 'f_S_bin0'          : (  4.4391e-01, 5.47e-02, -1 )
                     , 'f_S_bin1'          : (  6.0500e-02, 1.79e-02, -1 )
                     , 'f_S_bin2'          : (  1.0335e-02, 6.74e-03, -1 )
                     , 'f_S_bin3'          : (  8.7773e-03, 5.82e-03, -1 )
                     , 'f_S_bin4'          : (  4.8669e-02, 1.61e-02, -1 )
                     , 'f_S_bin5'          : (  1.9547e-01, 2.58e-02, -1 )
                     , 'lambdaCP'          : (  9.6630e-01, 1.78e-02, -1 )
                     , 'wTagDelP0OS'       : (  1.3722e-02, 1.20e-03, -1 )
                     , 'wTagDelP0SS'       : ( -1.5989e-02, 2.00e-03, -1 )
                     , 'wTagDelP1OS'       : (  6.9705e-02, 1.20e-02, -1 )
                     , 'wTagDelP1SS'       : (  1.5208e-02, 1.90e-02, -1 )
                     , 'wTagP0OS'          : (  3.9250e-01, 9.06e-03, -1 )
                     , 'wTagP0SS'          : (  4.3908e-01, 6.69e-03, -1 )
                     , 'wTagP1OS'          : (  1.0296e+00, 5.69e-02, -1 )
                     , 'wTagP1SS'          : (  9.3231e-01, 1.03e-01, -1 )
                    }

parValues2011 = {  'phiCP'           : (  6.9616e-02, 9.09e-02, -1. )
                 , 'lambdaCP'        : (  9.4368e-01, 3.04e-02, -1. )
                 , 'Gamma'           : (  6.7132e-01, 4.80e-03, -1. )
                 , 'dGamma'          : (  1.0053e-01, 1.62e-02, -1. )
                 , 'dM'              : (  1.7667e+01, 7.68e-02, -1. )
                 , 'A0Mag2'          : (  5.2122e-01, 6.08e-03, -1. )
                 , 'AperpMag2'       : (  2.4848e-01, 8.69e-03, -1. )
                 , 'f_S_bin0'        : (  2.2751e-01, 7.74e-02, -1. )
                 , 'f_S_bin1'        : (  6.7228e-02, 2.86e-02, -1. )
                 , 'f_S_bin2'        : (  7.9881e-03, 1.16e-02, -1. )
                 , 'f_S_bin3'        : (  1.6367e-02, 1.08e-02, -1. )
                 , 'f_S_bin4'        : (  5.4997e-02, 2.64e-02, -1. )
                 , 'f_S_bin5'        : (  1.6757e-01, 4.24e-02, -1. )
                 , 'AparPhase'       : (  3.3097e+00, 1.54e-01, -1. )
                 , 'AperpPhase'      : (  3.0759e+00, 2.17e-01, -1. )
                 , 'ASOddPhase_bin0' : (  1.3140e+00, 6.81e-01, -1. )
                 , 'ASOddPhase_bin1' : (  7.6351e-01, 2.76e-01, -1. )
                 , 'ASOddPhase_bin2' : (  4.9424e-01, 4.60e-01, -1. )
                 , 'ASOddPhase_bin3' : ( -5.0928e-01, 2.55e-01, -1. )
                 , 'ASOddPhase_bin4' : ( -4.5728e-01, 2.08e-01, -1. )
                 , 'ASOddPhase_bin5' : ( -6.5000e-01, 1.98e-01, -1. )
                 , 'm_sig_frac'      : (  7.5966e-01, 3.49e-02, -1. )
                 , 'm_sig_mean'      : (  5.3682e+03, 4.84e-02, -1. )
                 , 'm_sig_sigma_1'   : (  6.0797e+00, 1.33e-01, -1. )
                 , 'm_sig_sigma_sf'  : (  2.0652e+00, 8.87e-02, -1. )
                 , 'm_bkg_exp_bin0'  : ( -1.7313e-03, 1.47e-04, -1. )
                 , 'm_bkg_exp_bin1'  : ( -1.7486e-03, 1.77e-04, -1. )
                 , 'm_bkg_exp_bin2'  : ( -1.7360e-03, 1.98e-04, -1. )
                 , 'm_bkg_exp_bin3'  : ( -1.2316e-03, 1.98e-04, -1. )
                 , 'm_bkg_exp_bin4'  : ( -1.5366e-03, 1.55e-04, -1. )
                 , 'm_bkg_exp_bin5'  : ( -1.5859e-03, 9.95e-05, -1. )
                }


# trigger selection strings
triggerSelStrings = dict(  noSelection        = ''
                         , HLT1Unbiased       = 'hlt1_unbiased_dec==1 && hlt2_biased==1'
                         , HLT1Unbiased_tos   = 'hlt1_unbiased==1 && hlt2_biased==1'
                         , HLT1ExclBiased     = 'hlt1_excl_biased_dec==1 && hlt2_biased==1'
                         , HLT1ExclBiased_tos = 'hlt1_excl_biased==1 && hlt2_biased==1'
                         , paper2012          = '(hlt1_excl_biased_dec==1 || hlt1_unbiased_dec==1) && hlt2_biased==1'
                         , paper2012_tos      = '(hlt1_excl_biased==1 || hlt1_unbiased==1) && hlt2_biased==1'
                         , timeEffFit         = '(hlt1_excl_biased_dec==1 || hlt1_unbiased_dec==1) && (hlt2_biased==1 || hlt2_unbiased==1)'
                         , timeEffFit_tos     = '(hlt1_excl_biased==1 || hlt1_unbiased==1) && (hlt2_biased==1 || hlt2_unbiased==1)'
                         , unbiased           = 'hlt1_unbiased==1 && hlt2_unbiased==1'
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
