from math import sqrt
plotsFilePath = 'timeAcc.ps'
histsFilePath = 'timeAcc.root'
histsType = 'fit' # 'hists' # 'fit'
dataPath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
nTupleFilePath = 'fitNTuple_peakBkg_2011_2012_Reco14_TOS_20140215.root'
#nTupleFilePath = 'DEC/fitNTuple_peakBkg_2011_2012_Reco14_20140116.root'
nTupleName = 'DecayTree'
period2011 = '' # 'summer'
applyExclBScale = True
hlt1TOS = True
incHLT2UB = False
hlt1UnbEff = 0.7
plotMinMax = dict( UB = ( 0.85, 1.0 ), exclB = ( 0.3, 1.5 ) if applyExclBScale else ( 0.05, 0.22 ) )

#binBounds = [ 0.3, 0.5, 1.0, 1.7, 2.7, 14.0 ]
#binBounds = [ 0.3, 0.46, 0.63, 0.83, 1.05, 1.32, 1.65, 2.07, 2.67, 3.7, 14. ]
binBounds = [  0.30000, 0.33726, 0.37550, 0.41475, 0.45508, 0.49654, 0.53920, 0.58314, 0.62843, 0.67516, 0.72342, 0.77332, 0.82497, 0.87849
             , 0.93404, 0.99177, 1.05185, 1.11448, 1.17991, 1.24837, 1.32018, 1.39567, 1.47524, 1.55936, 1.64858, 1.74356, 1.84509, 1.95415
             , 2.07194, 2.19998, 2.34022, 2.49526, 2.66856, 2.86502, 3.09180, 3.36000, 3.68820, 4.11123, 4.70718, 5.72483, 14.00000
            ]

histNames = dict(  UB_2011    = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
                 , UB_2012    = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1DiMuon_Hlt2DiMuonDetached'
                 , exclB_2011 = 'Bs_HltPropertimeAcceptance_Data_2011_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                 , exclB_2012 = 'Bs_HltPropertimeAcceptance_Data_2012_40bins_Hlt1TrackAndTrackMuonExcl_Hlt2DiMuonDetached'
                )
accHists = { }
accHists = {  '2011' : dict(  file  = dataPath + 'Bs_HltPropertimeAcceptance_Data_2011_40bins_TOS.root'
                            , UB    = histNames['UB_2011']
                            , exclB = histNames['exclB_2011']
                           )
            , '2012' : dict(  file  = dataPath + 'Bs_HltPropertimeAcceptance_Data_2012_40bins_TOS.root'
                            , UB    = histNames['UB_2012']
                            , exclB = histNames['exclB_2012']
                           )
           }

# 5 bins, uniform UB (2011: 0.4813; 2012: 0.4388)
#fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.0692, 0.0036 ), ( 0.1056, 0.0031 ), ( 0.1328, 0.0036 ), ( 0.1440, 0.0041 )
#                                             , ( 0.1519, 0.0043 )
#                                            ]
#                , ('2012', 'hlt1_exclB' ) : [  ( 0.0781, 0.0027 ), ( 0.1429, 0.0026 ), ( 0.1650, 0.0028 ), ( 0.1814, 0.0032 )
#                                             , ( 0.2054, 0.0035 )
#                                            ]
#                , ('2011', 'hlt2_B' )     : [  ( 0.9648, 0.0077 ), ( 0.9832, 0.0027 ), ( 0.9647, 0.0033 ), ( 0.9581, 0.0040 )
#                                             , ( 0.9665, 0.0037 )
#                                            ]
#                , ('2012', 'hlt2_B' )     : [  ( 0.9715, 0.0050 ), ( 0.9862, 0.0019 ), ( 0.9796, 0.0020 ), ( 0.9813, 0.0020 )
#                                             , ( 0.9753, 0.0023 )
#                                            ]
#               }
# 10 bins, uniform UB (2011: 0.4813; 2012: 0.4387)
#fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.0648, 0.0039 ), ( 0.0878, 0.0043 ), ( 0.1077, 0.0046 ), ( 0.1230, 0.0050 )
#                                             , ( 0.1350, 0.0052 ), ( 0.1298, 0.0050 ), ( 0.1448, 0.0055 ), ( 0.1421, 0.0053 )
#                                             , ( 0.1578, 0.0057 ), ( 0.1456, 0.0055 )
#                                            ]
#                , ('2012', 'hlt1_exclB' ) : [  ( 0.0768, 0.0030 ), ( 0.1133, 0.0035 ), ( 0.1449, 0.0037 ), ( 0.1567, 0.0039 )
#                                             , ( 0.1636, 0.0039 ), ( 0.1699, 0.0040 ), ( 0.1728, 0.0041 ), ( 0.1866, 0.0043 )
#                                             , ( 0.2058, 0.0045 ), ( 0.2057, 0.0047 )
#                                            ]
#                , ('2011', 'hlt2_B' )     : [  ( 0.954,  0.010  ), ( 0.9874, 0.0041 ), ( 0.9868, 0.0037 ), ( 0.9778, 0.0046 )
#                                             , ( 0.9672, 0.0048 ), ( 0.9592, 0.0054 ), ( 0.9509, 0.0060 ), ( 0.9629, 0.0051 )
#                                             , ( 0.9642, 0.0053 ), ( 0.9696, 0.0049 )
#                                            ]
#                , ('2012', 'hlt2_B' )     : [  ( 0.9700, 0.0059 ), ( 0.9802, 0.0037 ), ( 0.9854, 0.0030 ), ( 0.9908, 0.0025 )
#                                             , ( 0.9784, 0.0031 ), ( 0.9792, 0.0030 ), ( 0.9838, 0.0026 ), ( 0.9798, 0.0029 )
#                                             , ( 0.9779, 0.0029 ), ( 0.9712, 0.0036 )
#                                            ]
#               }
# 10 bins, free UB
#fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.0648, 0.0040 ), ( 0.0878, 0.0043 ), ( 0.1077, 0.0046 ), ( 0.1230, 0.0050 )
#                                             , ( 0.1350, 0.0052 ), ( 0.1298, 0.0050 ), ( 0.1447, 0.0055 ), ( 0.1420, 0.0053 )
#                                             , ( 0.1579, 0.0057 ), ( 0.1456, 0.0055 )
#                                            ]
#                , ('2012', 'hlt1_exclB' ) : [  ( 0.0768, 0.0030 ), ( 0.1133, 0.0035 ), ( 0.1449, 0.0037 ), ( 0.1568, 0.0039 )
#                                             , ( 0.1636, 0.0039 ), ( 0.1699, 0.0040 ), ( 0.1728, 0.0041 ), ( 0.1866, 0.0043 )
#                                             , ( 0.2058, 0.0045 ), ( 0.2057, 0.0047 )
#                                            ]
#                , ('2011', 'hlt2_B' )     : [  ( 0.954,  0.010  ), ( 0.9871, 0.0042 ), ( 0.9869, 0.0037 ), ( 0.9780, 0.0045 )
#                                             , ( 0.9667, 0.0049 ), ( 0.9597, 0.0053 ), ( 0.9510, 0.0060 ), ( 0.9638, 0.0050 )
#                                             , ( 0.9636, 0.0054 ), ( 0.9692, 0.0050 )
#                                            ]
#                , ('2012', 'hlt2_B' )     : [  ( 0.9700, 0.0059 ), ( 0.9797, 0.0037 ), ( 0.9853, 0.0030 ), ( 0.9908, 0.0025 )
#                                             , ( 0.9785, 0.0031 ), ( 0.9789, 0.0031 ), ( 0.9840, 0.0026 ), ( 0.9796, 0.0030 )
#                                             , ( 0.9783, 0.0029 ), ( 0.9716, 0.0036 )
#                                            ]
#               }
# 40 bins, uniform UB (2011: 0.4814; 2012: 0.4390)
fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.0635, 0.0079 ), ( 0.0563, 0.0077 ), ( 0.0854, 0.0092 ), ( 0.0572, 0.0070 )
                                             , ( 0.0810, 0.0080 ), ( 0.0837, 0.0080 ), ( 0.0876, 0.0082 ), ( 0.1004, 0.0086 )
                                             , ( 0.0885, 0.0081 ), ( 0.1184, 0.0093 ), ( 0.1039, 0.0086 ), ( 0.1232, 0.0095 )
                                             , ( 0.1075, 0.0089 ), ( 0.1255, 0.0095 ), ( 0.1305, 0.0097 ), ( 0.1276, 0.0095 )
                                             , ( 0.148,  0.010  ), ( 0.1267, 0.0097 ), ( 0.1171, 0.0094 ), ( 0.153,  0.011  )
                                             , ( 0.148,  0.010  ), ( 0.1129, 0.0091 ), ( 0.1276, 0.0096 ), ( 0.1312, 0.0097 )
                                             , ( 0.141,  0.010  ), ( 0.135,  0.010  ), ( 0.162,  0.011  ), ( 0.140,  0.010  )
                                             , ( 0.159,  0.011  ), ( 0.145,  0.010  ), ( 0.1274, 0.0095 ), ( 0.140,  0.010  )
                                             , ( 0.162,  0.011  ), ( 0.147,  0.010  ), ( 0.172,  0.011  ), ( 0.155,  0.011  )
                                             , ( 0.159,  0.011  ), ( 0.145,  0.010  ), ( 0.134,  0.010  ), ( 0.149,  0.010  )
                                            ]
                , ('2012', 'hlt1_exclB' ) : [  ( 0.0533, 0.0053 ), ( 0.0694, 0.0057 ), ( 0.0914, 0.0064 ), ( 0.0910, 0.0062 )
                                             , ( 0.0811, 0.0057 ), ( 0.1118, 0.0066 ), ( 0.1204, 0.0068 ), ( 0.1363, 0.0072 )
                                             , ( 0.1284, 0.0069 ), ( 0.1380, 0.0071 ), ( 0.1372, 0.0070 ), ( 0.1758, 0.0080 )
                                             , ( 0.1570, 0.0073 ), ( 0.1408, 0.0070 ), ( 0.1789, 0.0080 ), ( 0.1514, 0.0072 )
                                             , ( 0.1483, 0.0072 ), ( 0.1611, 0.0074 ), ( 0.1823, 0.0081 ), ( 0.1634, 0.0076 )
                                             , ( 0.1677, 0.0077 ), ( 0.1522, 0.0073 ), ( 0.1745, 0.0078 ), ( 0.1855, 0.0081 )
                                             , ( 0.1764, 0.0078 ), ( 0.1777, 0.0078 ), ( 0.1712, 0.0077 ), ( 0.1662, 0.0076 )
                                             , ( 0.1822, 0.0080 ), ( 0.1857, 0.0081 ), ( 0.2021, 0.0087 ), ( 0.1769, 0.0080 )
                                             , ( 0.2073, 0.0086 ), ( 0.2095, 0.0088 ), ( 0.2182, 0.0089 ), ( 0.1887, 0.0082 )
                                             , ( 0.1909, 0.0084 ), ( 0.1922, 0.0084 ), ( 0.2178, 0.0091 ), ( 0.2245, 0.0095 )
                                            ]
                , ('2011', 'hlt2_B' )     : [  ( 0.947,  0.020  ), ( 0.905,  0.024  ), ( 0.939,  0.036  ), ( 0.978,  0.012  )
                                             , ( 0.9911, 0.0070 ), ( 0.9899, 0.0079 ), ( 0.975,  0.010  ), ( 0.998,  0.012  )
                                             , ( 0.9801, 0.0089 ), ( 0.9954, 0.0055 ), ( 0.9892, 0.0073 ), ( 0.9780, 0.0086 )
                                             , ( 0.9857, 0.0090 ), ( 0.9763, 0.0089 ), ( 0.965,  0.010  ), ( 0.9881, 0.0069 )
                                             , ( 0.967,  0.010  ), ( 0.961,  0.011  ), ( 0.9669, 0.0093 ), ( 0.9718, 0.0089 )
                                             , ( 0.9744, 0.0087 ), ( 0.958,  0.012  ), ( 0.955,  0.011  ), ( 0.948,  0.011  )
                                             , ( 0.960,  0.011  ), ( 0.941,  0.013  ), ( 0.938,  0.013  ), ( 0.964,  0.011  )
                                             , ( 0.959,  0.011  ), ( 0.9756, 0.0081 ), ( 0.9685, 0.0095 ), ( 0.952,  0.011  )
                                             , ( 0.953,  0.013  ), ( 0.9682, 0.0097 ), ( 0.960,  0.010  ), ( 0.969,  0.010  )
                                             , ( 0.963,  0.011  ), ( 0.954,  0.011  ), ( 0.9776, 0.0083 ), ( 0.9837, 0.0083 )
                                            ]
                , ('2012', 'hlt2_B' )     : [  ( 0.942 , 0.015  ), ( 0.968 , 0.011  ), ( 0.970 , 0.010  ), ( 0.990 , 0.017  )
                                             , ( 0.9839, 0.0066 ), ( 0.9821, 0.0079 ), ( 0.9819, 0.0063 ), ( 0.9784, 0.0074 )
                                             , ( 0.9903, 0.0070 ), ( 0.9796, 0.0065 ), ( 0.9885, 0.0052 ), ( 0.9818, 0.0063 )
                                             , ( 0.9928, 0.0041 ), ( 0.9886, 0.0052 ), ( 0.9937, 0.0056 ), ( 0.9882, 0.0050 )
                                             , ( 0.9836, 0.0057 ), ( 0.9813, 0.0056 ), ( 0.9785, 0.0062 ), ( 0.9696, 0.0069 )
                                             , ( 0.9678, 0.0069 ), ( 0.9879, 0.0052 ), ( 0.9785, 0.0060 ), ( 0.9823, 0.0061 )
                                             , ( 0.9803, 0.0057 ), ( 0.9841, 0.0050 ), ( 0.9880, 0.0048 ), ( 0.9841, 0.0049 )
                                             , ( 0.9699, 0.0065 ), ( 0.9876, 0.0052 ), ( 0.9797, 0.0065 ), ( 0.9809, 0.0053 )
                                             , ( 0.9794, 0.0055 ), ( 0.9814, 0.0058 ), ( 0.9694, 0.0064 ), ( 0.9821, 0.0053 )
                                             , ( 0.9732, 0.0068 ), ( 0.9743, 0.0066 ), ( 0.9774, 0.0065 ), ( 0.9613, 0.0084 )
                                            ]
               }


fitAcc = { }
fitAcc = {  ('2011',    'UB') : [ eff for eff in fitAccLevels[ ( '2011', 'hlt2_B' ) ] ]
          , ('2011', 'exclB') : [ ( eff1[0] * eff2[0], sqrt( eff2[0]**2 * eff1[1]**2 + eff1[0]**2 * eff2[1]**2 ) )\
                                  for eff1, eff2 in zip( fitAccLevels[ ( '2011', 'hlt1_exclB' ) ], fitAccLevels[ ( '2011', 'hlt2_B' ) ] )
                                ]
          , ('2012',  'UB')   : [ eff for eff in fitAccLevels[ ( '2012', 'hlt2_B' ) ] ]
          , ('2012', 'exclB') : [ ( eff1[0] * eff2[0], sqrt( eff2[0]**2 * eff1[1]**2 + eff1[0]**2 * eff2[1]**2 ) )\
                                  for eff1, eff2 in zip( fitAccLevels[ ( '2012', 'hlt1_exclB' ) ], fitAccLevels[ ( '2012', 'hlt2_B' ) ] )
                                ]
         }
#fitAcc = {  ('2011',    'UB') : [  ( 0.475428,   0.       ), ( 4.1717e-01, 2.68e-02 ), ( 4.3387e-01, 2.76e-02 ), ( 4.9658e-01, 3.05e-02 )
#                                 , ( 5.3039e-01, 3.21e-02 ), ( 5.0833e-01, 3.11e-02 ), ( 4.9299e-01, 3.04e-02 ), ( 5.0210e-01, 3.08e-02 )
#                                 , ( 5.0438e-01, 3.09e-02 ), ( 5.0947e-01, 3.11e-02 ), ( 5.1237e-01, 3.13e-02 ), ( 4.8296e-01, 2.99e-02 )
#                                 , ( 4.7940e-01, 2.97e-02 ), ( 4.9937e-01, 3.07e-02 ), ( 4.8477e-01, 3.00e-02 ), ( 5.0051e-01, 3.07e-02 )
#                                 , ( 4.7051e-01, 2.93e-02 ), ( 4.9831e-01, 3.06e-02 ), ( 5.3250e-01, 3.22e-02 ), ( 4.8741e-01, 3.01e-02 )
#                                 , ( 5.2250e-01, 3.18e-02 ), ( 4.7655e-01, 2.96e-02 ), ( 4.9040e-01, 3.03e-02 ), ( 5.0214e-01, 3.08e-02 )
#                                 , ( 4.8480e-01, 3.00e-02 ), ( 4.5579e-01, 2.87e-02 ), ( 4.6245e-01, 2.90e-02 ), ( 4.8362e-01, 3.00e-02 )
#                                 , ( 4.8704e-01, 3.01e-02 ), ( 5.3216e-01, 3.22e-02 ), ( 5.0839e-01, 3.11e-02 ), ( 5.0210e-01, 3.09e-02 )
#                                 , ( 4.3314e-01, 2.76e-02 ), ( 4.8788e-01, 3.02e-02 ), ( 5.0612e-01, 3.11e-02 ), ( 4.6090e-01, 2.90e-02 )
#                                 , ( 4.7729e-01, 2.97e-02 ), ( 4.8237e-01, 3.00e-02 ), ( 4.9961e-01, 3.08e-02 ), ( 5.0405e-01, 3.09e-02 )
#                                ]
#          , ('2011', 'exclB') : [  ( 0.173038,   0.       ), ( 1.4266e-01, 3.12e-02 ), ( 2.1043e-01, 4.17e-02 ), ( 1.7076e-01, 3.56e-02 )
#                                 , ( 2.4038e-01, 4.63e-02 ), ( 2.4343e-01, 4.68e-02 ), ( 2.4987e-01, 4.78e-02 ), ( 2.9090e-01, 5.40e-02 )
#                                 , ( 2.5457e-01, 4.85e-02 ), ( 3.4609e-01, 6.23e-02 ), ( 3.0412e-01, 5.60e-02 ), ( 3.5323e-01, 6.33e-02 )
#                                 , ( 3.0024e-01, 5.54e-02 ), ( 3.6169e-01, 6.46e-02 ), ( 3.6209e-01, 6.47e-02 ), ( 3.6641e-01, 6.53e-02 )
#                                 , ( 4.0317e-01, 7.08e-02 ), ( 3.6273e-01, 6.48e-02 ), ( 3.3338e-01, 6.04e-02 ), ( 4.3366e-01, 7.53e-02 )
#                                 , ( 4.3841e-01, 7.60e-02 ), ( 3.1395e-01, 5.75e-02 ), ( 3.5189e-01, 6.32e-02 ), ( 3.7816e-01, 6.71e-02 )
#                                 , ( 4.0036e-01, 7.04e-02 ), ( 3.6454e-01, 6.51e-02 ), ( 4.4240e-01, 7.66e-02 ), ( 3.9193e-01, 6.92e-02 )
#                                 , ( 4.4756e-01, 7.74e-02 ), ( 4.2747e-01, 7.45e-02 ), ( 3.6052e-01, 6.45e-02 ), ( 4.0084e-01, 7.05e-02 )
#                                 , ( 4.2715e-01, 7.44e-02 ), ( 4.1547e-01, 7.27e-02 ), ( 4.8562e-01, 8.31e-02 ), ( 4.2838e-01, 7.46e-02 )
#                                 , ( 4.3970e-01, 7.63e-02 ), ( 4.0479e-01, 7.11e-02 ), ( 3.8338e-01, 6.80e-02 ), ( 4.1666e-01, 7.29e-02 )
#                                ]
#          , ('2012',  'UB')   : [  ( 0.468739,   0.       ), ( 4.7971e-01, 2.04e-02 ), ( 4.7630e-01, 2.03e-02 ), ( 4.7721e-01, 2.04e-02 )
#                                 , ( 5.0380e-01, 2.12e-02 ), ( 4.8377e-01, 2.06e-02 ), ( 5.0551e-01, 2.13e-02 ), ( 4.7086e-01, 2.02e-02 )
#                                 , ( 4.7154e-01, 2.02e-02 ), ( 4.9224e-01, 2.09e-02 ), ( 4.8656e-01, 2.07e-02 ), ( 4.6919e-01, 2.01e-02 )
#                                 , ( 5.0785e-01, 2.14e-02 ), ( 4.8910e-01, 2.08e-02 ), ( 4.7011e-01, 2.02e-02 ), ( 4.9664e-01, 2.10e-02 )
#                                 , ( 4.8567e-01, 2.07e-02 ), ( 5.0071e-01, 2.11e-02 ), ( 4.7416e-01, 2.03e-02 ), ( 4.7433e-01, 2.03e-02 )
#                                 , ( 4.8371e-01, 2.06e-02 ), ( 4.8369e-01, 2.06e-02 ), ( 4.8215e-01, 2.06e-02 ), ( 4.7833e-01, 2.05e-02 )
#                                 , ( 4.9503e-01, 2.10e-02 ), ( 4.9613e-01, 2.10e-02 ), ( 4.8444e-01, 2.07e-02 ), ( 4.9790e-01, 2.11e-02 )
#                                 , ( 4.8077e-01, 2.06e-02 ), ( 4.9152e-01, 2.09e-02 ), ( 4.6201e-01, 2.00e-02 ), ( 4.9156e-01, 2.09e-02 )
#                                 , ( 4.9884e-01, 2.12e-02 ), ( 4.7740e-01, 2.05e-02 ), ( 4.9598e-01, 2.11e-02 ), ( 4.9286e-01, 2.10e-02 )
#                                 , ( 4.7467e-01, 2.04e-02 ), ( 4.8477e-01, 2.08e-02 ), ( 4.8031e-01, 2.07e-02 ), ( 4.6356e-01, 2.01e-02 )
#                                ]
#          , ('2012', 'exclB') : [  ( 0.200305,   0.       ), ( 2.8623e-01, 4.18e-02 ), ( 3.8356e-01, 5.29e-02 ), ( 3.8379e-01, 5.29e-02 )
#                                 , ( 3.5352e-01, 4.95e-02 ), ( 4.7159e-01, 6.28e-02 ), ( 5.2097e-01, 6.83e-02 ), ( 5.7154e-01, 7.39e-02 )
#                                 , ( 5.4616e-01, 7.11e-02 ), ( 5.8643e-01, 7.56e-02 ), ( 5.8760e-01, 7.57e-02 ), ( 7.3658e-01, 9.20e-02 )
#                                 , ( 6.9035e-01, 8.70e-02 ), ( 6.0130e-01, 7.72e-02 ), ( 7.5205e-01, 9.36e-02 ), ( 6.5083e-01, 8.27e-02 )
#                                 , ( 6.2923e-01, 8.03e-02 ), ( 6.9212e-01, 8.72e-02 ), ( 7.6639e-01, 9.52e-02 ), ( 6.8629e-01, 8.66e-02 )
#                                 , ( 7.0783e-01, 8.89e-02 ), ( 6.5467e-01, 8.31e-02 ), ( 7.3773e-01, 9.22e-02 ), ( 7.8406e-01, 9.71e-02 )
#                                 , ( 7.5833e-01, 9.44e-02 ), ( 7.5448e-01, 9.40e-02 ), ( 7.3442e-01, 9.18e-02 ), ( 7.2189e-01, 9.05e-02 )
#                                 , ( 7.6730e-01, 9.54e-02 ), ( 8.0526e-01, 9.94e-02 ), ( 8.4738e-01, 1.04e-01 ), ( 7.5294e-01, 9.39e-02 )
#                                 , ( 8.8758e-01, 1.07e-01 ), ( 8.8344e-01, 1.07e-01 ), ( 9.3000e-01, 1.10e-01 ), ( 8.0756e-01, 9.97e-02 )
#                                 , ( 7.9179e-01, 9.81e-02 ), ( 8.0093e-01, 9.91e-02 ), ( 9.1059e-01, 1.09e-01 ), ( 9.3188e-01, 1.11e-01 )
#                                ]
#         }

from ROOT import TFile
nTupleFile = TFile.Open( dataPath + nTupleFilePath )
nTuple = nTupleFile.Get(nTupleName)
sumW = { }
sumWSq = { }
for trigger in [ 'UB', 'exclB' ] :
    for period in [ '2011', '2012' ] :
        sumW[ ( period, trigger, ( 0, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumW[ ( period, trigger, ( 1, 0 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumW[ ( period, trigger, ( 1, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 0, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 1, 0 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 1, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        selStr = 'hlt1_excl_biased%s==%d && runPeriod==%s' % ( '' if hlt1TOS else '_dec', trigger == 'exclB', period )
        #selStr = 'hlt1_excl_biased==%d && runPeriod==%s' % ( trigger == 'exclB', period )
        if period == '2011' and period2011 == 'summer' : selStr += '&& runNumber > 87219 && runNumber < 94386'
        elif period == '2011' and period2011 == 'autumn' : selStr += '&& runNumber > 94386'

        dummyFile = TFile.Open( 'dummy.root', 'RECREATE' )
        nTupleSel = nTuple.CopyTree(selStr)
        print 'getting signal yields for %s - %s' % ( period, trigger )
        sumWTot = 0.
        sumWSqTot = 0.
        for ev in nTupleSel :
            timeBin = 0
            for bound in binBounds[ 1 : ] :
                if ev.time < bound : break
                timeBin += 1
            sumW[ ( period, trigger, ( ev.hlt2_unbiased, ev.hlt2_biased ) ) ][timeBin] += ev.sWeights_ipatia
            sumWSq[ ( period, trigger, ( ev.hlt2_unbiased, ev.hlt2_biased ) ) ][timeBin] += ev.sWeights_ipatia**2
            sumWTot += ev.sWeights_ipatia
            sumWSqTot += ev.sWeights_ipatia**2
        print 'found %.0f +/- %.0f signal events for %s - %s (%d entries):'\
              % ( sumWTot, sqrt(sumWSqTot), period, trigger, nTupleSel.GetEntries() )
        for HLT2Cat in [ ( 0, 1 ), ( 1, 0 ), ( 1, 1 ) ] :
            print '    (%d, %d): ' % HLT2Cat,
            for bin in range( len(binBounds) - 1 ) :
                print '%4.0f' % sumW[ ( period, trigger, HLT2Cat ) ][bin],
            print
        dummyFile.Close()

nTupleFile.Close()
import os
os.remove('dummy.root')

countLists = { ( '2011', 'UB' ) : [ ], ( '2012', 'UB' ) : [ ], ( '2011', 'exclB' ) : [ ], ( '2012', 'exclB' ) : [ ] }
countErrLists = { ( '2011', 'UB' ) : [ ], ( '2012', 'UB' ) : [ ], ( '2011', 'exclB' ) : [ ], ( '2012', 'exclB' ) : [ ] }
from array import array
boundsArr = array( 'd', binBounds )
from ROOT import TH1D
for period in [ '2011', '2012' ] :
    for binIt in range( len(binBounds) - 1 ) :
        N_hlt2ExclUnb = sumW[ ( period, 'UB', ( 1, 0 ) ) ][binIt] + sumW[ ( period, 'exclB', ( 1, 0 ) ) ][binIt]
        N_hlt2Both    = sumW[ ( period, 'UB', ( 1, 1 ) ) ][binIt] + sumW[ ( period, 'exclB', ( 1, 1 ) ) ][binIt]
        V_hlt2ExclUnb = sumWSq[ ( period, 'UB', ( 1, 0 ) ) ][binIt] + sumWSq[ ( period, 'exclB', ( 1, 0 ) ) ][binIt]
        V_hlt2Both    = sumWSq[ ( period, 'UB', ( 1, 1 ) ) ][binIt] + sumWSq[ ( period, 'exclB', ( 1, 1 ) ) ][binIt]
        hlt2BiasEff = N_hlt2Both / ( N_hlt2Both + N_hlt2ExclUnb )
        hlt2BiasEffErr = sqrt( N_hlt2ExclUnb**2 * V_hlt2Both + N_hlt2Both**2 * V_hlt2ExclUnb ) / ( N_hlt2Both + N_hlt2ExclUnb )**2
        countLists[ ( period, 'UB' ) ].append(hlt2BiasEff)
        countErrLists[ ( period, 'UB' ) ].append(hlt2BiasEffErr)

        N_hlt1ExclB = sumW[ ( period, 'exclB', ( 1, 1 ) ) ][binIt] + sumW[ ( period, 'exclB', ( 0, 1 ) ) ][binIt]
        if incHLT2UB : N_hlt1ExclB += sumW[ ( period, 'exclB', ( 1, 0 ) ) ][binIt]
        N_hlt1Unb = sumW[ ( period, 'UB', ( 1, 1 ) ) ][binIt] + sumW[ ( period, 'UB', ( 0, 1 ) ) ][binIt]
        if incHLT2UB : N_hlt1Unb += sumW[ ( period, 'UB', ( 1, 0 ) ) ][binIt]
        V_hlt1ExclB = sumWSq[ ( period, 'exclB', ( 1, 1 ) ) ][binIt] + sumWSq[ ( period, 'exclB', ( 0, 1 ) ) ][binIt]
        if incHLT2UB : V_hlt1ExclB += sumWSq[ ( period, 'exclB', ( 1, 0 ) ) ][binIt]
        V_hlt1Unb = sumWSq[ ( period, 'UB', ( 1, 1 ) ) ][binIt] + sumWSq[ ( period, 'UB', ( 0, 1 ) ) ][binIt]
        if incHLT2UB : V_hlt1Unb += sumWSq[ ( period, 'UB', ( 1, 0 ) ) ][binIt]
        hlt1ExclBEff = hlt1UnbEff * N_hlt1ExclB / N_hlt1Unb * hlt2BiasEff
        hlt1ExclBEffErr = hlt1UnbEff * sqrt( N_hlt1Unb**2 * hlt2BiasEff**2 * V_hlt1ExclB + N_hlt1ExclB**2 * hlt2BiasEff**2 * V_hlt1Unb\
                                            + N_hlt1ExclB**2 * N_hlt1Unb**2 * hlt2BiasEffErr**2 ) / N_hlt1Unb**2
        countLists[ ( period, 'exclB' ) ].append(hlt1ExclBEff)
        countErrLists[ ( period, 'exclB' ) ].append(hlt1ExclBEffErr)

print 'plotting acceptance bins'
from P2VV.Load import LHCbStyle
from ROOT import gStyle, TGraphErrors as Graph, TCanvas, kBlack, kBlue, kRed, kFullDotLarge
gStyle.SetEndErrorSize(3)
timeArr = array( 'd', range( 1, len(binBounds) ) )
timeErrArr = array( 'd', [0.] * ( len(binBounds) - 1 ) )
canv = TCanvas('dummy')
canv.Print( plotsFilePath + '[' )
for period in [ '2011', '2012' ] :
    filePathSplit = histsFilePath.split('.')
    histsFile = TFile.Open( '%s_%s.%s' % ( filePathSplit[0], period, filePathSplit[1] ), 'RECREATE' )
    unbName   = histNames[ 'UB_%s'    % period ]
    exclBName = histNames[ 'exclB_%s' % period ]
    hists = dict(  UB    = TH1D( unbName,   unbName,   len(boundsArr) - 1, boundsArr )
                 , exclB = TH1D( exclBName, exclBName, len(boundsArr) - 1, boundsArr )
                )
    for trigger in [ 'UB', 'exclB' ] :
        if accHists :
            accFile = TFile.Open( accHists[period]['file'] )
            hist = accFile.Get( accHists[period][trigger] )
            histList = [ hist.GetBinContent( int(bin) ) for bin in timeArr ]
            scale = float( len(histList) ) / sum(histList) if applyExclBScale and trigger == 'exclB' else 1.
            histArr = array( 'd', [ scale * val for val in histList ] )
            histGraph = Graph( len(histArr), timeArr, histArr )
            histGraph.SetMinimum( plotMinMax[trigger][0] )
            histGraph.SetMaximum( plotMinMax[trigger][1] )
            histGraph.SetMarkerStyle(kFullDotLarge)
            histGraph.SetLineColor(kBlack)
            histGraph.SetMarkerColor(kBlack)
            histGraph.SetLineWidth(3)
            histGraph.SetMarkerSize(1.0)
            histGraph.GetXaxis().SetTitle('time bin')
            histGraph.GetYaxis().SetTitle('acceptance (a.u.)')
        else :
            histGraph = None

        if countLists :
            countList = countLists[ ( period, trigger ) ]
            countErrList = countErrLists[ ( period, trigger ) ]
            scale = float( len(countList) ) / sum(countList) if applyExclBScale and trigger == 'exclB' else 1.
            countArr = array( 'd', [ scale * val for val in countList ] )
            countErrArr = array( 'd', [ scale * val for val in countErrList ] )
            countGraph = Graph( len(timeArr), timeArr, countArr, timeErrArr, countErrArr )
            countGraph.SetMinimum( plotMinMax[trigger][0] )
            countGraph.SetMaximum( plotMinMax[trigger][1] )
            countGraph.SetMarkerStyle(kFullDotLarge)
            countGraph.SetLineColor(kRed)
            countGraph.SetMarkerColor(kRed)
            countGraph.SetLineWidth(3)
            countGraph.SetMarkerSize(1.0)
            countGraph.GetXaxis().SetTitle('time bin')
            countGraph.GetYaxis().SetTitle('acceptance (a.u.)')
            if histsType != 'fit' :
                for binIt, val in enumerate(histArr) : hists[trigger].SetBinContent( binIt + 1, val )
        else :
            countGraph = None

        if fitAcc :
            fitList = fitAcc[ ( period, trigger ) ]
            scale = float( len(fitList) ) / sum( vals[0] for vals in fitList ) if applyExclBScale and trigger == 'exclB' else 1.
            fitArr = array( 'd', [ scale * vals[0] for vals in fitList ] )
            fitErrArr = array( 'd', [ scale * vals[1] for vals in fitList ] )
            fitGraph = Graph( len(timeArr), timeArr, fitArr, timeErrArr, fitErrArr )
            fitGraph.SetMinimum( plotMinMax[trigger][0] )
            fitGraph.SetMaximum( plotMinMax[trigger][1] )
            fitGraph.SetMarkerStyle(kFullDotLarge)
            fitGraph.SetLineColor(kBlue)
            fitGraph.SetMarkerColor(kBlue)
            fitGraph.SetLineWidth(3)
            fitGraph.SetMarkerSize(1.0)
            fitGraph.GetXaxis().SetTitle('time bin')
            fitGraph.GetYaxis().SetTitle('acceptance (a.u.)')
            if histsType == 'fit' :
                for binIt, val in enumerate(fitArr) : hists[trigger].SetBinContent( binIt + 1, val )
        else :
            fitGraph = None

        canv = TCanvas( '%s_%s' % ( period, trigger ) )
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)
        drawOpts = 'ALP'
        if histGraph :
            histGraph.Draw(drawOpts)
            drawOpts = 'LP SAMES'
        if countGraph :
            countGraph.Draw(drawOpts)
            drawOpts = 'LP SAMES'
        if fitGraph :
            fitGraph.Draw(drawOpts)
        canv.Print(plotsFilePath)

        canv = TCanvas( '%s_%s_hist' % ( period, trigger ) )
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)
        hists[trigger].Draw()
        canv.Print(plotsFilePath)

    histsFile.Write()
    histsFile.Close()

canv.Print( plotsFilePath + ']' )
