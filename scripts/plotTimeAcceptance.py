from math import sqrt
plotsFilePath = 'timeAcc.pdf'
histsFilePath = 'timeAcc.root'
histsType = 'fit' # 'hists' # 'fit'
dataPath = '/project/bfys/jleerdam/data/Bs2Jpsiphi/Reco14/'
nTupleFilePath = 'fitNTuple_peakBkg_2011_2012_Reco14_TOS_20140822_with_IPz_weights_v2_sigWeights.root' #'fitNTuple_peakBkg_2011_2012_Reco14_TOS_20140309.root'
#nTupleFilePath = 'DEC/fitNTuple_peakBkg_2011_2012_Reco14_20140116.root'
nTupleName = 'DecayTree'
period2011 = '' # 'summer'
weightName = 'sigWeight' #'sWeights_ipatia'
applyExclBScale = False
plotBins = False
hlt1TOS = True
incHLT2UB = False
hlt1UnbEff = 0.7
plotMinMax = dict( UB = ( 0.9, 1.0 ), exclB = ( 0.25, 1.55 ) if applyExclBScale else ( 0.05, 0.22 ) )
timeMinMax = ( 0., 7.1 )
divisions = dict( UB = ( 5, 5, 0 ), exclB = ( 8, 5, 0 ) )
LHCbLabel = '' # 'LHCb'
yTitle = dict(  UB    = 'Acceptance'
              , exclB = 'Acceptance (arb. scale)'
             )

drawHists = False
drawCounts = True
drawFit = True

#binBounds = [ 0.3, 14.0 ]
#binBounds = [ 0.3, 0.5, 1.0, 1.7, 2.7, 14.0 ]
#binBounds = [ 0.3, 0.46, 0.63, 0.83, 1.05, 1.32, 1.65, 2.07, 2.67, 3.7, 14. ]
#binBounds = [ 0.3, 0.403, 0.512, 0.636, 0.759, 0.896, 1.047, 1.211, 1.403, 1.622, 1.869, 2.170, 2.540, 3.033, 3.759, 5.321, 14. ]
binBounds = [  0.3,            0.337264386972, 0.37549669625,  0.414748554299, 0.455075831456, 0.496539120209, 0.53920428282
             , 0.583143080703, 0.628433900662, 0.675162596446, 0.72342346832,  0.773320408723, 0.824968249017, 0.878494351211
             , 0.934040500166, 0.991765167029, 1.05184623487,  1.11448430465,  1.17990673635,  1.24837263066,  1.32017902684
             , 1.39566869131,  1.47524001355,  1.55935973292,  1.64857952745,  1.74355796227,  1.84509002002,  1.95414759126
             , 2.07193619997,  2.19997646239,  2.34022446712,  2.49525577287,  2.6685581881,   2.86502098649,  3.09180358141
             , 3.36000284338,  3.68820308685,  4.11122609838,  4.70718114264,  5.724828229,    14.
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
# 40 bins, uniform UB (2011: 0.4868 +/- 0.0032; 2012: 0.4429 +/- 0.0021)
#fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.06553316 , 0.0080358395 ), ( 0.057065788, 0.007813602  ), ( 0.086664293, 0.009312648  )
#                                             , ( 0.063884553, 0.0075026504 ), ( 0.083185569, 0.0081742005 ), ( 0.085844445, 0.0081734907 )
#                                             , ( 0.096111367, 0.0087083272 ), ( 0.10697933 , 0.0089212501 ), ( 0.095740913, 0.008508349  )
#                                             , ( 0.13040634 , 0.0099084137 ), ( 0.10979424 , 0.0089189438 ), ( 0.13359363 , 0.010031347  )
#                                             , ( 0.11568291 , 0.0092381836 ), ( 0.1334907  , 0.0098341591 ), ( 0.13915294 , 0.010132405  )
#                                             , ( 0.13465431 , 0.0098295797 ), ( 0.15224258 , 0.01070236   ), ( 0.13470815 , 0.010056801  )
#                                             , ( 0.12770127 , 0.0097965399 ), ( 0.16070856 , 0.010926957  ), ( 0.15357037 , 0.010521502  )
#                                             , ( 0.11749946 , 0.0092924914 ), ( 0.1353639  , 0.01001682   ), ( 0.13755885 , 0.0099681597 )
#                                             , ( 0.14627491 , 0.010316067  ), ( 0.14223214 , 0.010552815  ), ( 0.16489676 , 0.011284534  )
#                                             , ( 0.14814479 , 0.010495572  ), ( 0.16253116 , 0.011124742  ), ( 0.15194141 , 0.010379752  )
#                                             , ( 0.1329831  , 0.0097480511 ), ( 0.14426854 , 0.010251041  ), ( 0.16384392 , 0.011373492  )
#                                             , ( 0.15408801 , 0.010708498  ), ( 0.17837116 , 0.011496694  ), ( 0.15903125 , 0.010975913  )
#                                             , ( 0.16355449 , 0.01114966   ), ( 0.14933061 , 0.010570292  ), ( 0.13072383 , 0.0098250496 )
#                                             , ( 0.153867   , 0.01066002   )
#                                            ]
#                , ('2012', 'hlt1_exclB' ) : [  ( 0.056782311, 0.0054992752 ), ( 0.072817631, 0.0059746117 ), ( 0.097233896, 0.0066510545 )
#                                             , ( 0.094286388, 0.0063634097 ), ( 0.086469098, 0.0060181831 ), ( 0.12677357 , 0.0071039412 )
#                                             , ( 0.13217053 , 0.0071898263 ), ( 0.14484951 , 0.0074707837 ), ( 0.13858122 , 0.0072334277 )
#                                             , ( 0.1453542  , 0.0073138873 ), ( 0.15037712 , 0.0073742835 ), ( 0.18477648 , 0.0082829213 )
#                                             , ( 0.17008072 , 0.0076894125 ), ( 0.15096312 , 0.0073328988 ), ( 0.19431242 , 0.0084115733 )
#                                             , ( 0.16323424 , 0.0075576594 ), ( 0.15949884 , 0.0075660974 ), ( 0.16974072 , 0.0076927342 )
#                                             , ( 0.19088565 , 0.0083247745 ), ( 0.17326526 , 0.0078776208 ), ( 0.17488658 , 0.0078759148 )
#                                             , ( 0.15775141 , 0.0074577549 ), ( 0.18393531 , 0.0081117827 ), ( 0.19457925 , 0.0083879636 )
#                                             , ( 0.18121628 , 0.0079797274 ), ( 0.18738883 , 0.0081090502 ), ( 0.18018204 , 0.0079538955 )
#                                             , ( 0.17754842 , 0.007932405  ), ( 0.18802795 , 0.0082168762 ), ( 0.19817746 , 0.008421722  )
#                                             , ( 0.2086808  , 0.008881189  ), ( 0.18630514 , 0.0082787984 ), ( 0.21260243 , 0.0087406632 )
#                                             , ( 0.21648677 , 0.0089808773 ), ( 0.22524265 , 0.0090796571 ), ( 0.19486311 , 0.0083775725 )
#                                             , ( 0.19634549 , 0.0085824185 ), ( 0.1972288  , 0.0085389328 ), ( 0.22093802 , 0.0091576646 )
#                                             , ( 0.2286749  , 0.0095578491 )
#                                            ]
#                , ('2011', 'hlt2_B' )     : [  ( 0.95646608, 0.021676572  ), ( 0.90043864, 0.023807196  ), ( 0.94663658, 0.032131384  )
#                                             , ( 0.97396273, 0.011952464  ), ( 0.9923056 , 0.0065942815 ), ( 0.9941027 , 0.0066649872 )
#                                             , ( 0.97397375, 0.010450416  ), ( 0.9975974 , 0.0076099974 ), ( 0.980143  , 0.0089766316 )
#                                             , ( 0.99476449, 0.0056561543 ), ( 0.98771929, 0.0073376447 ), ( 0.97768068, 0.0087038514 )
#                                             , ( 0.98561814, 0.0086703258 ), ( 0.97639175, 0.008986086  ), ( 0.96602462, 0.010332794  )
#                                             , ( 0.98787222, 0.006929685  ), ( 0.97041726, 0.0099562975 ), ( 0.96102078, 0.010493887  )
#                                             , ( 0.96696945, 0.0093416622 ), ( 0.97194125, 0.0089705556 ), ( 0.97315713, 0.008588489  )
#                                             , ( 0.95838225, 0.011442574  ), ( 0.95440078, 0.010958396  ), ( 0.95009583, 0.011172389  )
#                                             , ( 0.96094089, 0.010422217  ), ( 0.94115894, 0.012918544  ), ( 0.9408771 , 0.012754815  )
#                                             , ( 0.96399639, 0.010708339  ), ( 0.95725505, 0.01071411   ), ( 0.97615851, 0.0080561805 )
#                                             , ( 0.96916612, 0.0094056759 ), ( 0.95134192, 0.010779428  ), ( 0.95401083, 0.012784899  )
#                                             , ( 0.96821215, 0.0096667846 ), ( 0.96430947, 0.0095886325 ), ( 0.96927458, 0.010077415  )
#                                             , ( 0.96236765, 0.010541662  ), ( 0.95421669, 0.010899896  ), ( 0.980848  , 0.0077663318 )
#                                             , ( 0.98574418, 0.0079225785 )
#                                            ]
#                , ('2012', 'hlt2_B' )     : [  ( 0.94766515, 0.01506495   ), ( 0.96891688, 0.011531825  ), ( 0.97040945, 0.010622542  )
#                                             , ( 0.99020824, 0.0097190638 ), ( 0.98272311, 0.0066775051 ), ( 0.98309566, 0.0074680164 )
#                                             , ( 0.98137202, 0.0062818682 ), ( 0.97962784, 0.0072116676 ), ( 0.9906472 , 0.0071437605 )
#                                             , ( 0.97908444, 0.0064008249 ), ( 0.98729999, 0.0052721955 ), ( 0.9836951 , 0.0061692743 )
#                                             , ( 0.99196351, 0.0041901402 ), ( 0.98836702, 0.0051531016 ), ( 0.99321693, 0.0056001658 )
#                                             , ( 0.98766   , 0.0049837337 ), ( 0.98341904, 0.0056874619 ), ( 0.98078783, 0.0056062456 )
#                                             , ( 0.97812867, 0.0061474632 ), ( 0.97028256, 0.0068429306 ), ( 0.96864652, 0.0067250146 )
#                                             , ( 0.98835462, 0.004810783  ), ( 0.97849674, 0.0059224858 ), ( 0.98219987, 0.006105578  )
#                                             , ( 0.97971763, 0.0055234294 ), ( 0.98388981, 0.0050112735 ), ( 0.98758706, 0.0048331035 )
#                                             , ( 0.98353624, 0.0049725123 ), ( 0.9703552 , 0.0064593813 ), ( 0.98657511, 0.0051614835 )
#                                             , ( 0.97902852, 0.0063803135 ), ( 0.98010596, 0.0053496259 ), ( 0.97867627, 0.0054122423 )
#                                             , ( 0.98046719, 0.0056511278 ), ( 0.97064747, 0.0061893351 ), ( 0.98244466, 0.0053107924 )
#                                             , ( 0.97354332, 0.0067700998 ), ( 0.97610142, 0.00630721   ), ( 0.97513385, 0.0061309931 )
#                                             , ( 0.95969296, 0.0084767992 )
#                                            ]
#               }

# efficiency weights: 40 bins, uniform UB (2011: 0.4868 +/- 0.0032; 2012: 0.4429 +/- 0.0021)
fitAccLevels = {  ('2011', 'hlt1_exclB' ) : [  ( 0.06568702 , 0.0080700124 ), ( 0.05718789 , 0.0078451148 ), ( 0.086867609, 0.0093405674 )
                                             , ( 0.064021853, 0.0075319107 ), ( 0.083346913, 0.0082030808 ), ( 0.085975112, 0.0082003189 )
                                             , ( 0.096328097, 0.0087393356 ), ( 0.10711781 , 0.0089464313 ), ( 0.095897686, 0.0085328216 )
                                             , ( 0.13060023 , 0.0099341106 ), ( 0.10992871 , 0.0089386106 ), ( 0.13372579 , 0.010051205  )
                                             , ( 0.11577142 , 0.0092560339 ), ( 0.13353034 , 0.0098498801 ), ( 0.13920072 , 0.010148601  )
                                             , ( 0.13464458 , 0.0098413616 ), ( 0.15230611 , 0.010718706  ), ( 0.13470205 , 0.010066461  )
                                             , ( 0.12755927 , 0.0098010025 ), ( 0.16061039 , 0.010931461  ), ( 0.153431   , 0.010522391  )
                                             , ( 0.11736594 , 0.0092924022 ), ( 0.13520346 , 0.0100132    ), ( 0.13732295 , 0.0099604959 )
                                             , ( 0.1459978  , 0.010306094  ), ( 0.14193551 , 0.010541006  ), ( 0.16447774 , 0.011265075  )
                                             , ( 0.14772953 , 0.010473768  ), ( 0.16194456 , 0.011094377  ), ( 0.15150349 , 0.010352746  )
                                             , ( 0.13248423 , 0.0097151208 ), ( 0.14366409 , 0.010210596  ), ( 0.16306357 , 0.011318446  )
                                             , ( 0.15330573 , 0.010653215  ), ( 0.17719356 , 0.011415415  ), ( 0.15807181 , 0.010897596  )
                                             , ( 0.16232924 , 0.011050119  ), ( 0.14788981 , 0.010444581  ), ( 0.12928917 , 0.0096924523 )
                                             , ( 0.15748438 , 0.010811559  )
                                            ]
                , ('2012', 'hlt1_exclB' ) : [  ( 0.057014902, 0.0055381996 ), ( 0.073065463, 0.006013531  ), ( 0.097462053, 0.0066863859 )
                                             , ( 0.094529209, 0.0063958207 ), ( 0.086705848, 0.006048904  ), ( 0.12711356 , 0.0071364514 )
                                             , ( 0.13244126 , 0.0072168099 ), ( 0.14508011 , 0.0074989913 ), ( 0.13879574 , 0.007257255  )
                                             , ( 0.14560769 , 0.0073397809 ), ( 0.15061087 , 0.0073973352 ), ( 0.18494884 , 0.0083040836 )
                                             , ( 0.17021206 , 0.0077091493 ), ( 0.15114305 , 0.0073510386 ), ( 0.19449342 , 0.0084313873 )
                                             , ( 0.16331395 , 0.0075716059 ), ( 0.15960842 , 0.0075791539 ), ( 0.16978705 , 0.0077038652 )
                                             , ( 0.19080731 , 0.0083310282 ), ( 0.1730768  , 0.0078785556 ), ( 0.17481503 , 0.0078803299 )
                                             , ( 0.15757396 , 0.0074552816 ), ( 0.18374704 , 0.0081090055 ), ( 0.19434428 , 0.0083818311 )
                                             , ( 0.18098143 , 0.0079719015 ), ( 0.18706388 , 0.0080959185 ), ( 0.17968106 , 0.0079334868 )
                                             , ( 0.17709658 , 0.0079117197 ), ( 0.18747781 , 0.0081902828 ), ( 0.19758871 , 0.0083914602 )
                                             , ( 0.20778166 , 0.0088402167 ), ( 0.18538003 , 0.0082368473 ), ( 0.21147866 , 0.0086830961 )
                                             , ( 0.21513635 , 0.0089078253 ), ( 0.22396405 , 0.0090062568 ), ( 0.19341499 , 0.0082935519 )
                                             , ( 0.194378   , 0.0084747912 ), ( 0.19494576 , 0.0084094956 ), ( 0.21786814 , 0.0089742621 )
                                             , ( 0.22623403 , 0.0093296186 )
                                            ]
                , ('2011', 'hlt2_B' )     : [  
                                               ( 0.9574491 , 0.021364623  ), ( 0.90173573, 0.023657698  ), ( 0.94838472, 0.031801665  )
                                             , ( 0.97420944, 0.011853904  ), ( 0.99237827, 0.006538489  ), ( 0.99415915, 0.0066002425 )
                                             , ( 0.97410425, 0.010402811  ), ( 0.99763527, 0.0074961535 ), ( 0.98022397, 0.0089421237 )
                                             , ( 0.9947953 , 0.0056201286 ), ( 0.98776544, 0.0073073102 ), ( 0.97775932, 0.0086751143 )
                                             , ( 0.98565144, 0.0086475987 ), ( 0.97642907, 0.008963519  ), ( 0.9660683 , 0.010319535  )
                                             , ( 0.98787211, 0.006925313  ), ( 0.97044008, 0.0099509404 ), ( 0.96098795, 0.010498968  )
                                             , ( 0.96693011, 0.00935098   ), ( 0.97192624, 0.0089740285 ), ( 0.97312416, 0.0085956568 )
                                             , ( 0.95829584, 0.011461653  ), ( 0.95443068, 0.010951518  ), ( 0.95002339, 0.011183672  )
                                             , ( 0.96087559, 0.010434869  ), ( 0.94105835, 0.012936446  ), ( 0.94075938, 0.012773453  )
                                             , ( 0.96387574, 0.010737778  ), ( 0.95713323, 0.010739406  ), ( 0.97611534, 0.0080673242 )
                                             , ( 0.96911436, 0.0094131873 ), ( 0.95129886, 0.010783847  ), ( 0.95394263, 0.012809049  )
                                             , ( 0.96816245, 0.0096748103 ), ( 0.96432371, 0.0095741884 ), ( 0.96937416, 0.01004092   )
                                             , ( 0.96268374, 0.010460093  ), ( 0.95481479, 0.010755828  ), ( 0.98137772, 0.0075402622 )
                                             , ( 0.98752202, 0.006953996  )
                                            ]
                , ('2012', 'hlt2_B' )     : [  
                                               ( 0.94873951, 0.014877947  ), ( 0.96954962, 0.011339398  ), ( 0.97095117, 0.010458229  )
                                             , ( 0.99064251, 0.0093270589 ), ( 0.98289243, 0.0066136322 ), ( 0.9833428 , 0.0073635204 )
                                             , ( 0.98153266, 0.0062315187 ), ( 0.97987276, 0.0071331998 ), ( 0.99084545, 0.006988229  )
                                             , ( 0.97925456, 0.0063518127 ), ( 0.98738853, 0.005231654  ), ( 0.98387216, 0.0061066302 )
                                             , ( 0.99199519, 0.0041653176 ), ( 0.98845562, 0.0051138739 ), ( 0.99331235, 0.0055179126 )
                                             , ( 0.98772483, 0.0049551058 ), ( 0.98350479, 0.0056565265 ), ( 0.98085523, 0.0055837003 )
                                             , ( 0.9782768 , 0.0061081296 ), ( 0.97040362, 0.0068143826 ), ( 0.96874898, 0.0067030849 )
                                             , ( 0.98840995, 0.0047853757 ), ( 0.97859497, 0.0058969721 ), ( 0.98223646, 0.0060849064 )
                                             , ( 0.97975967, 0.0055098685 ), ( 0.98394265, 0.0049893566 ), ( 0.98764461, 0.004807003  )
                                             , ( 0.98355256, 0.004963001  ), ( 0.97041423, 0.0064440513 ), ( 0.98659401, 0.0051438241 )
                                             , ( 0.97912876, 0.0063472909 ), ( 0.98017423, 0.0053269682 ), ( 0.9787403 , 0.0053902274 )
                                             , ( 0.98050896, 0.0056277138 ), ( 0.97079313, 0.0061529402 ), ( 0.98249693, 0.0052807403 )
                                             , ( 0.97363055, 0.0067202703 ), ( 0.97641567, 0.0062154921 ), ( 0.97577487, 0.0059814767 )
                                             , ( 0.96261038, 0.0078636598 )
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
sumWBins = [ 0. ] * ( len(binBounds) - 1 )
timeVals = [ 0. ] * ( len(binBounds) - 1 )
for trigger in [ 'UB', 'exclB' ] :
    for period in [ '2011', '2012' ] :
        sumW[ ( period, trigger, ( 0, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumW[ ( period, trigger, ( 1, 0 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumW[ ( period, trigger, ( 1, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 0, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 1, 0 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        sumWSq[ ( period, trigger, ( 1, 1 ) ) ] = [ 0. ] * ( len(binBounds) - 1 )
        selStr = 'hlt1_excl_biased%s==%d && runPeriod==%s' % ( '' if hlt1TOS else '_dec', trigger == 'exclB', period )
        if period == '2011' and period2011 == 'summer' : selStr += '&& firstData == 1'
        elif period == '2011' and period2011 == 'autumn' : selStr += '&& firstData == 0'

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
            sigWeight = getattr( ev, weightName ) if weightName else 1.
            sumW[ ( period, trigger, ( ev.hlt2_unbiased, ev.hlt2_biased ) ) ][timeBin] += sigWeight
            sumWSq[ ( period, trigger, ( ev.hlt2_unbiased, ev.hlt2_biased ) ) ][timeBin] += sigWeight**2
            sumWTot += sigWeight
            sumWSqTot += sigWeight**2
            sumWBins[timeBin] += sigWeight
            timeVals[timeBin] += sigWeight * ev.time
        print 'found %.0f +/- %.0f signal events for %s - %s (%d entries):'\
              % ( sumWTot, sqrt(sumWSqTot), period, trigger, nTupleSel.GetEntries() )
        for HLT2Cat in [ ( 0, 1 ), ( 1, 0 ), ( 1, 1 ) ] :
            print '    (%d, %d): ' % HLT2Cat,
            for bin in range( len(binBounds) - 1 ) :
                print '%4.0f' % sumW[ ( period, trigger, HLT2Cat ) ][bin],
            print
        dummyFile.Close()
timeVals = [ val / sumWBin for val, sumWBin in zip( timeVals, sumWBins ) ]

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
        hlt2BiasEff = ( N_hlt2Both / ( N_hlt2Both + N_hlt2ExclUnb ) ) if N_hlt2Both + N_hlt2ExclUnb > 0. else 0.
        hlt2BiasEffErr = ( sqrt( N_hlt2ExclUnb**2 * V_hlt2Both + N_hlt2Both**2 * V_hlt2ExclUnb ) / ( N_hlt2Both + N_hlt2ExclUnb )**2 )\
                         if N_hlt2Both + N_hlt2ExclUnb > 0. else 0.
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
        hlt1ExclBEff = ( hlt1UnbEff * N_hlt1ExclB / N_hlt1Unb * hlt2BiasEff ) if N_hlt1Unb > 0. else 0.
        hlt1ExclBEffErr = ( hlt1UnbEff * sqrt( N_hlt1Unb**2 * hlt2BiasEff**2 * V_hlt1ExclB + N_hlt1ExclB**2 * hlt2BiasEff**2 * V_hlt1Unb\
                                              + N_hlt1ExclB**2 * N_hlt1Unb**2 * hlt2BiasEffErr**2 ) / N_hlt1Unb**2 )\
                          if N_hlt1Unb > 0. else 0.
        countLists[ ( period, 'exclB' ) ].append(hlt1ExclBEff)
        countErrLists[ ( period, 'exclB' ) ].append(hlt1ExclBEffErr)

print 'plotting acceptance bins'
from P2VV.Load import LHCbStyle
from ROOT import gStyle, TGraphAsymmErrors as Graph, TCanvas, kBlack, kBlue, kRed, kFullDotLarge
gStyle.SetColorModelPS(1)
gStyle.SetEndErrorSize( 3 if plotBins else 2 )

from ROOT import TLatex
label = TLatex()
label.SetTextAlign(32)
label.SetTextSize(0.072)

if plotBins :
    timeArr = array( 'd', range( 1, len(binBounds) ) )
    timeErrLoArr = array( 'd', [ 0. ] * ( len(binBounds) - 1 ) )
    timeErrHiArr = timeErrLoArr
else :
    timeArr = array( 'd', timeVals )
    timeErrLoArr = array( 'd', [ val - binBounds[it] for it, val in enumerate(timeVals) ] )
    timeErrHiArr = array( 'd', [ binBounds[ it + 1 ] - val for it, val in enumerate(timeVals) ] )

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
        if drawHists :
            accFile = TFile.Open( accHists[period]['file'] )
            hist = accFile.Get( accHists[period][trigger] )
            histList = [ hist.GetBinContent( int(bin) ) for bin in range( 1, len(binBounds) ) ]
            scale = float( len(histList) ) / sum(histList) if applyExclBScale and trigger == 'exclB' else 1.
            histArr = array( 'd', [ scale * val for val in histList ] )
            histErrArr = array( 'd', [ 0. ] * len(histArr) )
            histGraph = Graph( len(timeArr), timeArr, histArr, timeErrLoArr, timeErrHiArr, histErrArr, histErrArr )
            if not plotBins :
                histGraph.GetXaxis().SetLimits( timeMinMax[0], timeMinMax[1] )
            histGraph.SetMinimum( plotMinMax[trigger][0] )
            histGraph.SetMaximum( plotMinMax[trigger][1] )
            histGraph.SetMarkerStyle(kFullDotLarge)
            histGraph.SetLineColor(kBlack)
            histGraph.SetMarkerColor(kBlack)
            histGraph.SetLineWidth( 3 if plotBins else 2 )
            histGraph.SetMarkerSize( 1.0 if plotBins else 0.6 )
            histGraph.GetXaxis().SetTitleOffset(1.1)
            histGraph.GetYaxis().SetTitleOffset(1.0)
            histGraph.GetXaxis().SetTitle( 'Decay-time bin' if plotBins else 'Decay time [ps]' )
            histGraph.GetYaxis().SetTitle( yTitle[trigger] )
            histGraph.GetYaxis().SetNdivisions( divisions[trigger][0], divisions[trigger][1], divisions[trigger][2] )
        else :
            histGraph = None

        if drawCounts :
            countList = countLists[ ( period, trigger ) ]
            countErrList = countErrLists[ ( period, trigger ) ]
            scale = ( float( len(countList) ) / sum(countList) ) if applyExclBScale and trigger == 'exclB' and sum(countList) > 0. else 1.
            countArr = array( 'd', [ scale * val for val in countList ] )
            countErrArr = array( 'd', [ scale * val for val in countErrList ] )
            countGraph = Graph( len(timeArr), timeArr, countArr, timeErrLoArr, timeErrHiArr, countErrArr, countErrArr )
            if not plotBins :
                countGraph.GetXaxis().SetLimits( timeMinMax[0], timeMinMax[1] )
            countGraph.SetMinimum( plotMinMax[trigger][0] )
            countGraph.SetMaximum( plotMinMax[trigger][1] )
            countGraph.SetMarkerStyle(kFullDotLarge)
            color = kRed if drawFit else kBlue
            countGraph.SetLineColor(color)
            countGraph.SetMarkerColor(color)
            countGraph.SetLineWidth( 3 if plotBins else 2 )
            countGraph.SetMarkerSize( 1.0 if plotBins else 0.6 )
            countGraph.GetXaxis().SetTitleOffset(1.1)
            countGraph.GetYaxis().SetTitleOffset(1.0)
            countGraph.GetXaxis().SetTitle( 'Decay-time bin' if plotBins else 'Decay time [ps]' )
            countGraph.GetYaxis().SetTitle( yTitle[trigger] )
            countGraph.GetYaxis().SetNdivisions( divisions[trigger][0], divisions[trigger][1], divisions[trigger][2] )
            if histsType != 'fit' :
                for binIt, val in enumerate(countArr) : hists[trigger].SetBinContent( binIt + 1, val )
        else :
            countGraph = None

        if drawFit :
            fitList = fitAcc[ ( period, trigger ) ]
            scale = float( len(fitList) ) / sum( vals[0] for vals in fitList ) if applyExclBScale and trigger == 'exclB' else 1.
            #scale = 0.99 / max( vals[0] for vals in fitList ) if applyExclBScale and trigger == 'exclB' else 1.
            fitArr = array( 'd', [ scale * vals[0] for vals in fitList ] )
            fitErrArr = array( 'd', [ scale * vals[1] for vals in fitList ] )
            fitGraph = Graph( len(timeArr), timeArr, fitArr, timeErrLoArr, timeErrHiArr, fitErrArr, fitErrArr )
            if not plotBins :
                fitGraph.GetXaxis().SetLimits( timeMinMax[0], timeMinMax[1] )
            fitGraph.SetMinimum( plotMinMax[trigger][0] )
            fitGraph.SetMaximum( plotMinMax[trigger][1] )
            fitGraph.SetMarkerStyle(kFullDotLarge)
            fitGraph.SetLineColor(kBlue)
            fitGraph.SetMarkerColor(kBlue)
            fitGraph.SetLineWidth( 3 if plotBins else 2 )
            fitGraph.SetMarkerSize( 1.0 if plotBins else 0.6 )
            fitGraph.GetXaxis().SetTitleOffset(1.1)
            fitGraph.GetYaxis().SetTitleOffset(1.0)
            fitGraph.GetXaxis().SetTitle( 'Decay-time bin' if plotBins else 'Decay time [ps]' )
            fitGraph.GetYaxis().SetTitle( yTitle[trigger] )
            fitGraph.GetYaxis().SetNdivisions( divisions[trigger][0], divisions[trigger][1], divisions[trigger][2] )
            if histsType == 'fit' :
                for binIt, val in enumerate(fitArr) : hists[trigger].SetBinContent( binIt + 1, val )
        else :
            fitGraph = None

        canv = TCanvas( '%s_%s' % ( period, trigger ) )
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)
        drawOpts = 'ALP' if plotBins else 'APE1'
        if histGraph :
            histGraph.Draw(drawOpts)
            drawOpts = 'LP SAMES' if plotBins else 'PE1 SAMES'
        if countGraph :
            countGraph.Draw(drawOpts)
            drawOpts = 'LP SAMES' if plotBins else 'PE1 SAMES'
        if fitGraph :
            fitGraph.Draw(drawOpts)
        label.DrawLatexNDC( 0.88, 0.3, LHCbLabel )
        canv.Print(plotsFilePath)

        canv = TCanvas( '%s_%s_hist' % ( period, trigger ) )
        canv.SetLeftMargin(0.18)
        canv.SetRightMargin(0.05)
        canv.SetBottomMargin(0.18)
        canv.SetTopMargin(0.05)
        hists[trigger].Draw()
        label.DrawLatexNDC( 0.88, 0.3, LHCbLabel )
        canv.Print(plotsFilePath)

    histsFile.Write()
    histsFile.Close()

canv.Print( plotsFilePath + ']' )
