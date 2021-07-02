\ENCODING=ISO-8859-1
\Problem name: 

Maximize
 c7c4cc4dmda7em11ebm8574m50eb71d9aedf#0: BIOMASS_Ecoli_core_w_GAM#48
                                 - BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + 0 PFK_potential#196 + 0 PFL_potential#203
                                 + 0 PGI_potential#209 + 0 PGK_potential#212
                                 + 0 PGL_potential#216 + 0 ACALD_potential#221
                                 + 0 AKGt2r_potential#226 + 0 PGM_potential#231
                                 + 0 PIt2r_potential#234
                                 + 0 ALCD2x_potential#238
                                 + 0 ACALDt_potential#241
                                 + 0 ACKr_potential#244 + 0 PPC_potential#248
                                 + 0 ACONTa_potential#253
                                 + 0 ACONTb_potential#257
                                 + 0 ATPM_potential#260 + 0 PPCK_potential#262
                                 + 0 ACt2r_potential#264 + 0 PPS_potential#267
                                 + 0 ADK1_potential#270 + 0 AKGDH_potential#272
                                 + 0 ATPS4r_potential#275
                                 + 0 PTAr_potential#277 + 0 PYK_potential#279
                                 + 0 BIOMASS_Ecoli_core_w_GAM_potential#281
                                 + 0 PYRt2_potential#290 + 0 CO2t_potential#293
                                 + 0 RPE_potential#296 + 0 CS_potential#300
                                 + 0 RPI_potential#302
                                 + 0 SUCCt2_2_potential#304
                                 + 0 CYTBD_potential#308
                                 + 0 D_LACt2_potential#313
                                 + 0 _ENO_potential#317
                                 + 0 SUCCt3_potential#319
                                 + 0 _ETOHt2r_potential#321
                                 + 0 SUCDi_potential#324
                                 + 0 SUCOAS_potential#327
                                 + 0 TALA_potential#329 + 0 THD2_potential#332
                                 + 0 TKT1_potential#334 + 0 TKT2_potential#336
                                 + 0 TPI_potential#338
                                 + 0 _EX_ac_e_potential#341
                                 + 0 _EX_acald_e_potential#343
                                 + 0 _EX_akg_e_potential#345
                                 + 0 _EX_co2_e_potential#347
                                 + 0 _EX_etoh_e_potential#349
                                 + 0 _EX_for_e_potential#351
                                 + 0 _EX_fru_e_potential#354
                                 + 0 _EX_fum_e_potential#357
                                 + 0 _EX_glc__D_e_potential#360
                                 + 0 _EX_gln__L_e_potential#363
                                 + 0 _EX_glu__L_e_potential#366
                                 + 0 _EX_h_e_potential#369
                                 + 0 _EX_h2o_e_potential#371
                                 + 0 _EX_lac__D_e_potential#374
                                 + 0 _EX_mal__L_e_potential#376
                                 + 0 _EX_nh4_e_potential#379
                                 + 0 _EX_o2_e_potential#382
                                 + 0 _EX_pi_e_potential#385
                                 + 0 _EX_pyr_e_potential#387
                                 + 0 _EX_succ_e_potential#389
                                 + 0 FBA_potential#391 + 0 FBP_potential#393
                                 + 0 FORt2_potential#395 + 0 FORt_potential#397
                                 + 0 FRD7_potential#399
                                 + 0 FRUpts2_potential#401
                                 + 0 FUM_potential#403
                                 + 0 FUMt2_2_potential#406
                                 + 0 G6PDH2r_potential#408
                                 + 0 GAPD_potential#410
                                 + 0 GLCpts_potential#412
                                 + 0 GLNS_potential#414
                                 + 0 GLNabc_potential#417
                                 + 0 GLUDy_potential#419 + 0 GLUN_potential#421
                                 + 0 GLUSy_potential#423
                                 + 0 GLUt2r_potential#425 + 0 GND_potential#427
                                 + 0 H2Ot_potential#429
                                 + 0 ICDHyr_potential#431 + 0 ICL_potential#433
                                 + 0 LDH_D_potential#436 + 0 MALS_potential#438
                                 + 0 MALt2_2_potential#440
                                 + 0 MDH_potential#442 + 0 ME1_potential#444
                                 + 0 ME2_potential#446 + 0 NADH16_potential#448
                                 + 0 NADTRHD_potential#450
                                 + 0 NH4t_potential#452 + 0 O2t_potential#454
                                 + 0 PDH_potential#456
Subject To
 glc__D_e#0:                     - _EX_glc__D_e#102
                                 + _EX_glc__D_e_reverse_af641#103 - GLCpts#146
                                 + GLCpts_reverse_a52ae#147  = 0
 gln__L_c#1:                     - 0.2557 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.2557 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + GLNS#148 - GLNS_reverse_59581#149
                                 + GLNabc#150 - GLNabc_reverse_1d82a#151
                                 - GLUN#154 + GLUN_reverse_4ccdb#155
                                 - GLUSy#156 + GLUSy_reverse_6a00f#157  = 0
 gln__L_e#2:                     - _EX_gln__L_e#104
                                 + _EX_gln__L_e_reverse_6a1a1#105 - GLNabc#150
                                 + GLNabc_reverse_1d82a#151  = 0
 glu__L_c#3:                     - 4.9414 BIOMASS_Ecoli_core_w_GAM#48
                                 + 4.9414 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - GLNS#148 + GLNS_reverse_59581#149
                                 - GLUDy#152 + GLUDy_reverse_fa4e7#153
                                 + GLUN#154 - GLUN_reverse_4ccdb#155
                                 + 2 GLUSy#156 - 2 GLUSy_reverse_6a00f#157
                                 + GLUt2r#158 - GLUt2r_reverse_3e88a#159  = 0
 glu__L_e#4:                     - _EX_glu__L_e#106
                                 + _EX_glu__L_e_reverse_42f6c#107 - GLUt2r#158
                                 + GLUt2r_reverse_3e88a#159  = 0
 glx_c#5:                        ICL#166 - ICL_reverse_2f27e#167 - MALS#170
                                 + MALS_reverse_d7382#171  = 0
 h2o_c#6:                        - PGL#8 + PGL_reverse_2bb6b#9 - PPC#24
                                 + PPC_reverse_e854a#25 + ACONTa#26
                                 - ACONTa_reverse_cad6d#27 - ACONTb#28
                                 + ACONTb_reverse_e198a#29 - ATPM#30
                                 + ATPM_reverse_5b752#31 - PPS#36
                                 + PPS_reverse_1c319#37 + ATPS4r#42
                                 - ATPS4r_reverse_64306#43
                                 - 59.81 BIOMASS_Ecoli_core_w_GAM#48
                                 + 59.81 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - CS#56 + CS_reverse_8d7e9#57 + CYTBD#62
                                 - CYTBD_reverse_611ba#63 + _ENO#66
                                 - _ENO_reverse_40eea#67 - FBP#128
                                 + FBP_reverse_bf2c9#129 - FUM#138
                                 + FUM_reverse_d3642#139 - GLNabc#150
                                 + GLNabc_reverse_1d82a#151 - GLUDy#152
                                 + GLUDy_reverse_fa4e7#153 - GLUN#154
                                 + GLUN_reverse_4ccdb#155 + H2Ot#162
                                 - H2Ot_reverse_aa560#163 - MALS#170
                                 + MALS_reverse_d7382#171  = 0
 h2o_e#7:                        - _EX_h2o_e#110 + _EX_h2o_e_reverse_3ced4#111
                                 - H2Ot#162 + H2Ot_reverse_aa560#163  = 0
 h_c#8:                          PFK#0 - PFK_reverse_d24a6#1 + PGL#8
                                 - PGL_reverse_2bb6b#9 + ACALD#10
                                 - ACALD_reverse_fda2b#11 + AKGt2r#12
                                 - AKGt2r_reverse_5d500#13 + PIt2r#16
                                 - PIt2r_reverse_1cd61#17 + ALCD2x#18
                                 - ALCD2x_reverse_5d107#19 + PPC#24
                                 - PPC_reverse_e854a#25 + ATPM#30
                                 - ATPM_reverse_5b752#31 + ACt2r#34
                                 - ACt2r_reverse_64e05#35 + 2 PPS#36
                                 - 2 PPS_reverse_1c319#37 + 3 ATPS4r#42
                                 - 3 ATPS4r_reverse_64306#43 - PYK#46
                                 + PYK_reverse_bc8ff#47
                                 + 59.81 BIOMASS_Ecoli_core_w_GAM#48
                                 - 59.81 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + PYRt2#50 - PYRt2_reverse_a6e24#51 + CS#56
                                 - CS_reverse_8d7e9#57 + 2 SUCCt2_2#60
                                 - 2 SUCCt2_2_reverse_206ca#61 - 2 CYTBD#62
                                 + 2 CYTBD_reverse_611ba#63 + D_LACt2#64
                                 - D_LACt2_reverse_d0dc8#65 + SUCCt3#68
                                 - SUCCt3_reverse_7d82a#69 + _ETOHt2r#70
                                 - _ETOHt2r_reverse_b3d52#71 + 2 THD2#78
                                 - 2 THD2_reverse_f65dd#79 + FORt2#130
                                 - FORt2_reverse_89839#131 + 2 FUMt2_2#140
                                 - 2 FUMt2_2_reverse_7b4b6#141 + G6PDH2r#142
                                 - G6PDH2r_reverse_19ddf#143 + GAPD#144
                                 - GAPD_reverse_459c1#145 + GLNS#148
                                 - GLNS_reverse_59581#149 + GLNabc#150
                                 - GLNabc_reverse_1d82a#151 + GLUDy#152
                                 - GLUDy_reverse_fa4e7#153 - GLUSy#156
                                 + GLUSy_reverse_6a00f#157 + GLUt2r#158
                                 - GLUt2r_reverse_3e88a#159 + LDH_D#168
                                 - LDH_D_reverse_f8507#169 + MALS#170
                                 - MALS_reverse_d7382#171 + 2 MALt2_2#172
                                 - 2 MALt2_2_reverse_a635f#173 + MDH#174
                                 - MDH_reverse_ee52c#175 - 4 NADH16#180
                                 + 4 NADH16_reverse_330a2#181  = 0
 h_e#9:                          - AKGt2r#12 + AKGt2r_reverse_5d500#13
                                 - PIt2r#16 + PIt2r_reverse_1cd61#17 - ACt2r#34
                                 + ACt2r_reverse_64e05#35 - 4 ATPS4r#42
                                 + 4 ATPS4r_reverse_64306#43 - PYRt2#50
                                 + PYRt2_reverse_a6e24#51 - 2 SUCCt2_2#60
                                 + 2 SUCCt2_2_reverse_206ca#61 + 2 CYTBD#62
                                 - 2 CYTBD_reverse_611ba#63 - D_LACt2#64
                                 + D_LACt2_reverse_d0dc8#65 - SUCCt3#68
                                 + SUCCt3_reverse_7d82a#69 - _ETOHt2r#70
                                 + _ETOHt2r_reverse_b3d52#71 - 2 THD2#78
                                 + 2 THD2_reverse_f65dd#79 - _EX_h_e#108
                                 + _EX_h_e_reverse_3e0c5#109 - FORt2#130
                                 + FORt2_reverse_89839#131 - 2 FUMt2_2#140
                                 + 2 FUMt2_2_reverse_7b4b6#141 - GLUt2r#158
                                 + GLUt2r_reverse_3e88a#159 - 2 MALt2_2#172
                                 + 2 MALt2_2_reverse_a635f#173 + 3 NADH16#180
                                 - 3 NADH16_reverse_330a2#181  = 0
 icit_c#10:                      ACONTb#28 - ACONTb_reverse_e198a#29
                                 - ICDHyr#164 + ICDHyr_reverse_7f84b#165
                                 - ICL#166 + ICL_reverse_2f27e#167  = 0
 lac__D_c#11:                    D_LACt2#64 - D_LACt2_reverse_d0dc8#65
                                 - LDH_D#168 + LDH_D_reverse_f8507#169  = 0
 lac__D_e#12:                    - D_LACt2#64 + D_LACt2_reverse_d0dc8#65
                                 - _EX_lac__D_e#112
                                 + _EX_lac__D_e_reverse_f95b4#113  = 0
 mal__L_c#13:                    FUM#138 - FUM_reverse_d3642#139 + MALS#170
                                 - MALS_reverse_d7382#171 + MALt2_2#172
                                 - MALt2_2_reverse_a635f#173 - MDH#174
                                 + MDH_reverse_ee52c#175 - ME1#176
                                 + ME1_reverse_9736c#177 - ME2#178
                                 + ME2_reverse_2b0a2#179  = 0
 mal__L_e#14:                    - _EX_mal__L_e#114
                                 + _EX_mal__L_e_reverse_af154#115 - MALt2_2#172
                                 + MALt2_2_reverse_a635f#173  = 0
 nad_c#15:                       - ACALD#10 + ACALD_reverse_fda2b#11
                                 - ALCD2x#18 + ALCD2x_reverse_5d107#19
                                 - AKGDH#40 + AKGDH_reverse_08bdc#41
                                 - 3.547 BIOMASS_Ecoli_core_w_GAM#48
                                 + 3.547 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + THD2#78 - THD2_reverse_f65dd#79 - GAPD#144
                                 + GAPD_reverse_459c1#145 - LDH_D#168
                                 + LDH_D_reverse_f8507#169 - MDH#174
                                 + MDH_reverse_ee52c#175 - ME1#176
                                 + ME1_reverse_9736c#177 + NADH16#180
                                 - NADH16_reverse_330a2#181 - NADTRHD#182
                                 + NADTRHD_reverse_49725#183 - PDH#188
                                 + PDH_reverse_ca160#189  = 0
 nadh_c#16:                      ACALD#10 - ACALD_reverse_fda2b#11 + ALCD2x#18
                                 - ALCD2x_reverse_5d107#19 + AKGDH#40
                                 - AKGDH_reverse_08bdc#41
                                 + 3.547 BIOMASS_Ecoli_core_w_GAM#48
                                 - 3.547 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - THD2#78 + THD2_reverse_f65dd#79 + GAPD#144
                                 - GAPD_reverse_459c1#145 + LDH_D#168
                                 - LDH_D_reverse_f8507#169 + MDH#174
                                 - MDH_reverse_ee52c#175 + ME1#176
                                 - ME1_reverse_9736c#177 - NADH16#180
                                 + NADH16_reverse_330a2#181 + NADTRHD#182
                                 - NADTRHD_reverse_49725#183 + PDH#188
                                 - PDH_reverse_ca160#189  = 0
 nadp_c#17:                      13.0279 BIOMASS_Ecoli_core_w_GAM#48
                                 - 13.0279 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - THD2#78 + THD2_reverse_f65dd#79
                                 - G6PDH2r#142 + G6PDH2r_reverse_19ddf#143
                                 - GLUDy#152 + GLUDy_reverse_fa4e7#153
                                 + GLUSy#156 - GLUSy_reverse_6a00f#157
                                 - GND#160 + GND_reverse_eec5c#161 - ICDHyr#164
                                 + ICDHyr_reverse_7f84b#165 - ME2#178
                                 + ME2_reverse_2b0a2#179 + NADTRHD#182
                                 - NADTRHD_reverse_49725#183  = 0
 nadph_c#18:                     - 13.0279 BIOMASS_Ecoli_core_w_GAM#48
                                 + 13.0279 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + THD2#78 - THD2_reverse_f65dd#79
                                 + G6PDH2r#142 - G6PDH2r_reverse_19ddf#143
                                 + GLUDy#152 - GLUDy_reverse_fa4e7#153
                                 - GLUSy#156 + GLUSy_reverse_6a00f#157
                                 + GND#160 - GND_reverse_eec5c#161 + ICDHyr#164
                                 - ICDHyr_reverse_7f84b#165 + ME2#178
                                 - ME2_reverse_2b0a2#179 - NADTRHD#182
                                 + NADTRHD_reverse_49725#183  = 0
 nh4_c#19:                       - GLNS#148 + GLNS_reverse_59581#149
                                 + GLUDy#152 - GLUDy_reverse_fa4e7#153
                                 + GLUN#154 - GLUN_reverse_4ccdb#155 + NH4t#184
                                 - NH4t_reverse_551ee#185  = 0
 _13dpg_c#20:                    PGK#6 - PGK_reverse_02696#7 + GAPD#144
                                 - GAPD_reverse_459c1#145  = 0
 nh4_e#21:                       - _EX_nh4_e#116 + _EX_nh4_e_reverse_f9cc6#117
                                 - NH4t#184 + NH4t_reverse_551ee#185  = 0
 o2_c#22:                        - 0.5 CYTBD#62 + 0.5 CYTBD_reverse_611ba#63
                                 + O2t#186 - O2t_reverse_4d957#187  = 0
 _2pg_c#23:                      - PGM#14 + PGM_reverse_fc9af#15 - _ENO#66
                                 + _ENO_reverse_40eea#67  = 0
 o2_e#24:                        - _EX_o2_e#118 + _EX_o2_e_reverse_efa94#119
                                 - O2t#186 + O2t_reverse_4d957#187  = 0
 _3pg_c#25:                      - PGK#6 + PGK_reverse_02696#7 + PGM#14
                                 - PGM_reverse_fc9af#15
                                 - 1.496 BIOMASS_Ecoli_core_w_GAM#48
                                 + 1.496 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                  = 0
 oaa_c#26:                       PPC#24 - PPC_reverse_e854a#25 - PPCK#32
                                 + PPCK_reverse_2557d#33
                                 - 1.7867 BIOMASS_Ecoli_core_w_GAM#48
                                 + 1.7867 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - CS#56 + CS_reverse_8d7e9#57 + MDH#174
                                 - MDH_reverse_ee52c#175  = 0
 pep_c#27:                       - PPC#24 + PPC_reverse_e854a#25 + PPCK#32
                                 - PPCK_reverse_2557d#33 + PPS#36
                                 - PPS_reverse_1c319#37 - PYK#46
                                 + PYK_reverse_bc8ff#47
                                 - 0.5191 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.5191 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + _ENO#66 - _ENO_reverse_40eea#67
                                 - FRUpts2#136 + FRUpts2_reverse_58d3a#137
                                 - GLCpts#146 + GLCpts_reverse_a52ae#147  = 0
 _6pgc_c#28:                     PGL#8 - PGL_reverse_2bb6b#9 - GND#160
                                 + GND_reverse_eec5c#161  = 0
 pi_c#29:                        PIt2r#16 - PIt2r_reverse_1cd61#17 + PPC#24
                                 - PPC_reverse_e854a#25 + ATPM#30
                                 - ATPM_reverse_5b752#31 + PPS#36
                                 - PPS_reverse_1c319#37 - ATPS4r#42
                                 + ATPS4r_reverse_64306#43 - PTAr#44
                                 + PTAr_reverse_fce15#45
                                 + 59.81 BIOMASS_Ecoli_core_w_GAM#48
                                 - 59.81 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + SUCOAS#74 - SUCOAS_reverse_22958#75
                                 + FBP#128 - FBP_reverse_bf2c9#129 - GAPD#144
                                 + GAPD_reverse_459c1#145 + GLNS#148
                                 - GLNS_reverse_59581#149 + GLNabc#150
                                 - GLNabc_reverse_1d82a#151  = 0
 _6pgl_c#30:                     - PGL#8 + PGL_reverse_2bb6b#9 + G6PDH2r#142
                                 - G6PDH2r_reverse_19ddf#143  = 0
 pi_e#31:                        - PIt2r#16 + PIt2r_reverse_1cd61#17
                                 - _EX_pi_e#120 + _EX_pi_e_reverse_1fb09#121
                                  = 0
 ac_c#32:                        - ACKr#22 + ACKr_reverse_b49c0#23 + ACt2r#34
                                 - ACt2r_reverse_64e05#35  = 0
 pyr_c#33:                       - PFL#2 + PFL_reverse_af9ec#3 - PPS#36
                                 + PPS_reverse_1c319#37 + PYK#46
                                 - PYK_reverse_bc8ff#47
                                 - 2.8328 BIOMASS_Ecoli_core_w_GAM#48
                                 + 2.8328 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + PYRt2#50 - PYRt2_reverse_a6e24#51
                                 + FRUpts2#136 - FRUpts2_reverse_58d3a#137
                                 + GLCpts#146 - GLCpts_reverse_a52ae#147
                                 + LDH_D#168 - LDH_D_reverse_f8507#169
                                 + ME1#176 - ME1_reverse_9736c#177 + ME2#178
                                 - ME2_reverse_2b0a2#179 - PDH#188
                                 + PDH_reverse_ca160#189  = 0
 pyr_e#34:                       - PYRt2#50 + PYRt2_reverse_a6e24#51
                                 - _EX_pyr_e#122 + _EX_pyr_e_reverse_1f6de#123
                                  = 0
 q8_c#35:                        CYTBD#62 - CYTBD_reverse_611ba#63 - SUCDi#72
                                 + SUCDi_reverse_480f4#73 + FRD7#134
                                 - FRD7_reverse_e0172#135 - NADH16#180
                                 + NADH16_reverse_330a2#181  = 0
 q8h2_c#36:                      - CYTBD#62 + CYTBD_reverse_611ba#63 + SUCDi#72
                                 - SUCDi_reverse_480f4#73 - FRD7#134
                                 + FRD7_reverse_e0172#135 + NADH16#180
                                 - NADH16_reverse_330a2#181  = 0
 r5p_c#37:                       - 0.8977 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.8977 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - RPI#58 + RPI_reverse_853a1#59 - TKT1#80
                                 + TKT1_reverse_a1021#81  = 0
 ru5p__D_c#38:                   - RPE#54 + RPE_reverse_a1b04#55 + RPI#58
                                 - RPI_reverse_853a1#59 + GND#160
                                 - GND_reverse_eec5c#161  = 0
 ac_e#39:                        - ACt2r#34 + ACt2r_reverse_64e05#35
                                 - _EX_ac_e#86 + _EX_ac_e_reverse_0be96#87  = 0
 acald_c#40:                     - ACALD#10 + ACALD_reverse_fda2b#11
                                 + ALCD2x#18 - ALCD2x_reverse_5d107#19
                                 + ACALDt#20 - ACALDt_reverse_858fa#21  = 0
 s7p_c#41:                       - TALA#76 + TALA_reverse_adfda#77 + TKT1#80
                                 - TKT1_reverse_a1021#81  = 0
 acald_e#42:                     - ACALDt#20 + ACALDt_reverse_858fa#21
                                 - _EX_acald_e#88
                                 + _EX_acald_e_reverse_c096e#89  = 0
 accoa_c#43:                     PFL#2 - PFL_reverse_af9ec#3 + ACALD#10
                                 - ACALD_reverse_fda2b#11 - PTAr#44
                                 + PTAr_reverse_fce15#45
                                 - 3.7478 BIOMASS_Ecoli_core_w_GAM#48
                                 + 3.7478 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - CS#56 + CS_reverse_8d7e9#57 - MALS#170
                                 + MALS_reverse_d7382#171 + PDH#188
                                 - PDH_reverse_ca160#189  = 0
 succ_c#44:                      SUCCt2_2#60 - SUCCt2_2_reverse_206ca#61
                                 - SUCCt3#68 + SUCCt3_reverse_7d82a#69
                                 - SUCDi#72 + SUCDi_reverse_480f4#73
                                 - SUCOAS#74 + SUCOAS_reverse_22958#75
                                 + FRD7#134 - FRD7_reverse_e0172#135 + ICL#166
                                 - ICL_reverse_2f27e#167  = 0
 succ_e#45:                      - SUCCt2_2#60 + SUCCt2_2_reverse_206ca#61
                                 + SUCCt3#68 - SUCCt3_reverse_7d82a#69
                                 - _EX_succ_e#124
                                 + _EX_succ_e_reverse_a9039#125  = 0
 succoa_c#46:                    AKGDH#40 - AKGDH_reverse_08bdc#41 + SUCOAS#74
                                 - SUCOAS_reverse_22958#75  = 0
 acon_C_c#47:                    ACONTa#26 - ACONTa_reverse_cad6d#27
                                 - ACONTb#28 + ACONTb_reverse_e198a#29  = 0
 xu5p__D_c#48:                   RPE#54 - RPE_reverse_a1b04#55 - TKT1#80
                                 + TKT1_reverse_a1021#81 - TKT2#82
                                 + TKT2_reverse_7ebc7#83  = 0
 actp_c#49:                      ACKr#22 - ACKr_reverse_b49c0#23 + PTAr#44
                                 - PTAr_reverse_fce15#45  = 0
 adp_c#50:                       PFK#0 - PFK_reverse_d24a6#1 + PGK#6
                                 - PGK_reverse_02696#7 + ACKr#22
                                 - ACKr_reverse_b49c0#23 + ATPM#30
                                 - ATPM_reverse_5b752#31 + PPCK#32
                                 - PPCK_reverse_2557d#33 + 2 ADK1#38
                                 - 2 ADK1_reverse_a6f90#39 - ATPS4r#42
                                 + ATPS4r_reverse_64306#43 - PYK#46
                                 + PYK_reverse_bc8ff#47
                                 + 59.81 BIOMASS_Ecoli_core_w_GAM#48
                                 - 59.81 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + SUCOAS#74 - SUCOAS_reverse_22958#75
                                 + GLNS#148 - GLNS_reverse_59581#149
                                 + GLNabc#150 - GLNabc_reverse_1d82a#151  = 0
 akg_c#51:                       AKGt2r#12 - AKGt2r_reverse_5d500#13 - AKGDH#40
                                 + AKGDH_reverse_08bdc#41
                                 + 4.1182 BIOMASS_Ecoli_core_w_GAM#48
                                 - 4.1182 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + GLUDy#152 - GLUDy_reverse_fa4e7#153
                                 - GLUSy#156 + GLUSy_reverse_6a00f#157
                                 + ICDHyr#164 - ICDHyr_reverse_7f84b#165  = 0
 akg_e#52:                       - AKGt2r#12 + AKGt2r_reverse_5d500#13
                                 - _EX_akg_e#90 + _EX_akg_e_reverse_70d85#91
                                  = 0
 amp_c#53:                       PPS#36 - PPS_reverse_1c319#37 - ADK1#38
                                 + ADK1_reverse_a6f90#39  = 0
 atp_c#54:                       - PFK#0 + PFK_reverse_d24a6#1 - PGK#6
                                 + PGK_reverse_02696#7 - ACKr#22
                                 + ACKr_reverse_b49c0#23 - ATPM#30
                                 + ATPM_reverse_5b752#31 - PPCK#32
                                 + PPCK_reverse_2557d#33 - PPS#36
                                 + PPS_reverse_1c319#37 - ADK1#38
                                 + ADK1_reverse_a6f90#39 + ATPS4r#42
                                 - ATPS4r_reverse_64306#43 + PYK#46
                                 - PYK_reverse_bc8ff#47
                                 - 59.81 BIOMASS_Ecoli_core_w_GAM#48
                                 + 59.81 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - SUCOAS#74 + SUCOAS_reverse_22958#75
                                 - GLNS#148 + GLNS_reverse_59581#149
                                 - GLNabc#150 + GLNabc_reverse_1d82a#151  = 0
 cit_c#55:                       - ACONTa#26 + ACONTa_reverse_cad6d#27 + CS#56
                                 - CS_reverse_8d7e9#57  = 0
 co2_c#56:                       - PPC#24 + PPC_reverse_e854a#25 + PPCK#32
                                 - PPCK_reverse_2557d#33 + AKGDH#40
                                 - AKGDH_reverse_08bdc#41 + CO2t#52
                                 - CO2t_reverse_7c42f#53 + GND#160
                                 - GND_reverse_eec5c#161 + ICDHyr#164
                                 - ICDHyr_reverse_7f84b#165 + ME1#176
                                 - ME1_reverse_9736c#177 + ME2#178
                                 - ME2_reverse_2b0a2#179 + PDH#188
                                 - PDH_reverse_ca160#189  = 0
 co2_e#57:                       - CO2t#52 + CO2t_reverse_7c42f#53
                                 - _EX_co2_e#92 + _EX_co2_e_reverse_d0466#93
                                  = 0
 coa_c#58:                       - PFL#2 + PFL_reverse_af9ec#3 - ACALD#10
                                 + ACALD_reverse_fda2b#11 - AKGDH#40
                                 + AKGDH_reverse_08bdc#41 + PTAr#44
                                 - PTAr_reverse_fce15#45
                                 + 3.7478 BIOMASS_Ecoli_core_w_GAM#48
                                 - 3.7478 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + CS#56 - CS_reverse_8d7e9#57 - SUCOAS#74
                                 + SUCOAS_reverse_22958#75 + MALS#170
                                 - MALS_reverse_d7382#171 - PDH#188
                                 + PDH_reverse_ca160#189  = 0
 dhap_c#59:                      - TPI#84 + TPI_reverse_c2c3b#85 + FBA#126
                                 - FBA_reverse_84806#127  = 0
 _e4p_c#60:                      - 0.361 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.361 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + TALA#76 - TALA_reverse_adfda#77 - TKT2#82
                                 + TKT2_reverse_7ebc7#83  = 0
 _etoh_c#61:                     - ALCD2x#18 + ALCD2x_reverse_5d107#19
                                 + _ETOHt2r#70 - _ETOHt2r_reverse_b3d52#71  = 0
 _etoh_e#62:                     - _ETOHt2r#70 + _ETOHt2r_reverse_b3d52#71
                                 - _EX_etoh_e#94 + _EX_etoh_e_reverse_cc64f#95
                                  = 0
 f6p_c#63:                       - PFK#0 + PFK_reverse_d24a6#1 + PGI#4
                                 - PGI_reverse_27efc#5
                                 - 0.0709 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.0709 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + TALA#76 - TALA_reverse_adfda#77 + TKT2#82
                                 - TKT2_reverse_7ebc7#83 + FBP#128
                                 - FBP_reverse_bf2c9#129 + FRUpts2#136
                                 - FRUpts2_reverse_58d3a#137  = 0
 fdp_c#64:                       PFK#0 - PFK_reverse_d24a6#1 - FBA#126
                                 + FBA_reverse_84806#127 - FBP#128
                                 + FBP_reverse_bf2c9#129  = 0
 for_c#65:                       PFL#2 - PFL_reverse_af9ec#3 + FORt2#130
                                 - FORt2_reverse_89839#131 + FORt#132
                                 - FORt_reverse_40f9f#133  = 0
 for_e#66:                       - _EX_for_e#96 + _EX_for_e_reverse_23269#97
                                 - FORt2#130 + FORt2_reverse_89839#131
                                 - FORt#132 + FORt_reverse_40f9f#133  = 0
 fru_e#67:                       - _EX_fru_e#98 + _EX_fru_e_reverse_c3828#99
                                 - FRUpts2#136 + FRUpts2_reverse_58d3a#137  = 0
 fum_c#68:                       SUCDi#72 - SUCDi_reverse_480f4#73 - FRD7#134
                                 + FRD7_reverse_e0172#135 - FUM#138
                                 + FUM_reverse_d3642#139 + FUMt2_2#140
                                 - FUMt2_2_reverse_7b4b6#141  = 0
 fum_e#69:                       - _EX_fum_e#100 + _EX_fum_e_reverse_e3432#101
                                 - FUMt2_2#140 + FUMt2_2_reverse_7b4b6#141  = 0
 g3p_c#70:                       - 0.129 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.129 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - TALA#76 + TALA_reverse_adfda#77 + TKT1#80
                                 - TKT1_reverse_a1021#81 + TKT2#82
                                 - TKT2_reverse_7ebc7#83 + TPI#84
                                 - TPI_reverse_c2c3b#85 + FBA#126
                                 - FBA_reverse_84806#127 - GAPD#144
                                 + GAPD_reverse_459c1#145  = 0
 g6p_c#71:                       - PGI#4 + PGI_reverse_27efc#5
                                 - 0.205 BIOMASS_Ecoli_core_w_GAM#48
                                 + 0.205 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 - G6PDH2r#142 + G6PDH2r_reverse_19ddf#143
                                 + GLCpts#146 - GLCpts_reverse_a52ae#147  = 0
 PFK_revbinF#72:                 PFK#0 - 1000 PFK_revbin#195 <= 0
 PFK_revbinR#73:                 PFK_reverse_d24a6#1 + 1000 PFK_revbin#195
                                 <= 1000
 PFK_full_thermo#74:             adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 - f6p_c_ln_concentration#192
                                 + fdp_c_ln_concentration#193
                                 - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PFK_revbin#195  = 1000
 PFK_simple_thermo#75:           1000 PFK_revbin#195 + adp_c_potential#197
                                 - atp_c_potential#198 - f6p_c_potential#199
                                 + fdp_c_potential#200 + h_c_potential#201
                                 - RgPFK_simple_thermo#75  = 0
 PFL_revbinF#76:                 PFL#2 - 1000 PFL_revbin#202 <= 0
 PFL_revbinR#77:                 PFL_reverse_af9ec#3 + 1000 PFL_revbin#202
                                 <= 1000
 PFL_simple_thermo#78:           1000 PFL_revbin#202 + accoa_c_potential#204
                                 - coa_c_potential#205 + for_c_potential#206
                                 - pyr_c_potential#207 - RgPFL_simple_thermo#78
                                  = 0
 PGI_revbinF#79:                 PGI#4 - 1000 PGI_revbin#208 <= 0
 PGI_revbinR#80:                 PGI_reverse_27efc#5 + 1000 PGI_revbin#208
                                 <= 1000
 PGI_simple_thermo#81:           f6p_c_potential#199 + 1000 PGI_revbin#208
                                 - g6p_c_potential#210 - RgPGI_simple_thermo#81
                                  = 0
 PGK_revbinF#82:                 PGK#6 - 1000 PGK_revbin#211 <= 0
 PGK_revbinR#83:                 PGK_reverse_02696#7 + 1000 PGK_revbin#211
                                 <= 1000
 PGK_simple_thermo#84:           adp_c_potential#197 - atp_c_potential#198
                                 + 1000 PGK_revbin#211 + _13dpg_c_potential#213
                                 - _3pg_c_potential#214
                                 - RgPGK_simple_thermo#84  = 0
 PGL_revbinF#85:                 PGL#8 - 1000 PGL_revbin#215 <= 0
 PGL_revbinR#86:                 PGL_reverse_2bb6b#9 + 1000 PGL_revbin#215
                                 <= 1000
 PGL_simple_thermo#87:           h_c_potential#201 + 1000 PGL_revbin#215
                                 + _6pgc_c_potential#217
                                 - _6pgl_c_potential#218 - h2o_c_potential#219
                                 - RgPGL_simple_thermo#87  = 0
 ACALD_revbinF#88:               ACALD#10 - 1000 ACALD_revbin#220 <= 0
 ACALD_revbinR#89:               ACALD_reverse_fda2b#11 + 1000 ACALD_revbin#220
                                 <= 1000
 ACALD_simple_thermo#90:         h_c_potential#201 + accoa_c_potential#204
                                 - coa_c_potential#205 + 1000 ACALD_revbin#220
                                 - acald_c_potential#222 - nad_c_potential#223
                                 + nadh_c_potential#224
                                 - RgACALD_simple_thermo#90  = 0
 AKGt2r_revbinF#91:              AKGt2r#12 - 1000 AKGt2r_revbin#225 <= 0
 AKGt2r_revbinR#92:              AKGt2r_reverse_5d500#13
                                 + 1000 AKGt2r_revbin#225 <= 1000
 AKGt2r_simple_thermo#93:        h_c_potential#201 + 1000 AKGt2r_revbin#225
                                 + akg_c_potential#227 - akg_e_potential#228
                                 - h_e_potential#229
                                 - RgAKGt2r_simple_thermo#93  = 0
 PGM_revbinF#94:                 PGM#14 - 1000 PGM_revbin#230 <= 0
 PGM_revbinR#95:                 PGM_reverse_fc9af#15 + 1000 PGM_revbin#230
                                 <= 1000
 PGM_simple_thermo#96:           _3pg_c_potential#214 + 1000 PGM_revbin#230
                                 - _2pg_c_potential#232
                                 - RgPGM_simple_thermo#96  = 0
 PIt2r_revbinF#97:               PIt2r#16 - 1000 PIt2r_revbin#233 <= 0
 PIt2r_revbinR#98:               PIt2r_reverse_1cd61#17 + 1000 PIt2r_revbin#233
                                 <= 1000
 PIt2r_simple_thermo#99:         h_c_potential#201 - h_e_potential#229
                                 + 1000 PIt2r_revbin#233 + pi_c_potential#235
                                 - pi_e_potential#236
                                 - RgPIt2r_simple_thermo#99  = 0
 ALCD2x_revbinF#100:             ALCD2x#18 - 1000 ALCD2x_revbin#237 <= 0
 ALCD2x_revbinR#101:             ALCD2x_reverse_5d107#19
                                 + 1000 ALCD2x_revbin#237 <= 1000
 ALCD2x_simple_thermo#102:       h_c_potential#201 + acald_c_potential#222
                                 - nad_c_potential#223 + nadh_c_potential#224
                                 + 1000 ALCD2x_revbin#237
                                 - _etoh_c_potential#239
                                 - RgALCD2x_simple_thermo#102  = 0
 ACALDt_revbinF#103:             ACALDt#20 - 1000 ACALDt_revbin#240 <= 0
 ACALDt_revbinR#104:             ACALDt_reverse_858fa#21
                                 + 1000 ACALDt_revbin#240 <= 1000
 ACALDt_simple_thermo#105:       acald_c_potential#222 + 1000 ACALDt_revbin#240
                                 - acald_e_potential#242
                                 - RgACALDt_simple_thermo#105  = 0
 ACKr_revbinF#106:               ACKr#22 - 1000 ACKr_revbin#243 <= 0
 ACKr_revbinR#107:               ACKr_reverse_b49c0#23 + 1000 ACKr_revbin#243
                                 <= 1000
 ACKr_simple_thermo#108:         adp_c_potential#197 - atp_c_potential#198
                                 + 1000 ACKr_revbin#243 - ac_c_potential#245
                                 + actp_c_potential#246
                                 - RgACKr_simple_thermo#108  = 0
 PPC_revbinF#109:                PPC#24 - 1000 PPC_revbin#247 <= 0
 PPC_revbinR#110:                PPC_reverse_e854a#25 + 1000 PPC_revbin#247
                                 <= 1000
 PPC_simple_thermo#111:          h_c_potential#201 - h2o_c_potential#219
                                 + pi_c_potential#235 + 1000 PPC_revbin#247
                                 - co2_c_potential#249 + oaa_c_potential#250
                                 - pep_c_potential#251
                                 - RgPPC_simple_thermo#111  = 0
 ACONTa_revbinF#112:             ACONTa#26 - 1000 ACONTa_revbin#252 <= 0
 ACONTa_revbinR#113:             ACONTa_reverse_cad6d#27
                                 + 1000 ACONTa_revbin#252 <= 1000
 ACONTa_simple_thermo#114:       h2o_c_potential#219 + 1000 ACONTa_revbin#252
                                 + acon_C_c_potential#254 - cit_c_potential#255
                                 - RgACONTa_simple_thermo#114  = 0
 ACONTb_revbinF#115:             ACONTb#28 - 1000 ACONTb_revbin#256 <= 0
 ACONTb_revbinR#116:             ACONTb_reverse_e198a#29
                                 + 1000 ACONTb_revbin#256 <= 1000
 ACONTb_simple_thermo#117:       - h2o_c_potential#219 - acon_C_c_potential#254
                                 + 1000 ACONTb_revbin#256
                                 + icit_c_potential#258
                                 - RgACONTb_simple_thermo#117  = 0
 ATPM_revbinF#118:               ATPM#30 - 1000 ATPM_revbin#259 <= 0
 ATPM_revbinR#119:               ATPM_reverse_5b752#31 + 1000 ATPM_revbin#259
                                 <= 1000
 ATPM_simple_thermo#120:         adp_c_potential#197 - atp_c_potential#198
                                 + h_c_potential#201 - h2o_c_potential#219
                                 + pi_c_potential#235 + 1000 ATPM_revbin#259
                                 - RgATPM_simple_thermo#120  = 0
 PPCK_revbinF#121:               PPCK#32 - 1000 PPCK_revbin#261 <= 0
 PPCK_revbinR#122:               PPCK_reverse_2557d#33 + 1000 PPCK_revbin#261
                                 <= 1000
 PPCK_simple_thermo#123:         adp_c_potential#197 - atp_c_potential#198
                                 + co2_c_potential#249 - oaa_c_potential#250
                                 + pep_c_potential#251 + 1000 PPCK_revbin#261
                                 - RgPPCK_simple_thermo#123  = 0
 ACt2r_revbinF#124:              ACt2r#34 - 1000 ACt2r_revbin#263 <= 0
 ACt2r_revbinR#125:              ACt2r_reverse_64e05#35 + 1000 ACt2r_revbin#263
                                 <= 1000
 ACt2r_simple_thermo#126:        h_c_potential#201 - h_e_potential#229
                                 + ac_c_potential#245 + 1000 ACt2r_revbin#263
                                 - ac_e_potential#265
                                 - RgACt2r_simple_thermo#126  = 0
 PPS_revbinF#127:                PPS#36 - 1000 PPS_revbin#266 <= 0
 PPS_revbinR#128:                PPS_reverse_1c319#37 + 1000 PPS_revbin#266
                                 <= 1000
 PPS_simple_thermo#129:          - atp_c_potential#198 + 2 h_c_potential#201
                                 - pyr_c_potential#207 - h2o_c_potential#219
                                 + pi_c_potential#235 + pep_c_potential#251
                                 + 1000 PPS_revbin#266 + amp_c_potential#268
                                 - RgPPS_simple_thermo#129  = 0
 ADK1_revbinF#130:               ADK1#38 - 1000 ADK1_revbin#269 <= 0
 ADK1_revbinR#131:               ADK1_reverse_a6f90#39 + 1000 ADK1_revbin#269
                                 <= 1000
 ADK1_simple_thermo#132:         2 adp_c_potential#197 - atp_c_potential#198
                                 - amp_c_potential#268 + 1000 ADK1_revbin#269
                                 - RgADK1_simple_thermo#132  = 0
 AKGDH_revbinF#133:              AKGDH#40 - 1000 AKGDH_revbin#271 <= 0
 AKGDH_revbinR#134:              AKGDH_reverse_08bdc#41 + 1000 AKGDH_revbin#271
                                 <= 1000
 AKGDH_simple_thermo#135:        - coa_c_potential#205 - nad_c_potential#223
                                 + nadh_c_potential#224 - akg_c_potential#227
                                 + co2_c_potential#249 + 1000 AKGDH_revbin#271
                                 + succoa_c_potential#273
                                 - RgAKGDH_simple_thermo#135  = 0
 ATPS4r_revbinF#136:             ATPS4r#42 - 1000 ATPS4r_revbin#274 <= 0
 ATPS4r_revbinR#137:             ATPS4r_reverse_64306#43
                                 + 1000 ATPS4r_revbin#274 <= 1000
 ATPS4r_simple_thermo#138:       - adp_c_potential#197 + atp_c_potential#198
                                 + 3 h_c_potential#201 + h2o_c_potential#219
                                 - 4 h_e_potential#229 - pi_c_potential#235
                                 + 1000 ATPS4r_revbin#274
                                 - RgATPS4r_simple_thermo#138  = 0
 PTAr_revbinF#139:               PTAr#44 - 1000 PTAr_revbin#276 <= 0
 PTAr_revbinR#140:               PTAr_reverse_fce15#45 + 1000 PTAr_revbin#276
                                 <= 1000
 PTAr_simple_thermo#141:         - accoa_c_potential#204 + coa_c_potential#205
                                 - pi_c_potential#235 + actp_c_potential#246
                                 + 1000 PTAr_revbin#276
                                 - RgPTAr_simple_thermo#141  = 0
 PYK_revbinF#142:                PYK#46 - 1000 PYK_revbin#278 <= 0
 PYK_revbinR#143:                PYK_reverse_bc8ff#47 + 1000 PYK_revbin#278
                                 <= 1000
 PYK_simple_thermo#144:          - adp_c_potential#197 + atp_c_potential#198
                                 - h_c_potential#201 + pyr_c_potential#207
                                 - pep_c_potential#251 + 1000 PYK_revbin#278
                                 - RgPYK_simple_thermo#144  = 0
 BIOMASS_Ecoli_core_w_GAM_revbinF#145: BIOMASS_Ecoli_core_w_GAM#48
                                 - 1000 BIOMASS_Ecoli_core_w_GAM_revbin#280
                                 <= 0
 BIOMASS_Ecoli_core_w_GAM_revbinR#146:
                                 BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + 1000 BIOMASS_Ecoli_core_w_GAM_revbin#280
                                 <= 1000
 BIOMASS_Ecoli_core_w_GAM_simple_thermo#147: 59.81 adp_c_potential#197
                                 - 59.81 atp_c_potential#198
                                 - 0.0709 f6p_c_potential#199
                                 + 59.81 h_c_potential#201
                                 - 3.7478 accoa_c_potential#204
                                 + 3.7478 coa_c_potential#205
                                 - 2.8328 pyr_c_potential#207
                                 - 0.205 g6p_c_potential#210
                                 - 1.496 _3pg_c_potential#214
                                 - 59.81 h2o_c_potential#219
                                 - 3.547 nad_c_potential#223
                                 + 3.547 nadh_c_potential#224
                                 + 4.1182 akg_c_potential#227
                                 + 59.81 pi_c_potential#235
                                 - 1.7867 oaa_c_potential#250
                                 - 0.5191 pep_c_potential#251
                                 + 1000 BIOMASS_Ecoli_core_w_GAM_revbin#280
                                 - 0.361 _e4p_c_potential#282
                                 - 0.129 g3p_c_potential#283
                                 - 0.2557 gln__L_c_potential#284
                                 - 4.9414 glu__L_c_potential#285
                                 + 13.0279 nadp_c_potential#286
                                 - 13.0279 nadph_c_potential#287
                                 - 0.8977 r5p_c_potential#288
                                 - RgBIOMASS_Ecoli_core_w_GAM_simple_thermo#147
                                  = 0
 PYRt2_revbinF#148:              PYRt2#50 - 1000 PYRt2_revbin#289 <= 0
 PYRt2_revbinR#149:              PYRt2_reverse_a6e24#51 + 1000 PYRt2_revbin#289
                                 <= 1000
 PYRt2_simple_thermo#150:        h_c_potential#201 + pyr_c_potential#207
                                 - h_e_potential#229 + 1000 PYRt2_revbin#289
                                 - pyr_e_potential#291
                                 - RgPYRt2_simple_thermo#150  = 0
 CO2t_revbinF#151:               CO2t#52 - 1000 CO2t_revbin#292 <= 0
 CO2t_revbinR#152:               CO2t_reverse_7c42f#53 + 1000 CO2t_revbin#292
                                 <= 1000
 CO2t_simple_thermo#153:         co2_c_potential#249 + 1000 CO2t_revbin#292
                                 - co2_e_potential#294
                                 - RgCO2t_simple_thermo#153  = 0
 RPE_revbinF#154:                RPE#54 - 1000 RPE_revbin#295 <= 0
 RPE_revbinR#155:                RPE_reverse_a1b04#55 + 1000 RPE_revbin#295
                                 <= 1000
 RPE_simple_thermo#156:          1000 RPE_revbin#295 - ru5p__D_c_potential#297
                                 + xu5p__D_c_potential#298
                                 - RgRPE_simple_thermo#156  = 0
 CS_revbinF#157:                 CS#56 - 1000 CS_revbin#299 <= 0
 CS_revbinR#158:                 CS_reverse_8d7e9#57 + 1000 CS_revbin#299
                                 <= 1000
 CS_simple_thermo#159:           h_c_potential#201 - accoa_c_potential#204
                                 + coa_c_potential#205 - h2o_c_potential#219
                                 - oaa_c_potential#250 + cit_c_potential#255
                                 + 1000 CS_revbin#299 - RgCS_simple_thermo#159
                                  = 0
 RPI_revbinF#160:                RPI#58 - 1000 RPI_revbin#301 <= 0
 RPI_revbinR#161:                RPI_reverse_853a1#59 + 1000 RPI_revbin#301
                                 <= 1000
 RPI_simple_thermo#162:          - r5p_c_potential#288
                                 + ru5p__D_c_potential#297
                                 + 1000 RPI_revbin#301
                                 - RgRPI_simple_thermo#162  = 0
 SUCCt2_2_revbinF#163:           SUCCt2_2#60 - 1000 SUCCt2_2_revbin#303 <= 0
 SUCCt2_2_revbinR#164:           SUCCt2_2_reverse_206ca#61
                                 + 1000 SUCCt2_2_revbin#303 <= 1000
 SUCCt2_2_simple_thermo#165:     2 h_c_potential#201 - 2 h_e_potential#229
                                 + 1000 SUCCt2_2_revbin#303
                                 + succ_c_potential#305 - succ_e_potential#306
                                 - RgSUCCt2_2_simple_thermo#165  = 0
 CYTBD_revbinF#166:              CYTBD#62 - 1000 CYTBD_revbin#307 <= 0
 CYTBD_revbinR#167:              CYTBD_reverse_611ba#63 + 1000 CYTBD_revbin#307
                                 <= 1000
 CYTBD_simple_thermo#168:        - 2 h_c_potential#201 + h2o_c_potential#219
                                 + 2 h_e_potential#229 + 1000 CYTBD_revbin#307
                                 - 0.5 o2_c_potential#309 + q8_c_potential#310
                                 - q8h2_c_potential#311
                                 - RgCYTBD_simple_thermo#168  = 0
 D_LACt2_revbinF#169:            D_LACt2#64 - 1000 D_LACt2_revbin#312 <= 0
 D_LACt2_revbinR#170:            D_LACt2_reverse_d0dc8#65
                                 + 1000 D_LACt2_revbin#312 <= 1000
 D_LACt2_simple_thermo#171:      h_c_potential#201 - h_e_potential#229
                                 + 1000 D_LACt2_revbin#312
                                 + lac__D_c_potential#314
                                 - lac__D_e_potential#315
                                 - RgD_LACt2_simple_thermo#171  = 0
 _ENO_revbinF#172:               _ENO#66 - 1000 _ENO_revbin#316 <= 0
 _ENO_revbinR#173:               _ENO_reverse_40eea#67 + 1000 _ENO_revbin#316
                                 <= 1000
 _ENO_simple_thermo#174:         h2o_c_potential#219 - _2pg_c_potential#232
                                 + pep_c_potential#251 + 1000 _ENO_revbin#316
                                 - Rg_ENO_simple_thermo#174  = 0
 SUCCt3_revbinF#175:             SUCCt3#68 - 1000 SUCCt3_revbin#318 <= 0
 SUCCt3_revbinR#176:             SUCCt3_reverse_7d82a#69
                                 + 1000 SUCCt3_revbin#318 <= 1000
 SUCCt3_simple_thermo#177:       h_c_potential#201 - h_e_potential#229
                                 - succ_c_potential#305 + succ_e_potential#306
                                 + 1000 SUCCt3_revbin#318
                                 - RgSUCCt3_simple_thermo#177  = 0
 _ETOHt2r_revbinF#178:           _ETOHt2r#70 - 1000 _ETOHt2r_revbin#320 <= 0
 _ETOHt2r_revbinR#179:           _ETOHt2r_reverse_b3d52#71
                                 + 1000 _ETOHt2r_revbin#320 <= 1000
 _ETOHt2r_simple_thermo#180:     h_c_potential#201 - h_e_potential#229
                                 + _etoh_c_potential#239
                                 + 1000 _ETOHt2r_revbin#320
                                 - _etoh_e_potential#322
                                 - Rg_ETOHt2r_simple_thermo#180  = 0
 SUCDi_revbinF#181:              SUCDi#72 - 1000 SUCDi_revbin#323 <= 0
 SUCDi_revbinR#182:              SUCDi_reverse_480f4#73 + 1000 SUCDi_revbin#323
                                 <= 1000
 SUCDi_simple_thermo#183:        - succ_c_potential#305 - q8_c_potential#310
                                 + q8h2_c_potential#311 + 1000 SUCDi_revbin#323
                                 + fum_c_potential#325
                                 - RgSUCDi_simple_thermo#183  = 0
 SUCOAS_revbinF#184:             SUCOAS#74 - 1000 SUCOAS_revbin#326 <= 0
 SUCOAS_revbinR#185:             SUCOAS_reverse_22958#75
                                 + 1000 SUCOAS_revbin#326 <= 1000
 SUCOAS_simple_thermo#186:       adp_c_potential#197 - atp_c_potential#198
                                 - coa_c_potential#205 + pi_c_potential#235
                                 + succoa_c_potential#273
                                 - succ_c_potential#305
                                 + 1000 SUCOAS_revbin#326
                                 - RgSUCOAS_simple_thermo#186  = 0
 TALA_revbinF#187:               TALA#76 - 1000 TALA_revbin#328 <= 0
 TALA_revbinR#188:               TALA_reverse_adfda#77 + 1000 TALA_revbin#328
                                 <= 1000
 TALA_simple_thermo#189:         f6p_c_potential#199 + _e4p_c_potential#282
                                 - g3p_c_potential#283 + 1000 TALA_revbin#328
                                 - s7p_c_potential#330
                                 - RgTALA_simple_thermo#189  = 0
 THD2_revbinF#190:               THD2#78 - 1000 THD2_revbin#331 <= 0
 THD2_revbinR#191:               THD2_reverse_f65dd#79 + 1000 THD2_revbin#331
                                 <= 1000
 THD2_simple_thermo#192:         2 h_c_potential#201 + nad_c_potential#223
                                 - nadh_c_potential#224 - 2 h_e_potential#229
                                 - nadp_c_potential#286 + nadph_c_potential#287
                                 + 1000 THD2_revbin#331
                                 - RgTHD2_simple_thermo#192  = 0
 TKT1_revbinF#193:               TKT1#80 - 1000 TKT1_revbin#333 <= 0
 TKT1_revbinR#194:               TKT1_reverse_a1021#81 + 1000 TKT1_revbin#333
                                 <= 1000
 TKT1_simple_thermo#195:         g3p_c_potential#283 - r5p_c_potential#288
                                 - xu5p__D_c_potential#298
                                 + s7p_c_potential#330 + 1000 TKT1_revbin#333
                                 - RgTKT1_simple_thermo#195  = 0
 TKT2_revbinF#196:               TKT2#82 - 1000 TKT2_revbin#335 <= 0
 TKT2_revbinR#197:               TKT2_reverse_7ebc7#83 + 1000 TKT2_revbin#335
                                 <= 1000
 TKT2_simple_thermo#198:         f6p_c_potential#199 - _e4p_c_potential#282
                                 + g3p_c_potential#283
                                 - xu5p__D_c_potential#298
                                 + 1000 TKT2_revbin#335
                                 - RgTKT2_simple_thermo#198  = 0
 TPI_revbinF#199:                TPI#84 - 1000 TPI_revbin#337 <= 0
 TPI_revbinR#200:                TPI_reverse_c2c3b#85 + 1000 TPI_revbin#337
                                 <= 1000
 TPI_simple_thermo#201:          g3p_c_potential#283 + 1000 TPI_revbin#337
                                 - dhap_c_potential#339
                                 - RgTPI_simple_thermo#201  = 0
 _EX_ac_e_revbinF#202:           _EX_ac_e#86 - 1000 _EX_ac_e_revbin#340 <= 0
 _EX_ac_e_revbinR#203:           _EX_ac_e_reverse_0be96#87
                                 + 1000 _EX_ac_e_revbin#340 <= 1000
 _EX_ac_e_simple_thermo#204:     - ac_e_potential#265
                                 + 1000 _EX_ac_e_revbin#340
                                 - Rg_EX_ac_e_simple_thermo#204  = 0
 _EX_acald_e_revbinF#205:        _EX_acald_e#88 - 1000 _EX_acald_e_revbin#342
                                 <= 0
 _EX_acald_e_revbinR#206:        _EX_acald_e_reverse_c096e#89
                                 + 1000 _EX_acald_e_revbin#342 <= 1000
 _EX_acald_e_simple_thermo#207:  - acald_e_potential#242
                                 + 1000 _EX_acald_e_revbin#342
                                 - Rg_EX_acald_e_simple_thermo#207  = 0
 _EX_akg_e_revbinF#208:          _EX_akg_e#90 - 1000 _EX_akg_e_revbin#344 <= 0
 _EX_akg_e_revbinR#209:          _EX_akg_e_reverse_70d85#91
                                 + 1000 _EX_akg_e_revbin#344 <= 1000
 _EX_akg_e_simple_thermo#210:    - akg_e_potential#228
                                 + 1000 _EX_akg_e_revbin#344
                                 - Rg_EX_akg_e_simple_thermo#210  = 0
 _EX_co2_e_revbinF#211:          _EX_co2_e#92 - 1000 _EX_co2_e_revbin#346 <= 0
 _EX_co2_e_revbinR#212:          _EX_co2_e_reverse_d0466#93
                                 + 1000 _EX_co2_e_revbin#346 <= 1000
 _EX_co2_e_simple_thermo#213:    - co2_e_potential#294
                                 + 1000 _EX_co2_e_revbin#346
                                 - Rg_EX_co2_e_simple_thermo#213  = 0
 _EX_etoh_e_revbinF#214:         _EX_etoh_e#94 - 1000 _EX_etoh_e_revbin#348
                                 <= 0
 _EX_etoh_e_revbinR#215:         _EX_etoh_e_reverse_cc64f#95
                                 + 1000 _EX_etoh_e_revbin#348 <= 1000
 _EX_etoh_e_simple_thermo#216:   - _etoh_e_potential#322
                                 + 1000 _EX_etoh_e_revbin#348
                                 - Rg_EX_etoh_e_simple_thermo#216  = 0
 _EX_for_e_revbinF#217:          _EX_for_e#96 - 1000 _EX_for_e_revbin#350 <= 0
 _EX_for_e_revbinR#218:          _EX_for_e_reverse_23269#97
                                 + 1000 _EX_for_e_revbin#350 <= 1000
 _EX_for_e_simple_thermo#219:    1000 _EX_for_e_revbin#350
                                 - for_e_potential#352
                                 - Rg_EX_for_e_simple_thermo#219  = 0
 _EX_fru_e_revbinF#220:          _EX_fru_e#98 - 1000 _EX_fru_e_revbin#353 <= 0
 _EX_fru_e_revbinR#221:          _EX_fru_e_reverse_c3828#99
                                 + 1000 _EX_fru_e_revbin#353 <= 1000
 _EX_fru_e_simple_thermo#222:    1000 _EX_fru_e_revbin#353
                                 - fru_e_potential#355
                                 - Rg_EX_fru_e_simple_thermo#222  = 0
 _EX_fum_e_revbinF#223:          _EX_fum_e#100 - 1000 _EX_fum_e_revbin#356 <= 0
 _EX_fum_e_revbinR#224:          _EX_fum_e_reverse_e3432#101
                                 + 1000 _EX_fum_e_revbin#356 <= 1000
 _EX_fum_e_simple_thermo#225:    1000 _EX_fum_e_revbin#356
                                 - fum_e_potential#358
                                 - Rg_EX_fum_e_simple_thermo#225  = 0
 _EX_glc__D_e_revbinF#226:       _EX_glc__D_e#102
                                 - 1000 _EX_glc__D_e_revbin#359 <= 0
 _EX_glc__D_e_revbinR#227:       _EX_glc__D_e_reverse_af641#103
                                 + 1000 _EX_glc__D_e_revbin#359 <= 1000
 _EX_glc__D_e_simple_thermo#228: 1000 _EX_glc__D_e_revbin#359
                                 - glc__D_e_potential#361
                                 - Rg_EX_glc__D_e_simple_thermo#228  = 0
 _EX_gln__L_e_revbinF#229:       _EX_gln__L_e#104
                                 - 1000 _EX_gln__L_e_revbin#362 <= 0
 _EX_gln__L_e_revbinR#230:       _EX_gln__L_e_reverse_6a1a1#105
                                 + 1000 _EX_gln__L_e_revbin#362 <= 1000
 _EX_gln__L_e_simple_thermo#231: 1000 _EX_gln__L_e_revbin#362
                                 - gln__L_e_potential#364
                                 - Rg_EX_gln__L_e_simple_thermo#231  = 0
 _EX_glu__L_e_revbinF#232:       _EX_glu__L_e#106
                                 - 1000 _EX_glu__L_e_revbin#365 <= 0
 _EX_glu__L_e_revbinR#233:       _EX_glu__L_e_reverse_42f6c#107
                                 + 1000 _EX_glu__L_e_revbin#365 <= 1000
 _EX_glu__L_e_simple_thermo#234: 1000 _EX_glu__L_e_revbin#365
                                 - glu__L_e_potential#367
                                 - Rg_EX_glu__L_e_simple_thermo#234  = 0
 _EX_h_e_revbinF#235:            _EX_h_e#108 - 1000 _EX_h_e_revbin#368 <= 0
 _EX_h_e_revbinR#236:            _EX_h_e_reverse_3e0c5#109
                                 + 1000 _EX_h_e_revbin#368 <= 1000
 _EX_h_e_simple_thermo#237:      - h_e_potential#229 + 1000 _EX_h_e_revbin#368
                                 - Rg_EX_h_e_simple_thermo#237  = 0
 _EX_h2o_e_revbinF#238:          _EX_h2o_e#110 - 1000 _EX_h2o_e_revbin#370 <= 0
 _EX_h2o_e_revbinR#239:          _EX_h2o_e_reverse_3ced4#111
                                 + 1000 _EX_h2o_e_revbin#370 <= 1000
 _EX_h2o_e_simple_thermo#240:    1000 _EX_h2o_e_revbin#370
                                 - h2o_e_potential#372
                                 - Rg_EX_h2o_e_simple_thermo#240  = 0
 _EX_lac__D_e_revbinF#241:       _EX_lac__D_e#112
                                 - 1000 _EX_lac__D_e_revbin#373 <= 0
 _EX_lac__D_e_revbinR#242:       _EX_lac__D_e_reverse_f95b4#113
                                 + 1000 _EX_lac__D_e_revbin#373 <= 1000
 _EX_lac__D_e_simple_thermo#243: - lac__D_e_potential#315
                                 + 1000 _EX_lac__D_e_revbin#373
                                 - Rg_EX_lac__D_e_simple_thermo#243  = 0
 _EX_mal__L_e_revbinF#244:       _EX_mal__L_e#114
                                 - 1000 _EX_mal__L_e_revbin#375 <= 0
 _EX_mal__L_e_revbinR#245:       _EX_mal__L_e_reverse_af154#115
                                 + 1000 _EX_mal__L_e_revbin#375 <= 1000
 _EX_mal__L_e_simple_thermo#246: 1000 _EX_mal__L_e_revbin#375
                                 - mal__L_e_potential#377
                                 - Rg_EX_mal__L_e_simple_thermo#246  = 0
 _EX_nh4_e_revbinF#247:          _EX_nh4_e#116 - 1000 _EX_nh4_e_revbin#378 <= 0
 _EX_nh4_e_revbinR#248:          _EX_nh4_e_reverse_f9cc6#117
                                 + 1000 _EX_nh4_e_revbin#378 <= 1000
 _EX_nh4_e_simple_thermo#249:    1000 _EX_nh4_e_revbin#378
                                 - nh4_e_potential#380
                                 - Rg_EX_nh4_e_simple_thermo#249  = 0
 _EX_o2_e_revbinF#250:           _EX_o2_e#118 - 1000 _EX_o2_e_revbin#381 <= 0
 _EX_o2_e_revbinR#251:           _EX_o2_e_reverse_efa94#119
                                 + 1000 _EX_o2_e_revbin#381 <= 1000
 _EX_o2_e_simple_thermo#252:     1000 _EX_o2_e_revbin#381 - o2_e_potential#383
                                 - Rg_EX_o2_e_simple_thermo#252  = 0
 _EX_pi_e_revbinF#253:           _EX_pi_e#120 - 1000 _EX_pi_e_revbin#384 <= 0
 _EX_pi_e_revbinR#254:           _EX_pi_e_reverse_1fb09#121
                                 + 1000 _EX_pi_e_revbin#384 <= 1000
 _EX_pi_e_simple_thermo#255:     - pi_e_potential#236
                                 + 1000 _EX_pi_e_revbin#384
                                 - Rg_EX_pi_e_simple_thermo#255  = 0
 _EX_pyr_e_revbinF#256:          _EX_pyr_e#122 - 1000 _EX_pyr_e_revbin#386 <= 0
 _EX_pyr_e_revbinR#257:          _EX_pyr_e_reverse_1f6de#123
                                 + 1000 _EX_pyr_e_revbin#386 <= 1000
 _EX_pyr_e_simple_thermo#258:    - pyr_e_potential#291
                                 + 1000 _EX_pyr_e_revbin#386
                                 - Rg_EX_pyr_e_simple_thermo#258  = 0
 _EX_succ_e_revbinF#259:         _EX_succ_e#124 - 1000 _EX_succ_e_revbin#388
                                 <= 0
 _EX_succ_e_revbinR#260:         _EX_succ_e_reverse_a9039#125
                                 + 1000 _EX_succ_e_revbin#388 <= 1000
 _EX_succ_e_simple_thermo#261:   - succ_e_potential#306
                                 + 1000 _EX_succ_e_revbin#388
                                 - Rg_EX_succ_e_simple_thermo#261  = 0
 FBA_revbinF#262:                FBA#126 - 1000 FBA_revbin#390 <= 0
 FBA_revbinR#263:                FBA_reverse_84806#127 + 1000 FBA_revbin#390
                                 <= 1000
 FBA_simple_thermo#264:          - fdp_c_potential#200 + g3p_c_potential#283
                                 + dhap_c_potential#339 + 1000 FBA_revbin#390
                                 - RgFBA_simple_thermo#264  = 0
 FBP_revbinF#265:                FBP#128 - 1000 FBP_revbin#392 <= 0
 FBP_revbinR#266:                FBP_reverse_bf2c9#129 + 1000 FBP_revbin#392
                                 <= 1000
 FBP_simple_thermo#267:          f6p_c_potential#199 - fdp_c_potential#200
                                 - h2o_c_potential#219 + pi_c_potential#235
                                 + 1000 FBP_revbin#392
                                 - RgFBP_simple_thermo#267  = 0
 FORt2_revbinF#268:              FORt2#130 - 1000 FORt2_revbin#394 <= 0
 FORt2_revbinR#269:              FORt2_reverse_89839#131
                                 + 1000 FORt2_revbin#394 <= 1000
 FORt2_simple_thermo#270:        h_c_potential#201 + for_c_potential#206
                                 - h_e_potential#229 - for_e_potential#352
                                 + 1000 FORt2_revbin#394
                                 - RgFORt2_simple_thermo#270  = 0
 FORt_revbinF#271:               FORt#132 - 1000 FORt_revbin#396 <= 0
 FORt_revbinR#272:               FORt_reverse_40f9f#133 + 1000 FORt_revbin#396
                                 <= 1000
 FORt_simple_thermo#273:         for_c_potential#206 - for_e_potential#352
                                 + 1000 FORt_revbin#396
                                 - RgFORt_simple_thermo#273  = 0
 FRD7_revbinF#274:               FRD7#134 - 1000 FRD7_revbin#398 <= 0
 FRD7_revbinR#275:               FRD7_reverse_e0172#135 + 1000 FRD7_revbin#398
                                 <= 1000
 FRD7_simple_thermo#276:         succ_c_potential#305 + q8_c_potential#310
                                 - q8h2_c_potential#311 - fum_c_potential#325
                                 + 1000 FRD7_revbin#398
                                 - RgFRD7_simple_thermo#276  = 0
 FRUpts2_revbinF#277:            FRUpts2#136 - 1000 FRUpts2_revbin#400 <= 0
 FRUpts2_revbinR#278:            FRUpts2_reverse_58d3a#137
                                 + 1000 FRUpts2_revbin#400 <= 1000
 FRUpts2_simple_thermo#279:      f6p_c_potential#199 + pyr_c_potential#207
                                 - pep_c_potential#251 - fru_e_potential#355
                                 + 1000 FRUpts2_revbin#400
                                 - RgFRUpts2_simple_thermo#279  = 0
 FUM_revbinF#280:                FUM#138 - 1000 FUM_revbin#402 <= 0
 FUM_revbinR#281:                FUM_reverse_d3642#139 + 1000 FUM_revbin#402
                                 <= 1000
 FUM_simple_thermo#282:          - h2o_c_potential#219 - fum_c_potential#325
                                 + 1000 FUM_revbin#402 + mal__L_c_potential#404
                                 - RgFUM_simple_thermo#282  = 0
 FUMt2_2_revbinF#283:            FUMt2_2#140 - 1000 FUMt2_2_revbin#405 <= 0
 FUMt2_2_revbinR#284:            FUMt2_2_reverse_7b4b6#141
                                 + 1000 FUMt2_2_revbin#405 <= 1000
 FUMt2_2_simple_thermo#285:      2 h_c_potential#201 - 2 h_e_potential#229
                                 + fum_c_potential#325 - fum_e_potential#358
                                 + 1000 FUMt2_2_revbin#405
                                 - RgFUMt2_2_simple_thermo#285  = 0
 G6PDH2r_revbinF#286:            G6PDH2r#142 - 1000 G6PDH2r_revbin#407 <= 0
 G6PDH2r_revbinR#287:            G6PDH2r_reverse_19ddf#143
                                 + 1000 G6PDH2r_revbin#407 <= 1000
 G6PDH2r_simple_thermo#288:      h_c_potential#201 - g6p_c_potential#210
                                 + _6pgl_c_potential#218 - nadp_c_potential#286
                                 + nadph_c_potential#287
                                 + 1000 G6PDH2r_revbin#407
                                 - RgG6PDH2r_simple_thermo#288  = 0
 GAPD_revbinF#289:               GAPD#144 - 1000 GAPD_revbin#409 <= 0
 GAPD_revbinR#290:               GAPD_reverse_459c1#145 + 1000 GAPD_revbin#409
                                 <= 1000
 GAPD_simple_thermo#291:         h_c_potential#201 + _13dpg_c_potential#213
                                 - nad_c_potential#223 + nadh_c_potential#224
                                 - pi_c_potential#235 - g3p_c_potential#283
                                 + 1000 GAPD_revbin#409
                                 - RgGAPD_simple_thermo#291  = 0
 GLCpts_revbinF#292:             GLCpts#146 - 1000 GLCpts_revbin#411 <= 0
 GLCpts_revbinR#293:             GLCpts_reverse_a52ae#147
                                 + 1000 GLCpts_revbin#411 <= 1000
 GLCpts_simple_thermo#294:       pyr_c_potential#207 + g6p_c_potential#210
                                 - pep_c_potential#251 - glc__D_e_potential#361
                                 + 1000 GLCpts_revbin#411
                                 - RgGLCpts_simple_thermo#294  = 0
 GLNS_revbinF#295:               GLNS#148 - 1000 GLNS_revbin#413 <= 0
 GLNS_revbinR#296:               GLNS_reverse_59581#149 + 1000 GLNS_revbin#413
                                 <= 1000
 GLNS_simple_thermo#297:         adp_c_potential#197 - atp_c_potential#198
                                 + h_c_potential#201 + pi_c_potential#235
                                 + gln__L_c_potential#284
                                 - glu__L_c_potential#285
                                 + 1000 GLNS_revbin#413 - nh4_c_potential#415
                                 - RgGLNS_simple_thermo#297  = 0
 GLNabc_revbinF#298:             GLNabc#150 - 1000 GLNabc_revbin#416 <= 0
 GLNabc_revbinR#299:             GLNabc_reverse_1d82a#151
                                 + 1000 GLNabc_revbin#416 <= 1000
 GLNabc_simple_thermo#300:       adp_c_potential#197 - atp_c_potential#198
                                 + h_c_potential#201 - h2o_c_potential#219
                                 + pi_c_potential#235 + gln__L_c_potential#284
                                 - gln__L_e_potential#364
                                 + 1000 GLNabc_revbin#416
                                 - RgGLNabc_simple_thermo#300  = 0
 GLUDy_revbinF#301:              GLUDy#152 - 1000 GLUDy_revbin#418 <= 0
 GLUDy_revbinR#302:              GLUDy_reverse_fa4e7#153
                                 + 1000 GLUDy_revbin#418 <= 1000
 GLUDy_simple_thermo#303:        h_c_potential#201 - h2o_c_potential#219
                                 + akg_c_potential#227 - glu__L_c_potential#285
                                 - nadp_c_potential#286 + nadph_c_potential#287
                                 + nh4_c_potential#415 + 1000 GLUDy_revbin#418
                                 - RgGLUDy_simple_thermo#303  = 0
 GLUN_revbinF#304:               GLUN#154 - 1000 GLUN_revbin#420 <= 0
 GLUN_revbinR#305:               GLUN_reverse_4ccdb#155 + 1000 GLUN_revbin#420
                                 <= 1000
 GLUN_simple_thermo#306:         - h2o_c_potential#219 - gln__L_c_potential#284
                                 + glu__L_c_potential#285 + nh4_c_potential#415
                                 + 1000 GLUN_revbin#420
                                 - RgGLUN_simple_thermo#306  = 0
 GLUSy_revbinF#307:              GLUSy#156 - 1000 GLUSy_revbin#422 <= 0
 GLUSy_revbinR#308:              GLUSy_reverse_6a00f#157
                                 + 1000 GLUSy_revbin#422 <= 1000
 GLUSy_simple_thermo#309:        - h_c_potential#201 - akg_c_potential#227
                                 - gln__L_c_potential#284
                                 + 2 glu__L_c_potential#285
                                 + nadp_c_potential#286 - nadph_c_potential#287
                                 + 1000 GLUSy_revbin#422
                                 - RgGLUSy_simple_thermo#309  = 0
 GLUt2r_revbinF#310:             GLUt2r#158 - 1000 GLUt2r_revbin#424 <= 0
 GLUt2r_revbinR#311:             GLUt2r_reverse_3e88a#159
                                 + 1000 GLUt2r_revbin#424 <= 1000
 GLUt2r_simple_thermo#312:       h_c_potential#201 - h_e_potential#229
                                 + glu__L_c_potential#285
                                 - glu__L_e_potential#367
                                 + 1000 GLUt2r_revbin#424
                                 - RgGLUt2r_simple_thermo#312  = 0
 GND_revbinF#313:                GND#160 - 1000 GND_revbin#426 <= 0
 GND_revbinR#314:                GND_reverse_eec5c#161 + 1000 GND_revbin#426
                                 <= 1000
 GND_simple_thermo#315:          - _6pgc_c_potential#217 + co2_c_potential#249
                                 - nadp_c_potential#286 + nadph_c_potential#287
                                 + ru5p__D_c_potential#297
                                 + 1000 GND_revbin#426
                                 - RgGND_simple_thermo#315  = 0
 H2Ot_revbinF#316:               H2Ot#162 - 1000 H2Ot_revbin#428 <= 0
 H2Ot_revbinR#317:               H2Ot_reverse_aa560#163 + 1000 H2Ot_revbin#428
                                 <= 1000
 H2Ot_simple_thermo#318:         h2o_c_potential#219 - h2o_e_potential#372
                                 + 1000 H2Ot_revbin#428
                                 - RgH2Ot_simple_thermo#318  = 0
 ICDHyr_revbinF#319:             ICDHyr#164 - 1000 ICDHyr_revbin#430 <= 0
 ICDHyr_revbinR#320:             ICDHyr_reverse_7f84b#165
                                 + 1000 ICDHyr_revbin#430 <= 1000
 ICDHyr_simple_thermo#321:       akg_c_potential#227 + co2_c_potential#249
                                 - icit_c_potential#258 - nadp_c_potential#286
                                 + nadph_c_potential#287
                                 + 1000 ICDHyr_revbin#430
                                 - RgICDHyr_simple_thermo#321  = 0
 ICL_revbinF#322:                ICL#166 - 1000 ICL_revbin#432 <= 0
 ICL_revbinR#323:                ICL_reverse_2f27e#167 + 1000 ICL_revbin#432
                                 <= 1000
 ICL_simple_thermo#324:          - icit_c_potential#258 + succ_c_potential#305
                                 + 1000 ICL_revbin#432 + glx_c_potential#434
                                 - RgICL_simple_thermo#324  = 0
 LDH_D_revbinF#325:              LDH_D#168 - 1000 LDH_D_revbin#435 <= 0
 LDH_D_revbinR#326:              LDH_D_reverse_f8507#169
                                 + 1000 LDH_D_revbin#435 <= 1000
 LDH_D_simple_thermo#327:        h_c_potential#201 + pyr_c_potential#207
                                 - nad_c_potential#223 + nadh_c_potential#224
                                 - lac__D_c_potential#314
                                 + 1000 LDH_D_revbin#435
                                 - RgLDH_D_simple_thermo#327  = 0
 MALS_revbinF#328:               MALS#170 - 1000 MALS_revbin#437 <= 0
 MALS_revbinR#329:               MALS_reverse_d7382#171 + 1000 MALS_revbin#437
                                 <= 1000
 MALS_simple_thermo#330:         h_c_potential#201 - accoa_c_potential#204
                                 + coa_c_potential#205 - h2o_c_potential#219
                                 + mal__L_c_potential#404 - glx_c_potential#434
                                 + 1000 MALS_revbin#437
                                 - RgMALS_simple_thermo#330  = 0
 MALt2_2_revbinF#331:            MALt2_2#172 - 1000 MALt2_2_revbin#439 <= 0
 MALt2_2_revbinR#332:            MALt2_2_reverse_a635f#173
                                 + 1000 MALt2_2_revbin#439 <= 1000
 MALt2_2_simple_thermo#333:      2 h_c_potential#201 - 2 h_e_potential#229
                                 - mal__L_e_potential#377
                                 + mal__L_c_potential#404
                                 + 1000 MALt2_2_revbin#439
                                 - RgMALt2_2_simple_thermo#333  = 0
 MDH_revbinF#334:                MDH#174 - 1000 MDH_revbin#441 <= 0
 MDH_revbinR#335:                MDH_reverse_ee52c#175 + 1000 MDH_revbin#441
                                 <= 1000
 MDH_simple_thermo#336:          h_c_potential#201 - nad_c_potential#223
                                 + nadh_c_potential#224 + oaa_c_potential#250
                                 - mal__L_c_potential#404 + 1000 MDH_revbin#441
                                 - RgMDH_simple_thermo#336  = 0
 ME1_revbinF#337:                ME1#176 - 1000 ME1_revbin#443 <= 0
 ME1_revbinR#338:                ME1_reverse_9736c#177 + 1000 ME1_revbin#443
                                 <= 1000
 ME1_simple_thermo#339:          pyr_c_potential#207 - nad_c_potential#223
                                 + nadh_c_potential#224 + co2_c_potential#249
                                 - mal__L_c_potential#404 + 1000 ME1_revbin#443
                                 - RgME1_simple_thermo#339  = 0
 ME2_revbinF#340:                ME2#178 - 1000 ME2_revbin#445 <= 0
 ME2_revbinR#341:                ME2_reverse_2b0a2#179 + 1000 ME2_revbin#445
                                 <= 1000
 ME2_simple_thermo#342:          pyr_c_potential#207 + co2_c_potential#249
                                 - nadp_c_potential#286 + nadph_c_potential#287
                                 - mal__L_c_potential#404 + 1000 ME2_revbin#445
                                 - RgME2_simple_thermo#342  = 0
 NADH16_revbinF#343:             NADH16#180 - 1000 NADH16_revbin#447 <= 0
 NADH16_revbinR#344:             NADH16_reverse_330a2#181
                                 + 1000 NADH16_revbin#447 <= 1000
 NADH16_simple_thermo#345:       - 4 h_c_potential#201 + nad_c_potential#223
                                 - nadh_c_potential#224 + 3 h_e_potential#229
                                 - q8_c_potential#310 + q8h2_c_potential#311
                                 + 1000 NADH16_revbin#447
                                 - RgNADH16_simple_thermo#345  = 0
 NADTRHD_revbinF#346:            NADTRHD#182 - 1000 NADTRHD_revbin#449 <= 0
 NADTRHD_revbinR#347:            NADTRHD_reverse_49725#183
                                 + 1000 NADTRHD_revbin#449 <= 1000
 NADTRHD_simple_thermo#348:      - nad_c_potential#223 + nadh_c_potential#224
                                 + nadp_c_potential#286 - nadph_c_potential#287
                                 + 1000 NADTRHD_revbin#449
                                 - RgNADTRHD_simple_thermo#348  = 0
 NH4t_revbinF#349:               NH4t#184 - 1000 NH4t_revbin#451 <= 0
 NH4t_revbinR#350:               NH4t_reverse_551ee#185 + 1000 NH4t_revbin#451
                                 <= 1000
 NH4t_simple_thermo#351:         - nh4_e_potential#380 + nh4_c_potential#415
                                 + 1000 NH4t_revbin#451
                                 - RgNH4t_simple_thermo#351  = 0
 O2t_revbinF#352:                O2t#186 - 1000 O2t_revbin#453 <= 0
 O2t_revbinR#353:                O2t_reverse_4d957#187 + 1000 O2t_revbin#453
                                 <= 1000
 O2t_simple_thermo#354:          o2_c_potential#309 - o2_e_potential#383
                                 + 1000 O2t_revbin#453
                                 - RgO2t_simple_thermo#354  = 0
 PDH_revbinF#355:                PDH#188 - 1000 PDH_revbin#455 <= 0
 PDH_revbinR#356:                PDH_reverse_ca160#189 + 1000 PDH_revbin#455
                                 <= 1000
 PDH_simple_thermo#357:          accoa_c_potential#204 - coa_c_potential#205
                                 - pyr_c_potential#207 - nad_c_potential#223
                                 + nadh_c_potential#224 + co2_c_potential#249
                                 + 1000 PDH_revbin#455
                                 - RgPDH_simple_thermo#357  = 0
 PFL_full_thermo#358:            1000 PFL_revbin#202
                                 + accoa_c_ln_concentration#457
                                 - coa_c_ln_concentration#458
                                 + for_c_ln_concentration#459
                                 - pyr_c_ln_concentration#460  = 1000
 PGI_full_thermo#359:            f6p_c_ln_concentration#192
                                 + 1000 PGI_revbin#208
                                 - g6p_c_ln_concentration#461  = 1000
 PGK_full_thermo#360:            adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 + 1000 PGK_revbin#211
                                 + _13dpg_c_ln_concentration#462
                                 - _3pg_c_ln_concentration#463  = 1000
 PGL_full_thermo#361:            - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PGL_revbin#215
                                 + _6pgc_c_ln_concentration#464
                                 - _6pgl_c_ln_concentration#465
                                 - h2o_c_ln_concentration#466  = 1000
 ACALD_full_thermo#362:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 ACALD_revbin#220
                                 + accoa_c_ln_concentration#457
                                 - coa_c_ln_concentration#458
                                 - acald_c_ln_concentration#467
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469  = 1000
 AKGt2r_full_thermo#363:         - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 AKGt2r_revbin#225
                                 + akg_c_ln_concentration#470
                                 - akg_e_ln_concentration#471
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                  = 1000
 PGM_full_thermo#364:            1000 PGM_revbin#230
                                 + _3pg_c_ln_concentration#463
                                 - _2pg_c_ln_concentration#473  = 1000
 PIt2r_full_thermo#365:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PIt2r_revbin#233
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 + pi_c_ln_concentration#474
                                 - pi_e_ln_concentration#475  = 1000
 ALCD2x_full_thermo#366:         - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 ALCD2x_revbin#237
                                 + acald_c_ln_concentration#467
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - _etoh_c_ln_concentration#476  = 1000
 ACALDt_full_thermo#367:         1000 ACALDt_revbin#240
                                 + acald_c_ln_concentration#467
                                 - acald_e_ln_concentration#477  = 1000
 ACKr_full_thermo#368:           adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 + 1000 ACKr_revbin#243
                                 - ac_c_ln_concentration#478
                                 + actp_c_ln_concentration#479  = 1000
 PPC_full_thermo#369:            - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PPC_revbin#247
                                 - h2o_c_ln_concentration#466
                                 + pi_c_ln_concentration#474
                                 + 0.105360515657826 co2_c_ln_concentration#480
                                 + oaa_c_ln_concentration#481
                                 - pep_c_ln_concentration#482  = 1000
 ACONTa_full_thermo#370:         1000 ACONTa_revbin#252
                                 + h2o_c_ln_concentration#466
                                 + acon_C_c_ln_concentration#483
                                 - cit_c_ln_concentration#484  = 1000
 ACONTb_full_thermo#371:         1000 ACONTb_revbin#256
                                 - h2o_c_ln_concentration#466
                                 - acon_C_c_ln_concentration#483
                                 + icit_c_ln_concentration#485  = 1000
 ATPM_full_thermo#372:           adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 ATPM_revbin#259
                                 - h2o_c_ln_concentration#466
                                 + pi_c_ln_concentration#474  = 1000
 PPCK_full_thermo#373:           adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 + 1000 PPCK_revbin#261
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 - oaa_c_ln_concentration#481
                                 + pep_c_ln_concentration#482  = 1000
 ACt2r_full_thermo#374:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 ACt2r_revbin#263
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 + ac_c_ln_concentration#478
                                 - ac_e_ln_concentration#486  = 1000
 PPS_full_thermo#375:            - atp_c_ln_concentration#191
                                 - 0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 PPS_revbin#266
                                 - pyr_c_ln_concentration#460
                                 - h2o_c_ln_concentration#466
                                 + pi_c_ln_concentration#474
                                 + pep_c_ln_concentration#482
                                 + amp_c_ln_concentration#487  = 1000
 ADK1_full_thermo#376:           2 adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 + 1000 ADK1_revbin#269
                                 - amp_c_ln_concentration#487  = 1000
 AKGDH_full_thermo#377:          1000 AKGDH_revbin#271
                                 - coa_c_ln_concentration#458
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - akg_c_ln_concentration#470
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 + succoa_c_ln_concentration#488  = 1000
 ATPS4r_full_thermo#378:         - adp_c_ln_concentration#190
                                 + atp_c_ln_concentration#191
                                 - 0.0606081219525584 h_c_ln_concentration#194
                                 + 1000 ATPS4r_revbin#274
                                 + h2o_c_ln_concentration#466
                                 + 0.0808108292700779 h_e_ln_concentration#472
                                 - pi_c_ln_concentration#474  = 1000
 PTAr_full_thermo#379:           1000 PTAr_revbin#276
                                 - accoa_c_ln_concentration#457
                                 + coa_c_ln_concentration#458
                                 - pi_c_ln_concentration#474
                                 + actp_c_ln_concentration#479  = 1000
 PYK_full_thermo#380:            - adp_c_ln_concentration#190
                                 + atp_c_ln_concentration#191
                                 + 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PYK_revbin#278
                                 + pyr_c_ln_concentration#460
                                 - pep_c_ln_concentration#482  = 1000
 BIOMASS_Ecoli_core_w_GAM_full_thermo#381: 59.81 adp_c_ln_concentration#190
                                 - 59.81 atp_c_ln_concentration#191
                                 - 0.0709 f6p_c_ln_concentration#192
                                 - 1.20832392466084 h_c_ln_concentration#194
                                 + 1000 BIOMASS_Ecoli_core_w_GAM_revbin#280
                                 - 3.7478 accoa_c_ln_concentration#457
                                 + 3.7478 coa_c_ln_concentration#458
                                 - 2.8328 pyr_c_ln_concentration#460
                                 - 0.205 g6p_c_ln_concentration#461
                                 - 1.496 _3pg_c_ln_concentration#463
                                 - 59.81 h2o_c_ln_concentration#466
                                 - 3.547 nad_c_ln_concentration#468
                                 + 3.547 nadh_c_ln_concentration#469
                                 + 4.1182 akg_c_ln_concentration#470
                                 + 59.81 pi_c_ln_concentration#474
                                 - 1.7867 oaa_c_ln_concentration#481
                                 - 0.5191 pep_c_ln_concentration#482
                                 - 0.361 _e4p_c_ln_concentration#489
                                 - 0.129 g3p_c_ln_concentration#490
                                 - 0.2557 gln__L_c_ln_concentration#491
                                 - 4.9414 glu__L_c_ln_concentration#492
                                 + 13.0279 nadp_c_ln_concentration#493
                                 - 13.0279 nadph_c_ln_concentration#494
                                 - 0.8977 r5p_c_ln_concentration#495  = 1000
 PYRt2_full_thermo#382:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 PYRt2_revbin#289
                                 + pyr_c_ln_concentration#460
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 - pyr_e_ln_concentration#496  = 1000
 CO2t_full_thermo#383:           1000 CO2t_revbin#292
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 + 0.105360515657826 co2_e_ln_concentration#497
                                  = 1000
 RPE_full_thermo#384:            1000 RPE_revbin#295
                                 - ru5p__D_c_ln_concentration#498
                                 + xu5p__D_c_ln_concentration#499  = 1000
 CS_full_thermo#385:             - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 CS_revbin#299
                                 - accoa_c_ln_concentration#457
                                 + coa_c_ln_concentration#458
                                 - h2o_c_ln_concentration#466
                                 - oaa_c_ln_concentration#481
                                 + cit_c_ln_concentration#484  = 1000
 RPI_full_thermo#386:            1000 RPI_revbin#301
                                 - r5p_c_ln_concentration#495
                                 + ru5p__D_c_ln_concentration#498  = 1000
 SUCCt2_2_full_thermo#387:       - 0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 SUCCt2_2_revbin#303
                                 + 0.0404054146350389 h_e_ln_concentration#472
                                 + succ_c_ln_concentration#500
                                 - succ_e_ln_concentration#501  = 1000
 CYTBD_full_thermo#388:          0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 CYTBD_revbin#307
                                 + h2o_c_ln_concentration#466
                                 - 0.0404054146350389 h_e_ln_concentration#472
                                 + 0.0256466471937753 o2_c_ln_concentration#502
                                 + q8_c_ln_concentration#503
                                 - q8h2_c_ln_concentration#504  = 1000
 D_LACt2_full_thermo#389:        - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 D_LACt2_revbin#312
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 + lac__D_c_ln_concentration#505
                                 - lac__D_e_ln_concentration#506  = 1000
 _ENO_full_thermo#390:           1000 _ENO_revbin#316
                                 + h2o_c_ln_concentration#466
                                 - _2pg_c_ln_concentration#473
                                 + pep_c_ln_concentration#482  = 1000
 SUCCt3_full_thermo#391:         - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 SUCCt3_revbin#318
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 - succ_c_ln_concentration#500
                                 + succ_e_ln_concentration#501  = 1000
 _ETOHt2r_full_thermo#392:       - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 _ETOHt2r_revbin#320
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 + _etoh_c_ln_concentration#476
                                 - _etoh_e_ln_concentration#507  = 1000
 SUCDi_full_thermo#393:          1000 SUCDi_revbin#323
                                 - succ_c_ln_concentration#500
                                 - q8_c_ln_concentration#503
                                 + q8h2_c_ln_concentration#504
                                 + fum_c_ln_concentration#508  = 1000
 SUCOAS_full_thermo#394:         adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 + 1000 SUCOAS_revbin#326
                                 - coa_c_ln_concentration#458
                                 + pi_c_ln_concentration#474
                                 + succoa_c_ln_concentration#488
                                 - succ_c_ln_concentration#500  = 1000
 TALA_full_thermo#395:           f6p_c_ln_concentration#192
                                 + 1000 TALA_revbin#328
                                 + _e4p_c_ln_concentration#489
                                 - g3p_c_ln_concentration#490
                                 - s7p_c_ln_concentration#509  = 1000
 THD2_full_thermo#396:           - 0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 THD2_revbin#331
                                 + nad_c_ln_concentration#468
                                 - nadh_c_ln_concentration#469
                                 + 0.0404054146350389 h_e_ln_concentration#472
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494  = 1000
 TKT1_full_thermo#397:           1000 TKT1_revbin#333
                                 + g3p_c_ln_concentration#490
                                 - r5p_c_ln_concentration#495
                                 - xu5p__D_c_ln_concentration#499
                                 + s7p_c_ln_concentration#509  = 1000
 TKT2_full_thermo#398:           f6p_c_ln_concentration#192
                                 + 1000 TKT2_revbin#335
                                 - _e4p_c_ln_concentration#489
                                 + g3p_c_ln_concentration#490
                                 - xu5p__D_c_ln_concentration#499  = 1000
 TPI_full_thermo#399:            1000 TPI_revbin#337
                                 + g3p_c_ln_concentration#490
                                 - dhap_c_ln_concentration#510  = 1000
 _EX_ac_e_full_thermo#400:       1000 _EX_ac_e_revbin#340
                                 - ac_e_ln_concentration#486  = 1000
 _EX_acald_e_full_thermo#401:    1000 _EX_acald_e_revbin#342
                                 - acald_e_ln_concentration#477  = 1000
 _EX_akg_e_full_thermo#402:      1000 _EX_akg_e_revbin#344
                                 - akg_e_ln_concentration#471  = 1000
 _EX_co2_e_full_thermo#403:      1000 _EX_co2_e_revbin#346
                                 + 0.105360515657826 co2_e_ln_concentration#497
                                  = 1000
 _EX_etoh_e_full_thermo#404:     1000 _EX_etoh_e_revbin#348
                                 - _etoh_e_ln_concentration#507  = 1000
 _EX_for_e_full_thermo#405:      1000 _EX_for_e_revbin#350
                                 - for_e_ln_concentration#511  = 1000
 _EX_fru_e_full_thermo#406:      1000 _EX_fru_e_revbin#353
                                 - fru_e_ln_concentration#512  = 1000
 _EX_fum_e_full_thermo#407:      1000 _EX_fum_e_revbin#356
                                 - fum_e_ln_concentration#513  = 1000
 _EX_glc__D_e_full_thermo#408:   1000 _EX_glc__D_e_revbin#359
                                 + 0.0618754037180875 glc__D_e_ln_concentration#514
                                  = 1000
 _EX_gln__L_e_full_thermo#409:   1000 _EX_gln__L_e_revbin#362
                                 - gln__L_e_ln_concentration#515  = 1000
 _EX_glu__L_e_full_thermo#410:   1000 _EX_glu__L_e_revbin#365
                                 - glu__L_e_ln_concentration#516  = 1000
 _EX_h_e_full_thermo#411:        1000 _EX_h_e_revbin#368
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                  = 1000
 _EX_h2o_e_full_thermo#412:      1000 _EX_h2o_e_revbin#370
                                 - h2o_e_ln_concentration#517  = 1000
 _EX_lac__D_e_full_thermo#413:   1000 _EX_lac__D_e_revbin#373
                                 - lac__D_e_ln_concentration#506  = 1000
 _EX_mal__L_e_full_thermo#414:   1000 _EX_mal__L_e_revbin#375
                                 - mal__L_e_ln_concentration#518  = 1000
 _EX_nh4_e_full_thermo#415:      1000 _EX_nh4_e_revbin#378
                                 - nh4_e_ln_concentration#519  = 1000
 _EX_o2_e_full_thermo#416:       1000 _EX_o2_e_revbin#381
                                 + 0.0512932943875506 o2_e_ln_concentration#520
                                  = 1000
 _EX_pi_e_full_thermo#417:       1000 _EX_pi_e_revbin#384
                                 - pi_e_ln_concentration#475  = 1000
 _EX_pyr_e_full_thermo#418:      1000 _EX_pyr_e_revbin#386
                                 - pyr_e_ln_concentration#496  = 1000
 _EX_succ_e_full_thermo#419:     1000 _EX_succ_e_revbin#388
                                 - succ_e_ln_concentration#501  = 1000
 FBA_full_thermo#420:            - fdp_c_ln_concentration#193
                                 + 1000 FBA_revbin#390
                                 + g3p_c_ln_concentration#490
                                 + dhap_c_ln_concentration#510  = 1000
 FBP_full_thermo#421:            f6p_c_ln_concentration#192
                                 - fdp_c_ln_concentration#193
                                 + 1000 FBP_revbin#392
                                 - h2o_c_ln_concentration#466
                                 + pi_c_ln_concentration#474  = 1000
 FORt2_full_thermo#422:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 FORt2_revbin#394
                                 + for_c_ln_concentration#459
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 - for_e_ln_concentration#511  = 1000
 FORt_full_thermo#423:           1000 FORt_revbin#396
                                 + for_c_ln_concentration#459
                                 - for_e_ln_concentration#511  = 1000
 FRD7_full_thermo#424:           1000 FRD7_revbin#398
                                 + succ_c_ln_concentration#500
                                 + q8_c_ln_concentration#503
                                 - q8h2_c_ln_concentration#504
                                 - fum_c_ln_concentration#508  = 1000
 FRUpts2_full_thermo#425:        f6p_c_ln_concentration#192
                                 + 1000 FRUpts2_revbin#400
                                 + pyr_c_ln_concentration#460
                                 - pep_c_ln_concentration#482
                                 - fru_e_ln_concentration#512  = 1000
 FUM_full_thermo#426:            1000 FUM_revbin#402
                                 - h2o_c_ln_concentration#466
                                 - fum_c_ln_concentration#508
                                 + mal__L_c_ln_concentration#521  = 1000
 FUMt2_2_full_thermo#427:        - 0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 FUMt2_2_revbin#405
                                 + 0.0404054146350389 h_e_ln_concentration#472
                                 + fum_c_ln_concentration#508
                                 - fum_e_ln_concentration#513  = 1000
 G6PDH2r_full_thermo#428:        - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 G6PDH2r_revbin#407
                                 - g6p_c_ln_concentration#461
                                 + _6pgl_c_ln_concentration#465
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494  = 1000
 GAPD_full_thermo#429:           - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GAPD_revbin#409
                                 + _13dpg_c_ln_concentration#462
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - pi_c_ln_concentration#474
                                 - g3p_c_ln_concentration#490  = 1000
 GLCpts_full_thermo#430:         1000 GLCpts_revbin#411
                                 + pyr_c_ln_concentration#460
                                 + g6p_c_ln_concentration#461
                                 - pep_c_ln_concentration#482
                                 + 0.0618754037180875 glc__D_e_ln_concentration#514
                                  = 1000
 GLNS_full_thermo#431:           adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GLNS_revbin#413
                                 + pi_c_ln_concentration#474
                                 + gln__L_c_ln_concentration#491
                                 - glu__L_c_ln_concentration#492
                                 - nh4_c_ln_concentration#522  = 1000
 GLNabc_full_thermo#432:         adp_c_ln_concentration#190
                                 - atp_c_ln_concentration#191
                                 - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GLNabc_revbin#416
                                 - h2o_c_ln_concentration#466
                                 + pi_c_ln_concentration#474
                                 + gln__L_c_ln_concentration#491
                                 - gln__L_e_ln_concentration#515  = 1000
 GLUDy_full_thermo#433:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GLUDy_revbin#418
                                 - h2o_c_ln_concentration#466
                                 + akg_c_ln_concentration#470
                                 - glu__L_c_ln_concentration#492
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494
                                 + nh4_c_ln_concentration#522  = 1000
 GLUN_full_thermo#434:           1000 GLUN_revbin#420
                                 - h2o_c_ln_concentration#466
                                 - gln__L_c_ln_concentration#491
                                 + glu__L_c_ln_concentration#492
                                 + nh4_c_ln_concentration#522  = 1000
 GLUSy_full_thermo#435:          0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GLUSy_revbin#422
                                 - akg_c_ln_concentration#470
                                 - gln__L_c_ln_concentration#491
                                 + 2 glu__L_c_ln_concentration#492
                                 + nadp_c_ln_concentration#493
                                 - nadph_c_ln_concentration#494  = 1000
 GLUt2r_full_thermo#436:         - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 GLUt2r_revbin#424
                                 + 0.0202027073175195 h_e_ln_concentration#472
                                 + glu__L_c_ln_concentration#492
                                 - glu__L_e_ln_concentration#516  = 1000
 GND_full_thermo#437:            1000 GND_revbin#426
                                 - _6pgc_c_ln_concentration#464
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494
                                 + ru5p__D_c_ln_concentration#498  = 1000
 H2Ot_full_thermo#438:           1000 H2Ot_revbin#428
                                 + h2o_c_ln_concentration#466
                                 - h2o_e_ln_concentration#517  = 1000
 ICDHyr_full_thermo#439:         1000 ICDHyr_revbin#430
                                 + akg_c_ln_concentration#470
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 - icit_c_ln_concentration#485
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494  = 1000
 ICL_full_thermo#440:            1000 ICL_revbin#432
                                 - icit_c_ln_concentration#485
                                 + succ_c_ln_concentration#500
                                 + glx_c_ln_concentration#523  = 1000
 LDH_D_full_thermo#441:          - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 LDH_D_revbin#435
                                 + pyr_c_ln_concentration#460
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - lac__D_c_ln_concentration#505  = 1000
 MALS_full_thermo#442:           - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 MALS_revbin#437
                                 - accoa_c_ln_concentration#457
                                 + coa_c_ln_concentration#458
                                 - h2o_c_ln_concentration#466
                                 + mal__L_c_ln_concentration#521
                                 - glx_c_ln_concentration#523  = 1000
 MALt2_2_full_thermo#443:        - 0.0404054146350389 h_c_ln_concentration#194
                                 + 1000 MALt2_2_revbin#439
                                 + 0.0404054146350389 h_e_ln_concentration#472
                                 - mal__L_e_ln_concentration#518
                                 + mal__L_c_ln_concentration#521  = 1000
 MDH_full_thermo#444:            - 0.0202027073175195 h_c_ln_concentration#194
                                 + 1000 MDH_revbin#441
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 + oaa_c_ln_concentration#481
                                 - mal__L_c_ln_concentration#521  = 1000
 ME1_full_thermo#445:            1000 ME1_revbin#443
                                 + pyr_c_ln_concentration#460
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 - mal__L_c_ln_concentration#521  = 1000
 ME2_full_thermo#446:            1000 ME2_revbin#445
                                 + pyr_c_ln_concentration#460
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                 - nadp_c_ln_concentration#493
                                 + nadph_c_ln_concentration#494
                                 - mal__L_c_ln_concentration#521  = 1000
 NADH16_full_thermo#447:         0.0808108292700779 h_c_ln_concentration#194
                                 + 1000 NADH16_revbin#447
                                 + nad_c_ln_concentration#468
                                 - nadh_c_ln_concentration#469
                                 - 0.0606081219525584 h_e_ln_concentration#472
                                 - q8_c_ln_concentration#503
                                 + q8h2_c_ln_concentration#504  = 1000
 NADTRHD_full_thermo#448:        1000 NADTRHD_revbin#449
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 + nadp_c_ln_concentration#493
                                 - nadph_c_ln_concentration#494  = 1000
 NH4t_full_thermo#449:           1000 NH4t_revbin#451
                                 - nh4_e_ln_concentration#519
                                 + nh4_c_ln_concentration#522  = 1000
 O2t_full_thermo#450:            1000 O2t_revbin#453
                                 - 0.0512932943875506 o2_c_ln_concentration#502
                                 + 0.0512932943875506 o2_e_ln_concentration#520
                                  = 1000
 PDH_full_thermo#451:            1000 PDH_revbin#455
                                 + accoa_c_ln_concentration#457
                                 - coa_c_ln_concentration#458
                                 - pyr_c_ln_concentration#460
                                 - nad_c_ln_concentration#468
                                 + nadh_c_ln_concentration#469
                                 - 0.105360515657826 co2_c_ln_concentration#480
                                  = 1000
Bounds
      PFK#0 = 1000
      PFK_reverse_d24a6#1 = 0
 0 <= PFL#2 <= 1000
      PFL_reverse_af9ec#3 = 0
 0 <= PGI#4 <= 1000
 0 <= PGI_reverse_27efc#5 <= 1000
 0 <= PGK#6 <= 1000
 0 <= PGK_reverse_02696#7 <= 1000
 0 <= PGL#8 <= 1000
      PGL_reverse_2bb6b#9 = 0
 0 <= ACALD#10 <= 1000
 0 <= ACALD_reverse_fda2b#11 <= 1000
 0 <= AKGt2r#12 <= 1000
 0 <= AKGt2r_reverse_5d500#13 <= 1000
 0 <= PGM#14 <= 1000
 0 <= PGM_reverse_fc9af#15 <= 1000
 0 <= PIt2r#16 <= 1000
 0 <= PIt2r_reverse_1cd61#17 <= 1000
 0 <= ALCD2x#18 <= 1000
 0 <= ALCD2x_reverse_5d107#19 <= 1000
 0 <= ACALDt#20 <= 1000
 0 <= ACALDt_reverse_858fa#21 <= 1000
 0 <= ACKr#22 <= 1000
 0 <= ACKr_reverse_b49c0#23 <= 1000
 0 <= PPC#24 <= 1000
      PPC_reverse_e854a#25 = 0
 0 <= ACONTa#26 <= 1000
 0 <= ACONTa_reverse_cad6d#27 <= 1000
 0 <= ACONTb#28 <= 1000
 0 <= ACONTb_reverse_e198a#29 <= 1000
 8.39 <= ATPM#30 <= 1000
      ATPM_reverse_5b752#31 = 0
 0 <= PPCK#32 <= 1000
      PPCK_reverse_2557d#33 = 0
 0 <= ACt2r#34 <= 1000
 0 <= ACt2r_reverse_64e05#35 <= 1000
 0 <= PPS#36 <= 1000
      PPS_reverse_1c319#37 = 0
 0 <= ADK1#38 <= 1000
 0 <= ADK1_reverse_a6f90#39 <= 1000
 0 <= AKGDH#40 <= 1000
      AKGDH_reverse_08bdc#41 = 0
 0 <= ATPS4r#42 <= 1000
 0 <= ATPS4r_reverse_64306#43 <= 1000
 0 <= PTAr#44 <= 1000
 0 <= PTAr_reverse_fce15#45 <= 1000
 0 <= PYK#46 <= 1000
      PYK_reverse_bc8ff#47 = 0
 0 <= BIOMASS_Ecoli_core_w_GAM#48 <= 1000
      BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49 = 0
 0 <= PYRt2#50 <= 1000
 0 <= PYRt2_reverse_a6e24#51 <= 1000
 0 <= CO2t#52 <= 1000
 0 <= CO2t_reverse_7c42f#53 <= 1000
 0 <= RPE#54 <= 1000
 0 <= RPE_reverse_a1b04#55 <= 1000
 0 <= CS#56 <= 1000
      CS_reverse_8d7e9#57 = 0
 0 <= RPI#58 <= 1000
 0 <= RPI_reverse_853a1#59 <= 1000
 0 <= SUCCt2_2#60 <= 1000
      SUCCt2_2_reverse_206ca#61 = 0
 0 <= CYTBD#62 <= 1000
      CYTBD_reverse_611ba#63 = 0
 0 <= D_LACt2#64 <= 1000
 0 <= D_LACt2_reverse_d0dc8#65 <= 1000
 0 <= _ENO#66 <= 1000
 0 <= _ENO_reverse_40eea#67 <= 1000
 0 <= SUCCt3#68 <= 1000
      SUCCt3_reverse_7d82a#69 = 0
 0 <= _ETOHt2r#70 <= 1000
 0 <= _ETOHt2r_reverse_b3d52#71 <= 1000
 0 <= SUCDi#72 <= 1000
      SUCDi_reverse_480f4#73 = 0
 0 <= SUCOAS#74 <= 1000
 0 <= SUCOAS_reverse_22958#75 <= 1000
 0 <= TALA#76 <= 1000
 0 <= TALA_reverse_adfda#77 <= 1000
 0 <= THD2#78 <= 1000
      THD2_reverse_f65dd#79 = 0
 0 <= TKT1#80 <= 1000
 0 <= TKT1_reverse_a1021#81 <= 1000
 0 <= TKT2#82 <= 1000
 0 <= TKT2_reverse_7ebc7#83 <= 1000
 0 <= TPI#84 <= 1000
 0 <= TPI_reverse_c2c3b#85 <= 1000
 0 <= _EX_ac_e#86 <= 1000
      _EX_ac_e_reverse_0be96#87 = 0
 0 <= _EX_acald_e#88 <= 1000
      _EX_acald_e_reverse_c096e#89 = 0
 0 <= _EX_akg_e#90 <= 1000
      _EX_akg_e_reverse_70d85#91 = 0
 0 <= _EX_co2_e#92 <= 1000
 0 <= _EX_co2_e_reverse_d0466#93 <= 1000
 0 <= _EX_etoh_e#94 <= 1000
      _EX_etoh_e_reverse_cc64f#95 = 0
 0 <= _EX_for_e#96 <= 1000
      _EX_for_e_reverse_23269#97 = 0
 0 <= _EX_fru_e#98 <= 1000
      _EX_fru_e_reverse_c3828#99 = 0
 0 <= _EX_fum_e#100 <= 1000
      _EX_fum_e_reverse_e3432#101 = 0
 0 <= _EX_glc__D_e#102 <= 1000
 0 <= _EX_glc__D_e_reverse_af641#103 <= 10
 0 <= _EX_gln__L_e#104 <= 1000
      _EX_gln__L_e_reverse_6a1a1#105 = 0
 0 <= _EX_glu__L_e#106 <= 1000
      _EX_glu__L_e_reverse_42f6c#107 = 0
 0 <= _EX_h_e#108 <= 1000
 0 <= _EX_h_e_reverse_3e0c5#109 <= 1000
 0 <= _EX_h2o_e#110 <= 1000
 0 <= _EX_h2o_e_reverse_3ced4#111 <= 1000
 0 <= _EX_lac__D_e#112 <= 1000
      _EX_lac__D_e_reverse_f95b4#113 = 0
 0 <= _EX_mal__L_e#114 <= 1000
      _EX_mal__L_e_reverse_af154#115 = 0
 0 <= _EX_nh4_e#116 <= 1000
 0 <= _EX_nh4_e_reverse_f9cc6#117 <= 1000
 0 <= _EX_o2_e#118 <= 1000
 0 <= _EX_o2_e_reverse_efa94#119 <= 1000
 0 <= _EX_pi_e#120 <= 1000
 0 <= _EX_pi_e_reverse_1fb09#121 <= 1000
 0 <= _EX_pyr_e#122 <= 1000
      _EX_pyr_e_reverse_1f6de#123 = 0
 0 <= _EX_succ_e#124 <= 1000
      _EX_succ_e_reverse_a9039#125 = 0
 0 <= FBA#126 <= 1000
 0 <= FBA_reverse_84806#127 <= 1000
 0 <= FBP#128 <= 1000
      FBP_reverse_bf2c9#129 = 0
 0 <= FORt2#130 <= 1000
      FORt2_reverse_89839#131 = 0
      FORt#132 = 0
 0 <= FORt_reverse_40f9f#133 <= 1000
 0 <= FRD7#134 <= 1000
      FRD7_reverse_e0172#135 = 0
 0 <= FRUpts2#136 <= 1000
      FRUpts2_reverse_58d3a#137 = 0
 0 <= FUM#138 <= 1000
 0 <= FUM_reverse_d3642#139 <= 1000
 0 <= FUMt2_2#140 <= 1000
      FUMt2_2_reverse_7b4b6#141 = 0
 0 <= G6PDH2r#142 <= 1000
 0 <= G6PDH2r_reverse_19ddf#143 <= 1000
 0 <= GAPD#144 <= 1000
 0 <= GAPD_reverse_459c1#145 <= 1000
 0 <= GLCpts#146 <= 1000
      GLCpts_reverse_a52ae#147 = 0
 0 <= GLNS#148 <= 1000
      GLNS_reverse_59581#149 = 0
 0 <= GLNabc#150 <= 1000
      GLNabc_reverse_1d82a#151 = 0
 0 <= GLUDy#152 <= 1000
 0 <= GLUDy_reverse_fa4e7#153 <= 1000
 0 <= GLUN#154 <= 1000
      GLUN_reverse_4ccdb#155 = 0
 0 <= GLUSy#156 <= 1000
      GLUSy_reverse_6a00f#157 = 0
 0 <= GLUt2r#158 <= 1000
 0 <= GLUt2r_reverse_3e88a#159 <= 1000
 0 <= GND#160 <= 1000
      GND_reverse_eec5c#161 = 0
 0 <= H2Ot#162 <= 1000
 0 <= H2Ot_reverse_aa560#163 <= 1000
 0 <= ICDHyr#164 <= 1000
 0 <= ICDHyr_reverse_7f84b#165 <= 1000
 0 <= ICL#166 <= 1000
      ICL_reverse_2f27e#167 = 0
 0 <= LDH_D#168 <= 1000
 0 <= LDH_D_reverse_f8507#169 <= 1000
 0 <= MALS#170 <= 1000
      MALS_reverse_d7382#171 = 0
 0 <= MALt2_2#172 <= 1000
      MALt2_2_reverse_a635f#173 = 0
 0 <= MDH#174 <= 1000
 0 <= MDH_reverse_ee52c#175 <= 1000
 0 <= ME1#176 <= 1000
      ME1_reverse_9736c#177 = 0
 0 <= ME2#178 <= 1000
      ME2_reverse_2b0a2#179 = 0
 0 <= NADH16#180 <= 1000
      NADH16_reverse_330a2#181 = 0
 0 <= NADTRHD#182 <= 1000
      NADTRHD_reverse_49725#183 = 0
 0 <= NH4t#184 <= 1000
 0 <= NH4t_reverse_551ee#185 <= 1000
 0 <= O2t#186 <= 1000
 0 <= O2t_reverse_4d957#187 <= 1000
 0 <= PDH#188 <= 1000
      PDH_reverse_ca160#189 = 0
 0.001 <= adp_c_ln_concentration#190 <= 20
 0.001 <= atp_c_ln_concentration#191 <= 20
 0.001 <= f6p_c_ln_concentration#192 <= 20
 0.001 <= fdp_c_ln_concentration#193 <= 20
      h_c_ln_concentration#194 = 0.0053
 0 <= PFK_revbin#195 <= 1
 0 <= PFK_potential#196 <= 1000
 0 <= adp_c_potential#197 <= 1000
 0 <= atp_c_potential#198 <= 1000
 0 <= f6p_c_potential#199 <= 1000
 0 <= fdp_c_potential#200 <= 1000
 0 <= h_c_potential#201 <= 1000
 0 <= PFL_revbin#202 <= 1
 0 <= PFL_potential#203 <= 1000
 0 <= accoa_c_potential#204 <= 1000
 0 <= coa_c_potential#205 <= 1000
 0 <= for_c_potential#206 <= 1000
 0 <= pyr_c_potential#207 <= 1000
 0 <= PGI_revbin#208 <= 1
 0 <= PGI_potential#209 <= 1000
 0 <= g6p_c_potential#210 <= 1000
 0 <= PGK_revbin#211 <= 1
 0 <= PGK_potential#212 <= 1000
 0 <= _13dpg_c_potential#213 <= 1000
 0 <= _3pg_c_potential#214 <= 1000
 0 <= PGL_revbin#215 <= 1
 0 <= PGL_potential#216 <= 1000
 0 <= _6pgc_c_potential#217 <= 1000
 0 <= _6pgl_c_potential#218 <= 1000
 0 <= h2o_c_potential#219 <= 1000
 0 <= ACALD_revbin#220 <= 1
 0 <= ACALD_potential#221 <= 1000
 0 <= acald_c_potential#222 <= 1000
 0 <= nad_c_potential#223 <= 1000
 0 <= nadh_c_potential#224 <= 1000
 0 <= AKGt2r_revbin#225 <= 1
 0 <= AKGt2r_potential#226 <= 1000
 0 <= akg_c_potential#227 <= 1000
 0 <= akg_e_potential#228 <= 1000
 0 <= h_e_potential#229 <= 1000
 0 <= PGM_revbin#230 <= 1
 0 <= PGM_potential#231 <= 1000
 0 <= _2pg_c_potential#232 <= 1000
 0 <= PIt2r_revbin#233 <= 1
 0 <= PIt2r_potential#234 <= 1000
 0 <= pi_c_potential#235 <= 1000
 0 <= pi_e_potential#236 <= 1000
 0 <= ALCD2x_revbin#237 <= 1
 0 <= ALCD2x_potential#238 <= 1000
 0 <= _etoh_c_potential#239 <= 1000
 0 <= ACALDt_revbin#240 <= 1
 0 <= ACALDt_potential#241 <= 1000
 0 <= acald_e_potential#242 <= 1000
 0 <= ACKr_revbin#243 <= 1
 0 <= ACKr_potential#244 <= 1000
 0 <= ac_c_potential#245 <= 1000
 0 <= actp_c_potential#246 <= 1000
 0 <= PPC_revbin#247 <= 1
 0 <= PPC_potential#248 <= 1000
 0 <= co2_c_potential#249 <= 1000
 0 <= oaa_c_potential#250 <= 1000
 0 <= pep_c_potential#251 <= 1000
 0 <= ACONTa_revbin#252 <= 1
 0 <= ACONTa_potential#253 <= 1000
 0 <= acon_C_c_potential#254 <= 1000
 0 <= cit_c_potential#255 <= 1000
 0 <= ACONTb_revbin#256 <= 1
 0 <= ACONTb_potential#257 <= 1000
 0 <= icit_c_potential#258 <= 1000
 0 <= ATPM_revbin#259 <= 1
 0 <= ATPM_potential#260 <= 1000
 0 <= PPCK_revbin#261 <= 1
 0 <= PPCK_potential#262 <= 1000
 0 <= ACt2r_revbin#263 <= 1
 0 <= ACt2r_potential#264 <= 1000
 0 <= ac_e_potential#265 <= 1000
 0 <= PPS_revbin#266 <= 1
 0 <= PPS_potential#267 <= 1000
 0 <= amp_c_potential#268 <= 1000
 0 <= ADK1_revbin#269 <= 1
 0 <= ADK1_potential#270 <= 1000
 0 <= AKGDH_revbin#271 <= 1
 0 <= AKGDH_potential#272 <= 1000
 0 <= succoa_c_potential#273 <= 1000
 0 <= ATPS4r_revbin#274 <= 1
 0 <= ATPS4r_potential#275 <= 1000
 0 <= PTAr_revbin#276 <= 1
 0 <= PTAr_potential#277 <= 1000
 0 <= PYK_revbin#278 <= 1
 0 <= PYK_potential#279 <= 1000
 0 <= BIOMASS_Ecoli_core_w_GAM_revbin#280 <= 1
 0 <= BIOMASS_Ecoli_core_w_GAM_potential#281 <= 1000
 0 <= _e4p_c_potential#282 <= 1000
 0 <= g3p_c_potential#283 <= 1000
 0 <= gln__L_c_potential#284 <= 1000
 0 <= glu__L_c_potential#285 <= 1000
 0 <= nadp_c_potential#286 <= 1000
 0 <= nadph_c_potential#287 <= 1000
 0 <= r5p_c_potential#288 <= 1000
 0 <= PYRt2_revbin#289 <= 1
 0 <= PYRt2_potential#290 <= 1000
 0 <= pyr_e_potential#291 <= 1000
 0 <= CO2t_revbin#292 <= 1
 0 <= CO2t_potential#293 <= 1000
 0 <= co2_e_potential#294 <= 1000
 0 <= RPE_revbin#295 <= 1
 0 <= RPE_potential#296 <= 1000
 0 <= ru5p__D_c_potential#297 <= 1000
 0 <= xu5p__D_c_potential#298 <= 1000
 0 <= CS_revbin#299 <= 1
 0 <= CS_potential#300 <= 1000
 0 <= RPI_revbin#301 <= 1
 0 <= RPI_potential#302 <= 1000
 0 <= SUCCt2_2_revbin#303 <= 1
 0 <= SUCCt2_2_potential#304 <= 1000
 0 <= succ_c_potential#305 <= 1000
 0 <= succ_e_potential#306 <= 1000
 0 <= CYTBD_revbin#307 <= 1
 0 <= CYTBD_potential#308 <= 1000
 0 <= o2_c_potential#309 <= 1000
 0 <= q8_c_potential#310 <= 1000
 0 <= q8h2_c_potential#311 <= 1000
 0 <= D_LACt2_revbin#312 <= 1
 0 <= D_LACt2_potential#313 <= 1000
 0 <= lac__D_c_potential#314 <= 1000
 0 <= lac__D_e_potential#315 <= 1000
 0 <= _ENO_revbin#316 <= 1
 0 <= _ENO_potential#317 <= 1000
 0 <= SUCCt3_revbin#318 <= 1
 0 <= SUCCt3_potential#319 <= 1000
 0 <= _ETOHt2r_revbin#320 <= 1
 0 <= _ETOHt2r_potential#321 <= 1000
 0 <= _etoh_e_potential#322 <= 1000
 0 <= SUCDi_revbin#323 <= 1
 0 <= SUCDi_potential#324 <= 1000
 0 <= fum_c_potential#325 <= 1000
 0 <= SUCOAS_revbin#326 <= 1
 0 <= SUCOAS_potential#327 <= 1000
 0 <= TALA_revbin#328 <= 1
 0 <= TALA_potential#329 <= 1000
 0 <= s7p_c_potential#330 <= 1000
 0 <= THD2_revbin#331 <= 1
 0 <= THD2_potential#332 <= 1000
 0 <= TKT1_revbin#333 <= 1
 0 <= TKT1_potential#334 <= 1000
 0 <= TKT2_revbin#335 <= 1
 0 <= TKT2_potential#336 <= 1000
 0 <= TPI_revbin#337 <= 1
 0 <= TPI_potential#338 <= 1000
 0 <= dhap_c_potential#339 <= 1000
 0 <= _EX_ac_e_revbin#340 <= 1
 0 <= _EX_ac_e_potential#341 <= 1000
 0 <= _EX_acald_e_revbin#342 <= 1
 0 <= _EX_acald_e_potential#343 <= 1000
 0 <= _EX_akg_e_revbin#344 <= 1
 0 <= _EX_akg_e_potential#345 <= 1000
 0 <= _EX_co2_e_revbin#346 <= 1
 0 <= _EX_co2_e_potential#347 <= 1000
 0 <= _EX_etoh_e_revbin#348 <= 1
 0 <= _EX_etoh_e_potential#349 <= 1000
 0 <= _EX_for_e_revbin#350 <= 1
 0 <= _EX_for_e_potential#351 <= 1000
 0 <= for_e_potential#352 <= 1000
 0 <= _EX_fru_e_revbin#353 <= 1
 0 <= _EX_fru_e_potential#354 <= 1000
 0 <= fru_e_potential#355 <= 1000
 0 <= _EX_fum_e_revbin#356 <= 1
 0 <= _EX_fum_e_potential#357 <= 1000
 0 <= fum_e_potential#358 <= 1000
 0 <= _EX_glc__D_e_revbin#359 <= 1
 0 <= _EX_glc__D_e_potential#360 <= 1000
 0 <= glc__D_e_potential#361 <= 1000
 0 <= _EX_gln__L_e_revbin#362 <= 1
 0 <= _EX_gln__L_e_potential#363 <= 1000
 0 <= gln__L_e_potential#364 <= 1000
 0 <= _EX_glu__L_e_revbin#365 <= 1
 0 <= _EX_glu__L_e_potential#366 <= 1000
 0 <= glu__L_e_potential#367 <= 1000
 0 <= _EX_h_e_revbin#368 <= 1
 0 <= _EX_h_e_potential#369 <= 1000
 0 <= _EX_h2o_e_revbin#370 <= 1
 0 <= _EX_h2o_e_potential#371 <= 1000
 0 <= h2o_e_potential#372 <= 1000
 0 <= _EX_lac__D_e_revbin#373 <= 1
 0 <= _EX_lac__D_e_potential#374 <= 1000
 0 <= _EX_mal__L_e_revbin#375 <= 1
 0 <= _EX_mal__L_e_potential#376 <= 1000
 0 <= mal__L_e_potential#377 <= 1000
 0 <= _EX_nh4_e_revbin#378 <= 1
 0 <= _EX_nh4_e_potential#379 <= 1000
 0 <= nh4_e_potential#380 <= 1000
 0 <= _EX_o2_e_revbin#381 <= 1
 0 <= _EX_o2_e_potential#382 <= 1000
 0 <= o2_e_potential#383 <= 1000
 0 <= _EX_pi_e_revbin#384 <= 1
 0 <= _EX_pi_e_potential#385 <= 1000
 0 <= _EX_pyr_e_revbin#386 <= 1
 0 <= _EX_pyr_e_potential#387 <= 1000
 0 <= _EX_succ_e_revbin#388 <= 1
 0 <= _EX_succ_e_potential#389 <= 1000
 0 <= FBA_revbin#390 <= 1
 0 <= FBA_potential#391 <= 1000
 0 <= FBP_revbin#392 <= 1
 0 <= FBP_potential#393 <= 1000
 0 <= FORt2_revbin#394 <= 1
 0 <= FORt2_potential#395 <= 1000
 0 <= FORt_revbin#396 <= 1
 0 <= FORt_potential#397 <= 1000
 0 <= FRD7_revbin#398 <= 1
 0 <= FRD7_potential#399 <= 1000
 0 <= FRUpts2_revbin#400 <= 1
 0 <= FRUpts2_potential#401 <= 1000
 0 <= FUM_revbin#402 <= 1
 0 <= FUM_potential#403 <= 1000
 0 <= mal__L_c_potential#404 <= 1000
 0 <= FUMt2_2_revbin#405 <= 1
 0 <= FUMt2_2_potential#406 <= 1000
 0 <= G6PDH2r_revbin#407 <= 1
 0 <= G6PDH2r_potential#408 <= 1000
 0 <= GAPD_revbin#409 <= 1
 0 <= GAPD_potential#410 <= 1000
 0 <= GLCpts_revbin#411 <= 1
 0 <= GLCpts_potential#412 <= 1000
 0 <= GLNS_revbin#413 <= 1
 0 <= GLNS_potential#414 <= 1000
 0 <= nh4_c_potential#415 <= 1000
 0 <= GLNabc_revbin#416 <= 1
 0 <= GLNabc_potential#417 <= 1000
 0 <= GLUDy_revbin#418 <= 1
 0 <= GLUDy_potential#419 <= 1000
 0 <= GLUN_revbin#420 <= 1
 0 <= GLUN_potential#421 <= 1000
 0 <= GLUSy_revbin#422 <= 1
 0 <= GLUSy_potential#423 <= 1000
 0 <= GLUt2r_revbin#424 <= 1
 0 <= GLUt2r_potential#425 <= 1000
 0 <= GND_revbin#426 <= 1
 0 <= GND_potential#427 <= 1000
 0 <= H2Ot_revbin#428 <= 1
 0 <= H2Ot_potential#429 <= 1000
 0 <= ICDHyr_revbin#430 <= 1
 0 <= ICDHyr_potential#431 <= 1000
 0 <= ICL_revbin#432 <= 1
 0 <= ICL_potential#433 <= 1000
 0 <= glx_c_potential#434 <= 1000
 0 <= LDH_D_revbin#435 <= 1
 0 <= LDH_D_potential#436 <= 1000
 0 <= MALS_revbin#437 <= 1
 0 <= MALS_potential#438 <= 1000
 0 <= MALt2_2_revbin#439 <= 1
 0 <= MALt2_2_potential#440 <= 1000
 0 <= MDH_revbin#441 <= 1
 0 <= MDH_potential#442 <= 1000
 0 <= ME1_revbin#443 <= 1
 0 <= ME1_potential#444 <= 1000
 0 <= ME2_revbin#445 <= 1
 0 <= ME2_potential#446 <= 1000
 0 <= NADH16_revbin#447 <= 1
 0 <= NADH16_potential#448 <= 1000
 0 <= NADTRHD_revbin#449 <= 1
 0 <= NADTRHD_potential#450 <= 1000
 0 <= NH4t_revbin#451 <= 1
 0 <= NH4t_potential#452 <= 1000
 0 <= O2t_revbin#453 <= 1
 0 <= O2t_potential#454 <= 1000
 0 <= PDH_revbin#455 <= 1
 0 <= PDH_potential#456 <= 1000
 0.001 <= accoa_c_ln_concentration#457 <= 20
 0.001 <= coa_c_ln_concentration#458 <= 20
 0.001 <= for_c_ln_concentration#459 <= 20
 0.001 <= pyr_c_ln_concentration#460 <= 20
 0.001 <= g6p_c_ln_concentration#461 <= 20
 0.001 <= _13dpg_c_ln_concentration#462 <= 20
 0.001 <= _3pg_c_ln_concentration#463 <= 20
 0.001 <= _6pgc_c_ln_concentration#464 <= 20
 0.001 <= _6pgl_c_ln_concentration#465 <= 20
 0.001 <= h2o_c_ln_concentration#466 <= 20
 0.001 <= acald_c_ln_concentration#467 <= 20
 0.001 <= nad_c_ln_concentration#468 <= 20
 0.001 <= nadh_c_ln_concentration#469 <= 20
 0.001 <= akg_c_ln_concentration#470 <= 20
 0.001 <= akg_e_ln_concentration#471 <= 20
      h_e_ln_concentration#472 = 0.0053
 0.001 <= _2pg_c_ln_concentration#473 <= 20
 0.001 <= pi_c_ln_concentration#474 <= 20
 0.001 <= pi_e_ln_concentration#475 <= 20
 0.001 <= _etoh_c_ln_concentration#476 <= 20
 0.001 <= acald_e_ln_concentration#477 <= 20
 0.001 <= ac_c_ln_concentration#478 <= 20
 0.001 <= actp_c_ln_concentration#479 <= 20
 20 <= co2_c_ln_concentration#480 <= 24
 0.001 <= oaa_c_ln_concentration#481 <= 20
 0.001 <= pep_c_ln_concentration#482 <= 20
 0.001 <= acon_C_c_ln_concentration#483 <= 20
 0.001 <= cit_c_ln_concentration#484 <= 20
 0.001 <= icit_c_ln_concentration#485 <= 20
 0.001 <= ac_e_ln_concentration#486 <= 20
 0.001 <= amp_c_ln_concentration#487 <= 20
 0.001 <= succoa_c_ln_concentration#488 <= 20
 0.001 <= _e4p_c_ln_concentration#489 <= 20
 0.001 <= g3p_c_ln_concentration#490 <= 20
 0.001 <= gln__L_c_ln_concentration#491 <= 20
 0.001 <= glu__L_c_ln_concentration#492 <= 20
 0.001 <= nadp_c_ln_concentration#493 <= 20
 0.001 <= nadph_c_ln_concentration#494 <= 20
 0.001 <= r5p_c_ln_concentration#495 <= 20
 0.001 <= pyr_e_ln_concentration#496 <= 20
 20 <= co2_e_ln_concentration#497 <= 24
 0.001 <= ru5p__D_c_ln_concentration#498 <= 20
 0.001 <= xu5p__D_c_ln_concentration#499 <= 20
 0.001 <= succ_c_ln_concentration#500 <= 20
 0.001 <= succ_e_ln_concentration#501 <= 20
      o2_c_ln_concentration#502 = 0.672
 0.001 <= q8_c_ln_concentration#503 <= 20
 0.001 <= q8h2_c_ln_concentration#504 <= 20
 0.001 <= lac__D_c_ln_concentration#505 <= 20
 0.001 <= lac__D_e_ln_concentration#506 <= 20
 0.001 <= _etoh_e_ln_concentration#507 <= 20
 0.001 <= fum_c_ln_concentration#508 <= 20
 0.001 <= s7p_c_ln_concentration#509 <= 20
 0.001 <= dhap_c_ln_concentration#510 <= 20
 0.001 <= for_e_ln_concentration#511 <= 20
 0.001 <= fru_e_ln_concentration#512 <= 20
 0.001 <= fum_e_ln_concentration#513 <= 20
      glc__D_e_ln_concentration#514 = 10
 0.001 <= gln__L_e_ln_concentration#515 <= 20
 0.001 <= glu__L_e_ln_concentration#516 <= 20
 0.001 <= h2o_e_ln_concentration#517 <= 20
 0.001 <= mal__L_e_ln_concentration#518 <= 20
 0.001 <= nh4_e_ln_concentration#519 <= 20
      o2_e_ln_concentration#520 = 0.672
 0.001 <= mal__L_c_ln_concentration#521 <= 20
 0.001 <= nh4_c_ln_concentration#522 <= 20
 0.001 <= glx_c_ln_concentration#523 <= 20
 0 <= RgPFK_simple_thermo#75 <= 1000
 0 <= RgPFL_simple_thermo#78 <= 1000
 0 <= RgPGI_simple_thermo#81 <= 1000
 0 <= RgPGK_simple_thermo#84 <= 1000
 0 <= RgPGL_simple_thermo#87 <= 1000
 0 <= RgACALD_simple_thermo#90 <= 1000
 0 <= RgAKGt2r_simple_thermo#93 <= 1000
 0 <= RgPGM_simple_thermo#96 <= 1000
 0 <= RgPIt2r_simple_thermo#99 <= 1000
 0 <= RgALCD2x_simple_thermo#102 <= 1000
 0 <= RgACALDt_simple_thermo#105 <= 1000
 0 <= RgACKr_simple_thermo#108 <= 1000
 0 <= RgPPC_simple_thermo#111 <= 1000
 0 <= RgACONTa_simple_thermo#114 <= 1000
 0 <= RgACONTb_simple_thermo#117 <= 1000
 0 <= RgATPM_simple_thermo#120 <= 1000
 0 <= RgPPCK_simple_thermo#123 <= 1000
 0 <= RgACt2r_simple_thermo#126 <= 1000
 0 <= RgPPS_simple_thermo#129 <= 1000
 0 <= RgADK1_simple_thermo#132 <= 1000
 0 <= RgAKGDH_simple_thermo#135 <= 1000
 0 <= RgATPS4r_simple_thermo#138 <= 1000
 0 <= RgPTAr_simple_thermo#141 <= 1000
 0 <= RgPYK_simple_thermo#144 <= 1000
 0 <= RgBIOMASS_Ecoli_core_w_GAM_simple_thermo#147 <= 1000
 0 <= RgPYRt2_simple_thermo#150 <= 1000
 0 <= RgCO2t_simple_thermo#153 <= 1000
 0 <= RgRPE_simple_thermo#156 <= 1000
 0 <= RgCS_simple_thermo#159 <= 1000
 0 <= RgRPI_simple_thermo#162 <= 1000
 0 <= RgSUCCt2_2_simple_thermo#165 <= 1000
 0 <= RgCYTBD_simple_thermo#168 <= 1000
 0 <= RgD_LACt2_simple_thermo#171 <= 1000
 0 <= Rg_ENO_simple_thermo#174 <= 1000
 0 <= RgSUCCt3_simple_thermo#177 <= 1000
 0 <= Rg_ETOHt2r_simple_thermo#180 <= 1000
 0 <= RgSUCDi_simple_thermo#183 <= 1000
 0 <= RgSUCOAS_simple_thermo#186 <= 1000
 0 <= RgTALA_simple_thermo#189 <= 1000
 0 <= RgTHD2_simple_thermo#192 <= 1000
 0 <= RgTKT1_simple_thermo#195 <= 1000
 0 <= RgTKT2_simple_thermo#198 <= 1000
 0 <= RgTPI_simple_thermo#201 <= 1000
 0 <= Rg_EX_ac_e_simple_thermo#204 <= 1000
 0 <= Rg_EX_acald_e_simple_thermo#207 <= 1000
 0 <= Rg_EX_akg_e_simple_thermo#210 <= 1000
 0 <= Rg_EX_co2_e_simple_thermo#213 <= 1000
 0 <= Rg_EX_etoh_e_simple_thermo#216 <= 1000
 0 <= Rg_EX_for_e_simple_thermo#219 <= 1000
 0 <= Rg_EX_fru_e_simple_thermo#222 <= 1000
 0 <= Rg_EX_fum_e_simple_thermo#225 <= 1000
 0 <= Rg_EX_glc__D_e_simple_thermo#228 <= 1000
 0 <= Rg_EX_gln__L_e_simple_thermo#231 <= 1000
 0 <= Rg_EX_glu__L_e_simple_thermo#234 <= 1000
 0 <= Rg_EX_h_e_simple_thermo#237 <= 1000
 0 <= Rg_EX_h2o_e_simple_thermo#240 <= 1000
 0 <= Rg_EX_lac__D_e_simple_thermo#243 <= 1000
 0 <= Rg_EX_mal__L_e_simple_thermo#246 <= 1000
 0 <= Rg_EX_nh4_e_simple_thermo#249 <= 1000
 0 <= Rg_EX_o2_e_simple_thermo#252 <= 1000
 0 <= Rg_EX_pi_e_simple_thermo#255 <= 1000
 0 <= Rg_EX_pyr_e_simple_thermo#258 <= 1000
 0 <= Rg_EX_succ_e_simple_thermo#261 <= 1000
 0 <= RgFBA_simple_thermo#264 <= 1000
 0 <= RgFBP_simple_thermo#267 <= 1000
 0 <= RgFORt2_simple_thermo#270 <= 1000
 0 <= RgFORt_simple_thermo#273 <= 1000
 0 <= RgFRD7_simple_thermo#276 <= 1000
 0 <= RgFRUpts2_simple_thermo#279 <= 1000
 0 <= RgFUM_simple_thermo#282 <= 1000
 0 <= RgFUMt2_2_simple_thermo#285 <= 1000
 0 <= RgG6PDH2r_simple_thermo#288 <= 1000
 0 <= RgGAPD_simple_thermo#291 <= 1000
 0 <= RgGLCpts_simple_thermo#294 <= 1000
 0 <= RgGLNS_simple_thermo#297 <= 1000
 0 <= RgGLNabc_simple_thermo#300 <= 1000
 0 <= RgGLUDy_simple_thermo#303 <= 1000
 0 <= RgGLUN_simple_thermo#306 <= 1000
 0 <= RgGLUSy_simple_thermo#309 <= 1000
 0 <= RgGLUt2r_simple_thermo#312 <= 1000
 0 <= RgGND_simple_thermo#315 <= 1000
 0 <= RgH2Ot_simple_thermo#318 <= 1000
 0 <= RgICDHyr_simple_thermo#321 <= 1000
 0 <= RgICL_simple_thermo#324 <= 1000
 0 <= RgLDH_D_simple_thermo#327 <= 1000
 0 <= RgMALS_simple_thermo#330 <= 1000
 0 <= RgMALt2_2_simple_thermo#333 <= 1000
 0 <= RgMDH_simple_thermo#336 <= 1000
 0 <= RgME1_simple_thermo#339 <= 1000
 0 <= RgME2_simple_thermo#342 <= 1000
 0 <= RgNADH16_simple_thermo#345 <= 1000
 0 <= RgNADTRHD_simple_thermo#348 <= 1000
 0 <= RgNH4t_simple_thermo#351 <= 1000
 0 <= RgO2t_simple_thermo#354 <= 1000
 0 <= RgPDH_simple_thermo#357 <= 1000
Binaries
 PFK_revbin#195  PFL_revbin#202  PGI_revbin#208  PGK_revbin#211 
 PGL_revbin#215  ACALD_revbin#220  AKGt2r_revbin#225  PGM_revbin#230 
 PIt2r_revbin#233  ALCD2x_revbin#237  ACALDt_revbin#240  ACKr_revbin#243 
 PPC_revbin#247  ACONTa_revbin#252  ACONTb_revbin#256  ATPM_revbin#259 
 PPCK_revbin#261  ACt2r_revbin#263  PPS_revbin#266  ADK1_revbin#269 
 AKGDH_revbin#271  ATPS4r_revbin#274  PTAr_revbin#276  PYK_revbin#278 
 BIOMASS_Ecoli_core_w_GAM_revbin#280  PYRt2_revbin#289  CO2t_revbin#292 
 RPE_revbin#295  CS_revbin#299  RPI_revbin#301  SUCCt2_2_revbin#303 
 CYTBD_revbin#307  D_LACt2_revbin#312  _ENO_revbin#316  SUCCt3_revbin#318 
 _ETOHt2r_revbin#320  SUCDi_revbin#323  SUCOAS_revbin#326  TALA_revbin#328 
 THD2_revbin#331  TKT1_revbin#333  TKT2_revbin#335  TPI_revbin#337 
 _EX_ac_e_revbin#340  _EX_acald_e_revbin#342  _EX_akg_e_revbin#344 
 _EX_co2_e_revbin#346  _EX_etoh_e_revbin#348  _EX_for_e_revbin#350 
 _EX_fru_e_revbin#353  _EX_fum_e_revbin#356  _EX_glc__D_e_revbin#359 
 _EX_gln__L_e_revbin#362  _EX_glu__L_e_revbin#365  _EX_h_e_revbin#368 
 _EX_h2o_e_revbin#370  _EX_lac__D_e_revbin#373  _EX_mal__L_e_revbin#375 
 _EX_nh4_e_revbin#378  _EX_o2_e_revbin#381  _EX_pi_e_revbin#384 
 _EX_pyr_e_revbin#386  _EX_succ_e_revbin#388  FBA_revbin#390  FBP_revbin#392 
 FORt2_revbin#394  FORt_revbin#396  FRD7_revbin#398  FRUpts2_revbin#400 
 FUM_revbin#402  FUMt2_2_revbin#405  G6PDH2r_revbin#407  GAPD_revbin#409 
 GLCpts_revbin#411  GLNS_revbin#413  GLNabc_revbin#416  GLUDy_revbin#418 
 GLUN_revbin#420  GLUSy_revbin#422  GLUt2r_revbin#424  GND_revbin#426 
 H2Ot_revbin#428  ICDHyr_revbin#430  ICL_revbin#432  LDH_D_revbin#435 
 MALS_revbin#437  MALt2_2_revbin#439  MDH_revbin#441  ME1_revbin#443 
 ME2_revbin#445  NADH16_revbin#447  NADTRHD_revbin#449  NH4t_revbin#451 
 O2t_revbin#453  PDH_revbin#455 
End
