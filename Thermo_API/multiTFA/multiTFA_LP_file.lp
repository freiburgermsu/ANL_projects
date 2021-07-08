\ENCODING=ISO-8859-1
\Problem name: C:\Users\ANDREW~1\AppData\Local\Temp\tmp63argws2.sav

Maximize
 d03bf8f7mdf9fm11ebm80d0m50eb71d9aedf#0: BIOMASS_Ecoli_core_w_GAM#48
                                 - BIOMASS_Ecoli_core_w_GAM_reverse_712e5#49
                                 + 0 lnc_h2o_e#197 + 0 lnc_h_c#198
                                 + 0 lnc_h_e#199 + 0 lnc_o2_e#214
                                 + 0 dG_err_h2o_e#269 + 0 dG_err_h_c#270
                                 + 0 dG_err_h_e#271 + 0 dG_err_o2_e#286
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
 directionality_PFK#72:          PFK#0 - 1000 indicator_PFK#336 <= 0
 directionality_PFK_reverse_d24a6#73: PFK_reverse_d24a6#1
                                 - 1000 indicator_PFK_reverse_d24a6#337 <= 0
 ind_PFK#74:                     dG_PFK#334 + 100000000 indicator_PFK#336
                                 <= 100000000
 ind_PFK_reverse_d24a6#75:       dG_PFK_reverse_d24a6#335
                                 + 100000000 indicator_PFK_reverse_d24a6#337
                                 <= 100000000
 delG_PFK#76:                    - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244
                                 + 2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_fdp_c#254 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 + dG_err_f6p_c#325
                                 - dG_err_fdp_c#326 + dG_PFK#334
                                  = -19.6607111367925
 delG_PFK_reverse_d24a6#77:      2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244
                                 - 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_fdp_c#254 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 - dG_err_f6p_c#325
                                 + dG_err_fdp_c#326 + dG_PFK_reverse_d24a6#335
                                  = 19.6607111367925
 directionality_PFL#78:          PFL#2 - 1000 indicator_PFL#340 <= 0
 directionality_PFL_reverse_af9ec#79: PFL_reverse_af9ec#3
                                 - 1000 indicator_PFL_reverse_af9ec#341 <= 0
 ind_PFL#80:                     dG_PFL#338 + 100000000 indicator_PFL#340
                                 <= 100000000
 ind_PFL_reverse_af9ec#81:       dG_PFL_reverse_af9ec#339
                                 + 100000000 indicator_PFL_reverse_af9ec#341
                                 <= 100000000
 delG_PFL#82:                    2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_coa_c#248
                                 - 2.4776265 lnc_for_c#255 + dG_err_pyr_c#295
                                 - dG_err_accoa_c#305 + dG_err_coa_c#320
                                 - dG_err_for_c#327 + dG_PFL#338
                                  = -19.5617610410578
 delG_PFL_reverse_af9ec#83:      - 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_coa_c#248
                                 + 2.4776265 lnc_for_c#255 - dG_err_pyr_c#295
                                 + dG_err_accoa_c#305 - dG_err_coa_c#320
                                 + dG_err_for_c#327 + dG_PFL_reverse_af9ec#339
                                  = 19.5617610410578
 directionality_PGI#84:          PGI#4 - 1000 indicator_PGI#344 <= 0
 directionality_PGI_reverse_27efc#85: PGI_reverse_27efc#5
                                 - 1000 indicator_PGI_reverse_27efc#345 <= 0
 ind_PGI#86:                     dG_PGI#342 + 100000000 indicator_PGI#344
                                 <= 100000000
 ind_PGI_reverse_27efc#87:       dG_PGI_reverse_27efc#343
                                 + 100000000 indicator_PGI_reverse_27efc#345
                                 <= 100000000
 delG_PGI#88:                    - 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_g6p_c#261 - dG_err_f6p_c#325
                                 + dG_err_g6p_c#333 + dG_PGI#342
                                  = 2.92274522278399
 delG_PGI_reverse_27efc#89:      2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_g6p_c#261 + dG_err_f6p_c#325
                                 - dG_err_g6p_c#333 + dG_PGI_reverse_27efc#343
                                  = -2.92274522278399
 directionality_PGK#90:          PGK#6 - 1000 indicator_PGK#348 <= 0
 directionality_PGK_reverse_02696#91: PGK_reverse_02696#7
                                 - 1000 indicator_PGK_reverse_02696#349 <= 0
 ind_PGK#92:                     dG_PGK#346 + 100000000 indicator_PGK#348
                                 <= 100000000
 ind_PGK_reverse_02696#93:       dG_PGK_reverse_02696#347
                                 + 100000000 indicator_PGK_reverse_02696#349
                                 <= 100000000
 delG_PGK#94:                    - 2.4776265 lnc_13dpg_c#210
                                 + 2.4776265 lnc_3pg_c#215
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244 - dG_err_13dpg_c#282
                                 + dG_err_3pg_c#287 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 + dG_PGK#346
                                  = 17.902307221889
 delG_PGK_reverse_02696#95:      2.4776265 lnc_13dpg_c#210
                                 - 2.4776265 lnc_3pg_c#215
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244 + dG_err_13dpg_c#282
                                 - dG_err_3pg_c#287 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_PGK_reverse_02696#347
                                  = -17.902307221889
 directionality_PGL#96:          PGL#8 - 1000 indicator_PGL#352 <= 0
 directionality_PGL_reverse_2bb6b#97: PGL_reverse_2bb6b#9
                                 - 1000 indicator_PGL_reverse_2bb6b#353 <= 0
 ind_PGL#98:                     dG_PGL#350 + 100000000 indicator_PGL#352
                                 <= 100000000
 ind_PGL_reverse_2bb6b#99:       dG_PGL_reverse_2bb6b#351
                                 + 100000000 indicator_PGL_reverse_2bb6b#353
                                 <= 100000000
 delG_PGL#100:                   2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_6pgc_c#218
                                 + 2.4776265 lnc_6pgl_c#220 + dG_err_h2o_c#268
                                 - dG_err_6pgc_c#290 + dG_err_6pgl_c#292
                                 + dG_PGL#350  = -25.297541511457
 delG_PGL_reverse_2bb6b#101:     - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_6pgc_c#218
                                 - 2.4776265 lnc_6pgl_c#220 - dG_err_h2o_c#268
                                 + dG_err_6pgc_c#290 - dG_err_6pgl_c#292
                                 + dG_PGL_reverse_2bb6b#351  = 25.297541511457
 directionality_ACALD#102:       ACALD#10 - 1000 indicator_ACALD#356 <= 0
 directionality_ACALD_reverse_fda2b#103: ACALD_reverse_fda2b#11
                                 - 1000 indicator_ACALD_reverse_fda2b#357 <= 0
 ind_ACALD#104:                  dG_ACALD#354 + 100000000 indicator_ACALD#356
                                 <= 100000000
 ind_ACALD_reverse_fda2b#105:    dG_ACALD_reverse_fda2b#355
                                 + 100000000 indicator_ACALD_reverse_fda2b#357
                                 <= 100000000
 delG_ACALD#106:                 2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_acald_c#230
                                 - 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_coa_c#248 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 + dG_err_acald_c#302
                                 - dG_err_accoa_c#305 + dG_err_coa_c#320
                                 + dG_ACALD#354  = -20.3309049468082
 delG_ACALD_reverse_fda2b#107:   - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_acald_c#230
                                 + 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_coa_c#248 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 - dG_err_acald_c#302
                                 + dG_err_accoa_c#305 - dG_err_coa_c#320
                                 + dG_ACALD_reverse_fda2b#355
                                  = 20.3309049468082
 directionality_AKGt2r#108:      AKGt2r#12 - 1000 indicator_AKGt2r#360 <= 0
 directionality_AKGt2r_reverse_5d500#109: AKGt2r_reverse_5d500#13
                                 - 1000 indicator_AKGt2r_reverse_5d500#361 <= 0
 ind_AKGt2r#110:                 dG_AKGt2r#358 + 100000000 indicator_AKGt2r#360
                                 <= 100000000
 ind_AKGt2r_reverse_5d500#111:   dG_AKGt2r_reverse_5d500#359
                                 + 100000000 indicator_AKGt2r_reverse_5d500#361
                                 <= 100000000
 delG_AKGt2r#112:                - 2.4776265 lnc_akg_c#241
                                 + 2.4776265 lnc_akg_e#242 - dG_err_akg_c#313
                                 + dG_err_akg_e#314 + dG_AKGt2r#358
                                  = 5.6551079482812
 delG_AKGt2r_reverse_5d500#113:  2.4776265 lnc_akg_c#241
                                 - 2.4776265 lnc_akg_e#242 + dG_err_akg_c#313
                                 - dG_err_akg_e#314
                                 + dG_AKGt2r_reverse_5d500#359
                                  = -5.6551079482812
 directionality_PGM#114:         PGM#14 - 1000 indicator_PGM#364 <= 0
 directionality_PGM_reverse_fc9af#115: PGM_reverse_fc9af#15
                                 - 1000 indicator_PGM_reverse_fc9af#365 <= 0
 ind_PGM#116:                    dG_PGM#362 + 100000000 indicator_PGM#364
                                 <= 100000000
 ind_PGM_reverse_fc9af#117:      dG_PGM_reverse_fc9af#363
                                 + 100000000 indicator_PGM_reverse_fc9af#365
                                 <= 100000000
 delG_PGM#118:                   2.4776265 lnc_2pg_c#213
                                 - 2.4776265 lnc_3pg_c#215 + dG_err_2pg_c#285
                                 - dG_err_3pg_c#287 + dG_PGM#362
                                  = -4.69474453090197
 delG_PGM_reverse_fc9af#119:     - 2.4776265 lnc_2pg_c#213
                                 + 2.4776265 lnc_3pg_c#215 - dG_err_2pg_c#285
                                 + dG_err_3pg_c#287 + dG_PGM_reverse_fc9af#363
                                  = 4.69474453090197
 directionality_PIt2r#120:       PIt2r#16 - 1000 indicator_PIt2r#368 <= 0
 directionality_PIt2r_reverse_1cd61#121: PIt2r_reverse_1cd61#17
                                 - 1000 indicator_PIt2r_reverse_1cd61#369 <= 0
 ind_PIt2r#122:                  dG_PIt2r#366 + 100000000 indicator_PIt2r#368
                                 <= 100000000
 ind_PIt2r_reverse_1cd61#123:    dG_PIt2r_reverse_1cd61#367
                                 + 100000000 indicator_PIt2r_reverse_1cd61#369
                                 <= 100000000
 delG_PIt2r#124:                 - 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_pi_e#221 - dG_err_pi_c#291
                                 + dG_err_pi_e#293 + dG_PIt2r#366
                                  = -0.482485756200172
 delG_PIt2r_reverse_1cd61#125:   2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_pi_e#221 + dG_err_pi_c#291
                                 - dG_err_pi_e#293 + dG_PIt2r_reverse_1cd61#367
                                  = 0.482485756200172
 directionality_ALCD2x#126:      ALCD2x#18 - 1000 indicator_ALCD2x#372 <= 0
 directionality_ALCD2x_reverse_5d107#127: ALCD2x_reverse_5d107#19
                                 - 1000 indicator_ALCD2x_reverse_5d107#373 <= 0
 ind_ALCD2x#128:                 dG_ALCD2x#370 + 100000000 indicator_ALCD2x#372
                                 <= 100000000
 ind_ALCD2x_reverse_5d107#129:   dG_ALCD2x_reverse_5d107#371
                                 + 100000000 indicator_ALCD2x_reverse_5d107#373
                                 <= 100000000
 delG_ALCD2x#130:                2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_acald_c#230
                                 + 2.4776265 lnc_etoh_c#251 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_acald_c#302
                                 + dG_err_etoh_c#323 + dG_ALCD2x#370
                                  = 17.9788594927938
 delG_ALCD2x_reverse_5d107#131:  - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_acald_c#230
                                 - 2.4776265 lnc_etoh_c#251 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_acald_c#302
                                 - dG_err_etoh_c#323
                                 + dG_ALCD2x_reverse_5d107#371
                                  = -17.9788594927938
 directionality_ACALDt#132:      ACALDt#20 - 1000 indicator_ACALDt#376 <= 0
 directionality_ACALDt_reverse_858fa#133: ACALDt_reverse_858fa#21
                                 - 1000 indicator_ACALDt_reverse_858fa#377 <= 0
 ind_ACALDt#134:                 dG_ACALDt#374 + 100000000 indicator_ACALDt#376
                                 <= 100000000
 ind_ACALDt_reverse_858fa#135:   dG_ACALDt_reverse_858fa#375
                                 + 100000000 indicator_ACALDt_reverse_858fa#377
                                 <= 100000000
 delG_ACALDt#136:                - 2.4776265 lnc_acald_c#230
                                 + 2.4776265 lnc_acald_e#232
                                 - dG_err_acald_c#302 + dG_err_acald_e#304
                                 + dG_ACALDt#374  = 10.0434871056239
 delG_ACALDt_reverse_858fa#137:  2.4776265 lnc_acald_c#230
                                 - 2.4776265 lnc_acald_e#232
                                 + dG_err_acald_c#302 - dG_err_acald_e#304
                                 + dG_ACALDt_reverse_858fa#375
                                  = -10.0434871056239
 directionality_ACKr#138:        ACKr#22 - 1000 indicator_ACKr#380 <= 0
 directionality_ACKr_reverse_b49c0#139: ACKr_reverse_b49c0#23
                                 - 1000 indicator_ACKr_reverse_b49c0#381 <= 0
 ind_ACKr#140:                   dG_ACKr#378 + 100000000 indicator_ACKr#380
                                 <= 100000000
 ind_ACKr_reverse_b49c0#141:     dG_ACKr_reverse_b49c0#379
                                 + 100000000 indicator_ACKr_reverse_b49c0#381
                                 <= 100000000
 delG_ACKr#142:                  2.4776265 lnc_ac_c#222
                                 - 2.4776265 lnc_actp_c#239
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244 + dG_err_ac_c#294
                                 - dG_err_actp_c#311 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 + dG_ACKr#378
                                  = 10.8315050400761
 delG_ACKr_reverse_b49c0#143:    - 2.4776265 lnc_ac_c#222
                                 + 2.4776265 lnc_actp_c#239
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244 - dG_err_ac_c#294
                                 + dG_err_actp_c#311 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_ACKr_reverse_b49c0#379
                                  = -10.8315050400761
 directionality_PPC#144:         PPC#24 - 1000 indicator_PPC#384 <= 0
 directionality_PPC_reverse_e854a#145: PPC_reverse_e854a#25
                                 - 1000 indicator_PPC_reverse_e854a#385 <= 0
 ind_PPC#146:                    dG_PPC#382 + 100000000 indicator_PPC#384
                                 <= 100000000
 ind_PPC_reverse_e854a#147:      dG_PPC_reverse_e854a#383
                                 + 100000000 indicator_PPC_reverse_e854a#385
                                 <= 100000000
 delG_PPC#148:                   2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_oaa_c#216
                                 + 2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_co2_c#246 + dG_err_h2o_c#268
                                 - dG_err_oaa_c#288 + dG_err_pep_c#289
                                 - dG_err_pi_c#291 + dG_err_co2_c#318
                                 + dG_PPC#382  = -40.2561391085733
 delG_PPC_reverse_e854a#149:     - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_oaa_c#216
                                 - 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_co2_c#246 - dG_err_h2o_c#268
                                 + dG_err_oaa_c#288 - dG_err_pep_c#289
                                 + dG_err_pi_c#291 - dG_err_co2_c#318
                                 + dG_PPC_reverse_e854a#383  = 40.2561391085733
 directionality_ACONTa#150:      ACONTa#26 - 1000 indicator_ACONTa#388 <= 0
 directionality_ACONTa_reverse_cad6d#151: ACONTa_reverse_cad6d#27
                                 - 1000 indicator_ACONTa_reverse_cad6d#389 <= 0
 ind_ACONTa#152:                 dG_ACONTa#386 + 100000000 indicator_ACONTa#388
                                 <= 100000000
 ind_ACONTa_reverse_cad6d#153:   dG_ACONTa_reverse_cad6d#387
                                 + 100000000 indicator_ACONTa_reverse_cad6d#389
                                 <= 100000000
 delG_ACONTa#154:                - 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_acon_C_c#237
                                 + 2.4776265 lnc_cit_c#245 - dG_err_h2o_c#268
                                 - dG_err_acon_C_c#309 + dG_err_cit_c#317
                                 + dG_ACONTa#386  = 7.00479011591341
 delG_ACONTa_reverse_cad6d#155:  2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_acon_C_c#237
                                 - 2.4776265 lnc_cit_c#245 + dG_err_h2o_c#268
                                 + dG_err_acon_C_c#309 - dG_err_cit_c#317
                                 + dG_ACONTa_reverse_cad6d#387
                                  = -7.00479011591341
 directionality_ACONTb#156:      ACONTb#28 - 1000 indicator_ACONTb#392 <= 0
 directionality_ACONTb_reverse_e198a#157: ACONTb_reverse_e198a#29
                                 - 1000 indicator_ACONTb_reverse_e198a#393 <= 0
 ind_ACONTb#158:                 dG_ACONTb#390 + 100000000 indicator_ACONTb#392
                                 <= 100000000
 ind_ACONTb_reverse_e198a#159:   dG_ACONTb_reverse_e198a#391
                                 + 100000000 indicator_ACONTb_reverse_e198a#393
                                 <= 100000000
 delG_ACONTb#160:                2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_icit_c#200
                                 + 2.4776265 lnc_acon_C_c#237
                                 + dG_err_h2o_c#268 - dG_err_icit_c#272
                                 + dG_err_acon_C_c#309 + dG_ACONTb#390
                                  = -1.86877955305329
 delG_ACONTb_reverse_e198a#161:  - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_icit_c#200
                                 - 2.4776265 lnc_acon_C_c#237
                                 - dG_err_h2o_c#268 + dG_err_icit_c#272
                                 - dG_err_acon_C_c#309
                                 + dG_ACONTb_reverse_e198a#391
                                  = 1.86877955305329
 directionality_ATPM#162:        ATPM#30 - 1000 indicator_ATPM#396 <= 0
 directionality_ATPM_reverse_5b752#163: ATPM_reverse_5b752#31
                                 - 1000 indicator_ATPM_reverse_5b752#397 <= 0
 ind_ATPM#164:                   dG_ATPM#394 + 100000000 indicator_ATPM#396
                                 <= 100000000
 ind_ATPM_reverse_5b752#165:     dG_ATPM_reverse_5b752#395
                                 + 100000000 indicator_ATPM_reverse_5b752#397
                                 <= 100000000
 delG_ATPM#166:                  2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244 + dG_err_h2o_c#268
                                 - dG_err_pi_c#291 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 + dG_ATPM#394
                                  = -32.3960040940278
 delG_ATPM_reverse_5b752#167:    - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244 - dG_err_h2o_c#268
                                 + dG_err_pi_c#291 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_ATPM_reverse_5b752#395
                                  = 32.3960040940278
 directionality_PPCK#168:        PPCK#32 - 1000 indicator_PPCK#400 <= 0
 directionality_PPCK_reverse_2557d#169: PPCK_reverse_2557d#33
                                 - 1000 indicator_PPCK_reverse_2557d#401 <= 0
 ind_PPCK#170:                   dG_PPCK#398 + 100000000 indicator_PPCK#400
                                 <= 100000000
 ind_PPCK_reverse_2557d#171:     dG_PPCK_reverse_2557d#399
                                 + 100000000 indicator_PPCK_reverse_2557d#401
                                 <= 100000000
 delG_PPCK#172:                  2.4776265 lnc_oaa_c#216
                                 - 2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244
                                 - 2.4776265 lnc_co2_c#246 + dG_err_oaa_c#288
                                 - dG_err_pep_c#289 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 - dG_err_co2_c#318
                                 + dG_PPCK#398  = 7.86013501454545
 delG_PPCK_reverse_2557d#173:    - 2.4776265 lnc_oaa_c#216
                                 + 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244
                                 + 2.4776265 lnc_co2_c#246 - dG_err_oaa_c#288
                                 + dG_err_pep_c#289 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_err_co2_c#318
                                 + dG_PPCK_reverse_2557d#399
                                  = -7.86013501454545
 directionality_ACt2r#174:       ACt2r#34 - 1000 indicator_ACt2r#404 <= 0
 directionality_ACt2r_reverse_64e05#175: ACt2r_reverse_64e05#35
                                 - 1000 indicator_ACt2r_reverse_64e05#405 <= 0
 ind_ACt2r#176:                  dG_ACt2r#402 + 100000000 indicator_ACt2r#404
                                 <= 100000000
 ind_ACt2r_reverse_64e05#177:    dG_ACt2r_reverse_64e05#403
                                 + 100000000 indicator_ACt2r_reverse_64e05#405
                                 <= 100000000
 delG_ACt2r#178:                 - 2.4776265 lnc_ac_c#222
                                 + 2.4776265 lnc_ac_e#229 - dG_err_ac_c#294
                                 + dG_err_ac_e#301 + dG_ACt2r#402
                                  = 5.5787946489177
 delG_ACt2r_reverse_64e05#179:   2.4776265 lnc_ac_c#222
                                 - 2.4776265 lnc_ac_e#229 + dG_err_ac_c#294
                                 - dG_err_ac_e#301 + dG_ACt2r_reverse_64e05#403
                                  = -5.5787946489177
 directionality_PPS#180:         PPS#36 - 1000 indicator_PPS#408 <= 0
 directionality_PPS_reverse_1c319#181: PPS_reverse_1c319#37
                                 - 1000 indicator_PPS_reverse_1c319#409 <= 0
 ind_PPS#182:                    dG_PPS#406 + 100000000 indicator_PPS#408
                                 <= 100000000
 ind_PPS_reverse_1c319#183:      dG_PPS_reverse_1c319#407
                                 + 100000000 indicator_PPS_reverse_1c319#409
                                 <= 100000000
 delG_PPS#184:                   2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_amp_c#243
                                 + 2.4776265 lnc_atp_c#244 + dG_err_h2o_c#268
                                 - dG_err_pep_c#289 - dG_err_pi_c#291
                                 + dG_err_pyr_c#295 - dG_err_amp_c#315
                                 + dG_err_atp_c#316 + dG_PPS#406
                                  = -8.52688867273173
 delG_PPS_reverse_1c319#185:     - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_amp_c#243
                                 - 2.4776265 lnc_atp_c#244 - dG_err_h2o_c#268
                                 + dG_err_pep_c#289 + dG_err_pi_c#291
                                 - dG_err_pyr_c#295 + dG_err_amp_c#315
                                 - dG_err_atp_c#316 + dG_PPS_reverse_1c319#407
                                  = 8.52688867273173
 directionality_ADK1#186:        ADK1#38 - 1000 indicator_ADK1#412 <= 0
 directionality_ADK1_reverse_a6f90#187: ADK1_reverse_a6f90#39
                                 - 1000 indicator_ADK1_reverse_a6f90#413 <= 0
 ind_ADK1#188:                   dG_ADK1#410 + 100000000 indicator_ADK1#412
                                 <= 100000000
 ind_ADK1_reverse_a6f90#189:     dG_ADK1_reverse_a6f90#411
                                 + 100000000 indicator_ADK1_reverse_a6f90#413
                                 <= 100000000
 delG_ADK1#190:                  - 4.955253 lnc_adp_c#240
                                 + 2.4776265 lnc_amp_c#243
                                 + 2.4776265 lnc_atp_c#244 - 2 dG_err_adp_c#312
                                 + dG_err_amp_c#315 + dG_err_atp_c#316
                                 + dG_ADK1#410  = -1.73015083791415
 delG_ADK1_reverse_a6f90#191:    4.955253 lnc_adp_c#240
                                 - 2.4776265 lnc_amp_c#243
                                 - 2.4776265 lnc_atp_c#244 + 2 dG_err_adp_c#312
                                 - dG_err_amp_c#315 - dG_err_atp_c#316
                                 + dG_ADK1_reverse_a6f90#411
                                  = 1.73015083791415
 directionality_AKGDH#192:       AKGDH#40 - 1000 indicator_AKGDH#416 <= 0
 directionality_AKGDH_reverse_08bdc#193: AKGDH_reverse_08bdc#41
                                 - 1000 indicator_AKGDH_reverse_08bdc#417 <= 0
 ind_AKGDH#194:                  dG_AKGDH#414 + 100000000 indicator_AKGDH#416
                                 <= 100000000
 ind_AKGDH_reverse_08bdc#195:    dG_AKGDH_reverse_08bdc#415
                                 + 100000000 indicator_AKGDH_reverse_08bdc#417
                                 <= 100000000
 delG_AKGDH#196:                 2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_succoa_c#236
                                 + 2.4776265 lnc_akg_c#241
                                 - 2.4776265 lnc_co2_c#246
                                 + 2.4776265 lnc_coa_c#248 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_succoa_c#308
                                 + dG_err_akg_c#313 - dG_err_co2_c#318
                                 + dG_err_coa_c#320 + dG_AKGDH#414
                                  = -26.3706432027718
 delG_AKGDH_reverse_08bdc#197:   - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_succoa_c#236
                                 - 2.4776265 lnc_akg_c#241
                                 + 2.4776265 lnc_co2_c#246
                                 - 2.4776265 lnc_coa_c#248 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_succoa_c#308
                                 - dG_err_akg_c#313 + dG_err_co2_c#318
                                 - dG_err_coa_c#320
                                 + dG_AKGDH_reverse_08bdc#415
                                  = 26.3706432027718
 directionality_ATPS4r#198:      ATPS4r#42 - 1000 indicator_ATPS4r#420 <= 0
 directionality_ATPS4r_reverse_64306#199: ATPS4r_reverse_64306#43
                                 - 1000 indicator_ATPS4r_reverse_64306#421 <= 0
 ind_ATPS4r#200:                 dG_ATPS4r#418 + 100000000 indicator_ATPS4r#420
                                 <= 100000000
 ind_ATPS4r_reverse_64306#201:   dG_ATPS4r_reverse_64306#419
                                 + 100000000 indicator_ATPS4r_reverse_64306#421
                                 <= 100000000
 delG_ATPS4r#202:                - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244 - dG_err_h2o_c#268
                                 + dG_err_pi_c#291 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_ATPS4r#418
                                  = 27.7908339080396
 delG_ATPS4r_reverse_64306#203:  2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244 + dG_err_h2o_c#268
                                 - dG_err_pi_c#291 - dG_err_adp_c#312
                                 + dG_err_atp_c#316
                                 + dG_ATPS4r_reverse_64306#419
                                  = -27.7908339080396
 directionality_PTAr#204:        PTAr#44 - 1000 indicator_PTAr#424 <= 0
 directionality_PTAr_reverse_fce15#205: PTAr_reverse_fce15#45
                                 - 1000 indicator_PTAr_reverse_fce15#425 <= 0
 ind_PTAr#206:                   dG_PTAr#422 + 100000000 indicator_PTAr#424
                                 <= 100000000
 ind_PTAr_reverse_fce15#207:     dG_PTAr_reverse_fce15#423
                                 + 100000000 indicator_PTAr_reverse_fce15#425
                                 <= 100000000
 delG_PTAr#208:                  2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_actp_c#239
                                 - 2.4776265 lnc_coa_c#248 + dG_err_pi_c#291
                                 + dG_err_accoa_c#305 - dG_err_actp_c#311
                                 - dG_err_coa_c#320 + dG_PTAr#422
                                  = 7.84455498813668
 delG_PTAr_reverse_fce15#209:    - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_actp_c#239
                                 + 2.4776265 lnc_coa_c#248 - dG_err_pi_c#291
                                 - dG_err_accoa_c#305 + dG_err_actp_c#311
                                 + dG_err_coa_c#320 + dG_PTAr_reverse_fce15#423
                                  = -7.84455498813668
 directionality_PYK#210:         PYK#46 - 1000 indicator_PYK#428 <= 0
 directionality_PYK_reverse_bc8ff#211: PYK_reverse_bc8ff#47
                                 - 1000 indicator_PYK_reverse_bc8ff#429 <= 0
 ind_PYK#212:                    dG_PYK#426 + 100000000 indicator_PYK#428
                                 <= 100000000
 ind_PYK_reverse_bc8ff#213:      dG_PYK_reverse_bc8ff#427
                                 + 100000000 indicator_PYK_reverse_bc8ff#429
                                 <= 100000000
 delG_PYK#214:                   2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244 + dG_err_pep_c#289
                                 - dG_err_pyr_c#295 + dG_err_adp_c#312
                                 - dG_err_atp_c#316 + dG_PYK#426
                                  = -22.1389645833819
 delG_PYK_reverse_bc8ff#215:     - 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244 - dG_err_pep_c#289
                                 + dG_err_pyr_c#295 - dG_err_adp_c#312
                                 + dG_err_atp_c#316 + dG_PYK_reverse_bc8ff#427
                                  = 22.1389645833819
 directionality_PYRt2#216:       PYRt2#50 - 1000 indicator_PYRt2#432 <= 0
 directionality_PYRt2_reverse_a6e24#217: PYRt2_reverse_a6e24#51
                                 - 1000 indicator_PYRt2_reverse_a6e24#433 <= 0
 ind_PYRt2#218:                  dG_PYRt2#430 + 100000000 indicator_PYRt2#432
                                 <= 100000000
 ind_PYRt2_reverse_a6e24#219:    dG_PYRt2_reverse_a6e24#431
                                 + 100000000 indicator_PYRt2_reverse_a6e24#433
                                 <= 100000000
 delG_PYRt2#220:                 - 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_pyr_e#224 - dG_err_pyr_c#295
                                 + dG_err_pyr_e#296 + dG_PYRt2#430
                                  = 5.57180761758205
 delG_PYRt2_reverse_a6e24#221:   2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_pyr_e#224 + dG_err_pyr_c#295
                                 - dG_err_pyr_e#296
                                 + dG_PYRt2_reverse_a6e24#431
                                  = -5.57180761758205
 directionality_CO2t#222:        CO2t#52 - 1000 indicator_CO2t#436 <= 0
 directionality_CO2t_reverse_7c42f#223: CO2t_reverse_7c42f#53
                                 - 1000 indicator_CO2t_reverse_7c42f#437 <= 0
 ind_CO2t#224:                   dG_CO2t#434 + 100000000 indicator_CO2t#436
                                 <= 100000000
 ind_CO2t_reverse_7c42f#225:     dG_CO2t_reverse_7c42f#435
                                 + 100000000 indicator_CO2t_reverse_7c42f#437
                                 <= 100000000
 delG_CO2t#226:                  - 2.4776265 lnc_co2_c#246
                                 + 2.4776265 lnc_co2_e#247 - dG_err_co2_c#318
                                 + dG_err_co2_e#319 + dG_CO2t#434  = 0
 delG_CO2t_reverse_7c42f#227:    2.4776265 lnc_co2_c#246
                                 - 2.4776265 lnc_co2_e#247 + dG_err_co2_c#318
                                 - dG_err_co2_e#319 + dG_CO2t_reverse_7c42f#435
                                  = 0
 directionality_RPE#228:         RPE#54 - 1000 indicator_RPE#440 <= 0
 directionality_RPE_reverse_a1b04#229: RPE_reverse_a1b04#55
                                 - 1000 indicator_RPE_reverse_a1b04#441 <= 0
 ind_RPE#230:                    dG_RPE#438 + 100000000 indicator_RPE#440
                                 <= 100000000
 ind_RPE_reverse_a1b04#231:      dG_RPE_reverse_a1b04#439
                                 + 100000000 indicator_RPE_reverse_a1b04#441
                                 <= 100000000
 delG_RPE#232:                   2.4776265 lnc_ru5p__D_c#228
                                 - 2.4776265 lnc_xu5p__D_c#238
                                 + dG_err_ru5p__D_c#300 - dG_err_xu5p__D_c#310
                                 + dG_RPE#438  = -3.38541847606643
 delG_RPE_reverse_a1b04#233:     - 2.4776265 lnc_ru5p__D_c#228
                                 + 2.4776265 lnc_xu5p__D_c#238
                                 - dG_err_ru5p__D_c#300 + dG_err_xu5p__D_c#310
                                 + dG_RPE_reverse_a1b04#439  = 3.38541847606643
 directionality_CS#234:          CS#56 - 1000 indicator_CS#444 <= 0
 directionality_CS_reverse_8d7e9#235: CS_reverse_8d7e9#57
                                 - 1000 indicator_CS_reverse_8d7e9#445 <= 0
 ind_CS#236:                     dG_CS#442 + 100000000 indicator_CS#444
                                 <= 100000000
 ind_CS_reverse_8d7e9#237:       dG_CS_reverse_8d7e9#443
                                 + 100000000 indicator_CS_reverse_8d7e9#445
                                 <= 100000000
 delG_CS#238:                    2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_oaa_c#216
                                 + 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_cit_c#245
                                 - 2.4776265 lnc_coa_c#248 + dG_err_h2o_c#268
                                 + dG_err_oaa_c#288 + dG_err_accoa_c#305
                                 - dG_err_cit_c#317 - dG_err_coa_c#320
                                 + dG_CS#442  = -37.522257555946
 delG_CS_reverse_8d7e9#239:      - 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_oaa_c#216
                                 - 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_cit_c#245
                                 + 2.4776265 lnc_coa_c#248 - dG_err_h2o_c#268
                                 - dG_err_oaa_c#288 - dG_err_accoa_c#305
                                 + dG_err_cit_c#317 + dG_err_coa_c#320
                                 + dG_CS_reverse_8d7e9#443  = 37.522257555946
 directionality_RPI#240:         RPI#58 - 1000 indicator_RPI#448 <= 0
 directionality_RPI_reverse_853a1#241: RPI_reverse_853a1#59
                                 - 1000 indicator_RPI_reverse_853a1#449 <= 0
 ind_RPI#242:                    dG_RPI#446 + 100000000 indicator_RPI#448
                                 <= 100000000
 ind_RPI_reverse_853a1#243:      dG_RPI_reverse_853a1#447
                                 + 100000000 indicator_RPI_reverse_853a1#449
                                 <= 100000000
 delG_RPI#244:                   2.4776265 lnc_r5p_c#227
                                 - 2.4776265 lnc_ru5p__D_c#228
                                 + dG_err_r5p_c#299 - dG_err_ru5p__D_c#300
                                 + dG_RPI#446  = 4.47815672235834
 delG_RPI_reverse_853a1#245:     - 2.4776265 lnc_r5p_c#227
                                 + 2.4776265 lnc_ru5p__D_c#228
                                 - dG_err_r5p_c#299 + dG_err_ru5p__D_c#300
                                 + dG_RPI_reverse_853a1#447
                                  = -4.47815672235834
 directionality_SUCCt2_2#246:    SUCCt2_2#60 - 1000 indicator_SUCCt2_2#452 <= 0
 directionality_SUCCt2_2_reverse_206ca#247: SUCCt2_2_reverse_206ca#61
                                 - 1000 indicator_SUCCt2_2_reverse_206ca#453
                                 <= 0
 ind_SUCCt2_2#248:               dG_SUCCt2_2#450
                                 + 100000000 indicator_SUCCt2_2#452
                                 <= 100000000
 ind_SUCCt2_2_reverse_206ca#249: dG_SUCCt2_2_reverse_206ca#451
                                 + 100000000 indicator_SUCCt2_2_reverse_206ca#453
                                 <= 100000000
 delG_SUCCt2_2#250:              - 2.4776265 lnc_succ_c#234
                                 + 2.4776265 lnc_succ_e#235 - dG_err_succ_c#306
                                 + dG_err_succ_e#307 + dG_SUCCt2_2#450
                                  = 4.6102863572564
 delG_SUCCt2_2_reverse_206ca#251: 2.4776265 lnc_succ_c#234
                                 - 2.4776265 lnc_succ_e#235 + dG_err_succ_c#306
                                 - dG_err_succ_e#307
                                 + dG_SUCCt2_2_reverse_206ca#451
                                  = -4.6102863572564
 directionality_CYTBD#252:       CYTBD#62 - 1000 indicator_CYTBD#456 <= 0
 directionality_CYTBD_reverse_611ba#253: CYTBD_reverse_611ba#63
                                 - 1000 indicator_CYTBD_reverse_611ba#457 <= 0
 ind_CYTBD#254:                  dG_CYTBD#454 + 100000000 indicator_CYTBD#456
                                 <= 100000000
 ind_CYTBD_reverse_611ba#255:    dG_CYTBD_reverse_611ba#455
                                 + 100000000 indicator_CYTBD_reverse_611ba#457
                                 <= 100000000
 delG_CYTBD#256:                 - 2.4776265 lnc_h2o_c#196
                                 + 1.23881325 lnc_o2_c#212
                                 - 2.4776265 lnc_q8_c#225
                                 + 2.4776265 lnc_q8h2_c#226 - dG_err_h2o_c#268
                                 + 0.5 dG_err_o2_c#284 - dG_err_q8_c#297
                                 + dG_err_q8h2_c#298 + dG_CYTBD#454
                                  = -91.2616478639669
 delG_CYTBD_reverse_611ba#257:   2.4776265 lnc_h2o_c#196
                                 - 1.23881325 lnc_o2_c#212
                                 + 2.4776265 lnc_q8_c#225
                                 - 2.4776265 lnc_q8h2_c#226 + dG_err_h2o_c#268
                                 - 0.5 dG_err_o2_c#284 + dG_err_q8_c#297
                                 - dG_err_q8h2_c#298
                                 + dG_CYTBD_reverse_611ba#455
                                  = 91.2616478639669
 directionality_D_LACt2#258:     D_LACt2#64 - 1000 indicator_D_LACt2#460 <= 0
 directionality_D_LACt2_reverse_d0dc8#259: D_LACt2_reverse_d0dc8#65
                                 - 1000 indicator_D_LACt2_reverse_d0dc8#461
                                 <= 0
 ind_D_LACt2#260:                dG_D_LACt2#458
                                 + 100000000 indicator_D_LACt2#460 <= 100000000
 ind_D_LACt2_reverse_d0dc8#261:  dG_D_LACt2_reverse_d0dc8#459
                                 + 100000000 indicator_D_LACt2_reverse_d0dc8#461
                                 <= 100000000
 delG_D_LACt2#262:               - 2.4776265 lnc_lac__D_c#201
                                 + 2.4776265 lnc_lac__D_e#202
                                 - dG_err_lac__D_c#273 + dG_err_lac__D_e#274
                                 + dG_D_LACt2#458  = 10.5946218476061
 delG_D_LACt2_reverse_d0dc8#263: 2.4776265 lnc_lac__D_c#201
                                 - 2.4776265 lnc_lac__D_e#202
                                 + dG_err_lac__D_c#273 - dG_err_lac__D_e#274
                                 + dG_D_LACt2_reverse_d0dc8#459
                                  = -10.5946218476061
 directionality_ENO#264:         _ENO#66 - 1000 indicator_ENO#464 <= 0
 directionality_ENO_reverse_40eea#265: _ENO_reverse_40eea#67
                                 - 1000 indicator_ENO_reverse_40eea#465 <= 0
 ind_ENO#266:                    dG_ENO#462 + 100000000 indicator_ENO#464
                                 <= 100000000
 ind_ENO_reverse_40eea#267:      dG_ENO_reverse_40eea#463
                                 + 100000000 indicator_ENO_reverse_40eea#465
                                 <= 100000000
 delG_ENO#268:                   - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_2pg_c#213
                                 - 2.4776265 lnc_pep_c#217 - dG_err_h2o_c#268
                                 + dG_err_2pg_c#285 - dG_err_pep_c#289
                                 + dG_ENO#462  = -3.70011313012651
 delG_ENO_reverse_40eea#269:     2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_2pg_c#213
                                 + 2.4776265 lnc_pep_c#217 + dG_err_h2o_c#268
                                 - dG_err_2pg_c#285 + dG_err_pep_c#289
                                 + dG_ENO_reverse_40eea#463  = 3.70011313012651
 directionality_SUCCt3#270:      SUCCt3#68 - 1000 indicator_SUCCt3#468 <= 0
 directionality_SUCCt3_reverse_7d82a#271: SUCCt3_reverse_7d82a#69
                                 - 1000 indicator_SUCCt3_reverse_7d82a#469 <= 0
 ind_SUCCt3#272:                 dG_SUCCt3#466 + 100000000 indicator_SUCCt3#468
                                 <= 100000000
 ind_SUCCt3_reverse_7d82a#273:   dG_SUCCt3_reverse_7d82a#467
                                 + 100000000 indicator_SUCCt3_reverse_7d82a#469
                                 <= 100000000
 delG_SUCCt3#274:                2.4776265 lnc_succ_c#234
                                 - 2.4776265 lnc_succ_e#235 + dG_err_succ_c#306
                                 - dG_err_succ_e#307 + dG_SUCCt3#466
                                  = -8.06416399674747
 delG_SUCCt3_reverse_7d82a#275:  - 2.4776265 lnc_succ_c#234
                                 + 2.4776265 lnc_succ_e#235 - dG_err_succ_c#306
                                 + dG_err_succ_e#307
                                 + dG_SUCCt3_reverse_7d82a#467
                                  = 8.06416399674747
 directionality_ETOHt2r#276:     _ETOHt2r#70 - 1000 indicator_ETOHt2r#472 <= 0
 directionality_ETOHt2r_reverse_b3d52#277: _ETOHt2r_reverse_b3d52#71
                                 - 1000 indicator_ETOHt2r_reverse_b3d52#473
                                 <= 0
 ind_ETOHt2r#278:                dG_ETOHt2r#470
                                 + 100000000 indicator_ETOHt2r#472 <= 100000000
 ind_ETOHt2r_reverse_b3d52#279:  dG_ETOHt2r_reverse_b3d52#471
                                 + 100000000 indicator_ETOHt2r_reverse_b3d52#473
                                 <= 100000000
 delG_ETOHt2r#280:               - 2.4776265 lnc_etoh_c#251
                                 + 2.4776265 lnc_etoh_e#252 - dG_err_etoh_c#323
                                 + dG_err_etoh_e#324 + dG_ETOHt2r#470
                                  = 13.9139381119389
 delG_ETOHt2r_reverse_b3d52#281: 2.4776265 lnc_etoh_c#251
                                 - 2.4776265 lnc_etoh_e#252 + dG_err_etoh_c#323
                                 - dG_err_etoh_e#324
                                 + dG_ETOHt2r_reverse_b3d52#471
                                  = -13.9139381119389
 directionality_SUCDi#282:       SUCDi#72 - 1000 indicator_SUCDi#476 <= 0
 directionality_SUCDi_reverse_480f4#283: SUCDi_reverse_480f4#73
                                 - 1000 indicator_SUCDi_reverse_480f4#477 <= 0
 ind_SUCDi#284:                  dG_SUCDi#474 + 100000000 indicator_SUCDi#476
                                 <= 100000000
 ind_SUCDi_reverse_480f4#285:    dG_SUCDi_reverse_480f4#475
                                 + 100000000 indicator_SUCDi_reverse_480f4#477
                                 <= 100000000
 delG_SUCDi#286:                 2.4776265 lnc_q8_c#225
                                 - 2.4776265 lnc_q8h2_c#226
                                 + 2.4776265 lnc_succ_c#234
                                 - 2.4776265 lnc_fum_c#258 + dG_err_q8_c#297
                                 - dG_err_q8h2_c#298 + dG_err_succ_c#306
                                 - dG_err_fum_c#330 + dG_SUCDi#474
                                  = -73.6344007799448
 delG_SUCDi_reverse_480f4#287:   - 2.4776265 lnc_q8_c#225
                                 + 2.4776265 lnc_q8h2_c#226
                                 - 2.4776265 lnc_succ_c#234
                                 + 2.4776265 lnc_fum_c#258 - dG_err_q8_c#297
                                 + dG_err_q8h2_c#298 - dG_err_succ_c#306
                                 + dG_err_fum_c#330
                                 + dG_SUCDi_reverse_480f4#475
                                  = 73.6344007799448
 directionality_SUCOAS#288:      SUCOAS#74 - 1000 indicator_SUCOAS#480 <= 0
 directionality_SUCOAS_reverse_22958#289: SUCOAS_reverse_22958#75
                                 - 1000 indicator_SUCOAS_reverse_22958#481 <= 0
 ind_SUCOAS#290:                 dG_SUCOAS#478 + 100000000 indicator_SUCOAS#480
                                 <= 100000000
 ind_SUCOAS_reverse_22958#291:   dG_SUCOAS_reverse_22958#479
                                 + 100000000 indicator_SUCOAS_reverse_22958#481
                                 <= 100000000
 delG_SUCOAS#292:                - 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_succ_c#234
                                 - 2.4776265 lnc_succoa_c#236
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244
                                 + 2.4776265 lnc_coa_c#248 - dG_err_pi_c#291
                                 + dG_err_succ_c#306 - dG_err_succoa_c#308
                                 - dG_err_adp_c#312 + dG_err_atp_c#316
                                 + dG_err_coa_c#320 + dG_SUCOAS#478
                                  = -2.03735147359816
 delG_SUCOAS_reverse_22958#293:  2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_succ_c#234
                                 + 2.4776265 lnc_succoa_c#236
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244
                                 - 2.4776265 lnc_coa_c#248 + dG_err_pi_c#291
                                 - dG_err_succ_c#306 + dG_err_succoa_c#308
                                 + dG_err_adp_c#312 - dG_err_atp_c#316
                                 - dG_err_coa_c#320
                                 + dG_SUCOAS_reverse_22958#479
                                  = 2.03735147359816
 directionality_TALA#294:        TALA#76 - 1000 indicator_TALA#484 <= 0
 directionality_TALA_reverse_adfda#295: TALA_reverse_adfda#77
                                 - 1000 indicator_TALA_reverse_adfda#485 <= 0
 ind_TALA#296:                   dG_TALA#482 + 100000000 indicator_TALA#484
                                 <= 100000000
 ind_TALA_reverse_adfda#297:     dG_TALA_reverse_adfda#483
                                 + 100000000 indicator_TALA_reverse_adfda#485
                                 <= 100000000
 delG_TALA#298:                  2.4776265 lnc_s7p_c#231
                                 - 2.4776265 lnc_e4p_c#250
                                 - 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_g3p_c#260 + dG_err_s7p_c#303
                                 - dG_err_e4p_c#322 - dG_err_f6p_c#325
                                 + dG_err_g3p_c#332 + dG_TALA#482
                                  = -0.600104512370081
 delG_TALA_reverse_adfda#299:    - 2.4776265 lnc_s7p_c#231
                                 + 2.4776265 lnc_e4p_c#250
                                 + 2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_g3p_c#260 - dG_err_s7p_c#303
                                 + dG_err_e4p_c#322 + dG_err_f6p_c#325
                                 - dG_err_g3p_c#332 + dG_TALA_reverse_adfda#483
                                  = 0.600104512370081
 directionality_THD2#300:        THD2#78 - 1000 indicator_THD2#488 <= 0
 directionality_THD2_reverse_f65dd#301: THD2_reverse_f65dd#79
                                 - 1000 indicator_THD2_reverse_f65dd#489 <= 0
 ind_THD2#302:                   dG_THD2#486 + 100000000 indicator_THD2#488
                                 <= 100000000
 ind_THD2_reverse_f65dd#303:     dG_THD2_reverse_f65dd#487
                                 + 100000000 indicator_THD2_reverse_f65dd#489
                                 <= 100000000
 delG_THD2#304:                  - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_nadp_c#279
                                 - dG_err_nadph_c#280 + dG_THD2#486
                                  = -2.03564395831553
 delG_THD2_reverse_f65dd#305:    2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_nadp_c#279
                                 + dG_err_nadph_c#280
                                 + dG_THD2_reverse_f65dd#487
                                  = 2.03564395831553
 directionality_TKT1#306:        TKT1#80 - 1000 indicator_TKT1#492 <= 0
 directionality_TKT1_reverse_a1021#307: TKT1_reverse_a1021#81
                                 - 1000 indicator_TKT1_reverse_a1021#493 <= 0
 ind_TKT1#308:                   dG_TKT1#490 + 100000000 indicator_TKT1#492
                                 <= 100000000
 ind_TKT1_reverse_a1021#309:     dG_TKT1_reverse_a1021#491
                                 + 100000000 indicator_TKT1_reverse_a1021#493
                                 <= 100000000
 delG_TKT1#310:                  2.4776265 lnc_r5p_c#227
                                 - 2.4776265 lnc_s7p_c#231
                                 + 2.4776265 lnc_xu5p__D_c#238
                                 - 2.4776265 lnc_g3p_c#260 + dG_err_r5p_c#299
                                 - dG_err_s7p_c#303 + dG_err_xu5p__D_c#310
                                 - dG_err_g3p_c#332 + dG_TKT1#490
                                  = -1.50702649191294
 delG_TKT1_reverse_a1021#311:    - 2.4776265 lnc_r5p_c#227
                                 + 2.4776265 lnc_s7p_c#231
                                 - 2.4776265 lnc_xu5p__D_c#238
                                 + 2.4776265 lnc_g3p_c#260 - dG_err_r5p_c#299
                                 + dG_err_s7p_c#303 - dG_err_xu5p__D_c#310
                                 + dG_err_g3p_c#332 + dG_TKT1_reverse_a1021#491
                                  = 1.50702649191294
 directionality_TKT2#312:        TKT2#82 - 1000 indicator_TKT2#496 <= 0
 directionality_TKT2_reverse_7ebc7#313: TKT2_reverse_7ebc7#83
                                 - 1000 indicator_TKT2_reverse_7ebc7#497 <= 0
 ind_TKT2#314:                   dG_TKT2#494 + 100000000 indicator_TKT2#496
                                 <= 100000000
 ind_TKT2_reverse_7ebc7#315:     dG_TKT2_reverse_7ebc7#495
                                 + 100000000 indicator_TKT2_reverse_7ebc7#497
                                 <= 100000000
 delG_TKT2#316:                  2.4776265 lnc_xu5p__D_c#238
                                 + 2.4776265 lnc_e4p_c#250
                                 - 2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_g3p_c#260
                                 + dG_err_xu5p__D_c#310 + dG_err_e4p_c#322
                                 - dG_err_f6p_c#325 - dG_err_g3p_c#332
                                 + dG_TKT2#494  = -9.97298988800389
 delG_TKT2_reverse_7ebc7#317:    - 2.4776265 lnc_xu5p__D_c#238
                                 - 2.4776265 lnc_e4p_c#250
                                 + 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_g3p_c#260
                                 - dG_err_xu5p__D_c#310 - dG_err_e4p_c#322
                                 + dG_err_f6p_c#325 + dG_err_g3p_c#332
                                 + dG_TKT2_reverse_7ebc7#495
                                  = 9.97298988800389
 directionality_TPI#318:         TPI#84 - 1000 indicator_TPI#500 <= 0
 directionality_TPI_reverse_c2c3b#319: TPI_reverse_c2c3b#85
                                 - 1000 indicator_TPI_reverse_c2c3b#501 <= 0
 ind_TPI#320:                    dG_TPI#498 + 100000000 indicator_TPI#500
                                 <= 100000000
 ind_TPI_reverse_c2c3b#321:      dG_TPI_reverse_c2c3b#499
                                 + 100000000 indicator_TPI_reverse_c2c3b#501
                                 <= 100000000
 delG_TPI#322:                   2.4776265 lnc_dhap_c#249
                                 - 2.4776265 lnc_g3p_c#260 + dG_err_dhap_c#321
                                 - dG_err_g3p_c#332 + dG_TPI#498
                                  = 5.57555793012284
 delG_TPI_reverse_c2c3b#323:     - 2.4776265 lnc_dhap_c#249
                                 + 2.4776265 lnc_g3p_c#260 - dG_err_dhap_c#321
                                 + dG_err_g3p_c#332 + dG_TPI_reverse_c2c3b#499
                                  = -5.57555793012284
 directionality_FBA#324:         FBA#126 - 1000 indicator_FBA#504 <= 0
 directionality_FBA_reverse_84806#325: FBA_reverse_84806#127
                                 - 1000 indicator_FBA_reverse_84806#505 <= 0
 ind_FBA#326:                    dG_FBA#502 + 100000000 indicator_FBA#504
                                 <= 100000000
 ind_FBA_reverse_84806#327:      dG_FBA_reverse_84806#503
                                 + 100000000 indicator_FBA_reverse_84806#505
                                 <= 100000000
 delG_FBA#328:                   - 2.4776265 lnc_dhap_c#249
                                 + 2.4776265 lnc_fdp_c#254
                                 - 2.4776265 lnc_g3p_c#260 - dG_err_dhap_c#321
                                 + dG_err_fdp_c#326 - dG_err_g3p_c#332
                                 + dG_FBA#502  = 21.5331182375942
 delG_FBA_reverse_84806#329:     2.4776265 lnc_dhap_c#249
                                 - 2.4776265 lnc_fdp_c#254
                                 + 2.4776265 lnc_g3p_c#260 + dG_err_dhap_c#321
                                 - dG_err_fdp_c#326 + dG_err_g3p_c#332
                                 + dG_FBA_reverse_84806#503
                                  = -21.5331182375942
 directionality_FBP#330:         FBP#128 - 1000 indicator_FBP#508 <= 0
 directionality_FBP_reverse_bf2c9#331: FBP_reverse_bf2c9#129
                                 - 1000 indicator_FBP_reverse_bf2c9#509 <= 0
 ind_FBP#332:                    dG_FBP#506 + 100000000 indicator_FBP#508
                                 <= 100000000
 ind_FBP_reverse_bf2c9#333:      dG_FBP_reverse_bf2c9#507
                                 + 100000000 indicator_FBP_reverse_bf2c9#509
                                 <= 100000000
 delG_FBP#334:                   2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_fdp_c#254 + dG_err_h2o_c#268
                                 - dG_err_pi_c#291 - dG_err_f6p_c#325
                                 + dG_err_fdp_c#326 + dG_FBP#506
                                  = -12.7352929572355
 delG_FBP_reverse_bf2c9#335:     - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_fdp_c#254 - dG_err_h2o_c#268
                                 + dG_err_pi_c#291 + dG_err_f6p_c#325
                                 - dG_err_fdp_c#326 + dG_FBP_reverse_bf2c9#507
                                  = 12.7352929572355
 directionality_FORt2#336:       FORt2#130 - 1000 indicator_FORt2#512 <= 0
 directionality_FORt2_reverse_89839#337: FORt2_reverse_89839#131
                                 - 1000 indicator_FORt2_reverse_89839#513 <= 0
 ind_FORt2#338:                  dG_FORt2#510 + 100000000 indicator_FORt2#512
                                 <= 100000000
 ind_FORt2_reverse_89839#339:    dG_FORt2_reverse_89839#511
                                 + 100000000 indicator_FORt2_reverse_89839#513
                                 <= 100000000
 delG_FORt2#340:                 - 2.4776265 lnc_for_c#255
                                 + 2.4776265 lnc_for_e#256 - dG_err_for_c#327
                                 + dG_err_for_e#328 + dG_FORt2#510
                                  = 0.553738322059918
 delG_FORt2_reverse_89839#341:   2.4776265 lnc_for_c#255
                                 - 2.4776265 lnc_for_e#256 + dG_err_for_c#327
                                 - dG_err_for_e#328
                                 + dG_FORt2_reverse_89839#511
                                  = -0.553738322059918
 directionality_FORt#342:        FORt#132 - 1000 indicator_FORt#516 <= 0
 directionality_FORt_reverse_40f9f#343: FORt_reverse_40f9f#133
                                 - 1000 indicator_FORt_reverse_40f9f#517 <= 0
 ind_FORt#344:                   dG_FORt#514 + 100000000 indicator_FORt#516
                                 <= 100000000
 ind_FORt_reverse_40f9f#345:     dG_FORt_reverse_40f9f#515
                                 + 100000000 indicator_FORt_reverse_40f9f#517
                                 <= 100000000
 delG_FORt#346:                  - 2.4776265 lnc_for_c#255
                                 + 2.4776265 lnc_for_e#256 - dG_err_for_c#327
                                 + dG_err_for_e#328 + dG_FORt#514
                                  = 1.70503086855694
 delG_FORt_reverse_40f9f#347:    2.4776265 lnc_for_c#255
                                 - 2.4776265 lnc_for_e#256 + dG_err_for_c#327
                                 - dG_err_for_e#328 + dG_FORt_reverse_40f9f#515
                                  = -1.70503086855694
 directionality_FRD7#348:        FRD7#134 - 1000 indicator_FRD7#520 <= 0
 directionality_FRD7_reverse_e0172#349: FRD7_reverse_e0172#135
                                 - 1000 indicator_FRD7_reverse_e0172#521 <= 0
 ind_FRD7#350:                   dG_FRD7#518 + 100000000 indicator_FRD7#520
                                 <= 100000000
 ind_FRD7_reverse_e0172#351:     dG_FRD7_reverse_e0172#519
                                 + 100000000 indicator_FRD7_reverse_e0172#521
                                 <= 100000000
 delG_FRD7#352:                  - 2.4776265 lnc_q8_c#225
                                 + 2.4776265 lnc_q8h2_c#226
                                 - 2.4776265 lnc_succ_c#234
                                 + 2.4776265 lnc_fum_c#258 - dG_err_q8_c#297
                                 + dG_err_q8h2_c#298 - dG_err_succ_c#306
                                 + dG_err_fum_c#330 + dG_FRD7#518
                                  = 73.6344007799448
 delG_FRD7_reverse_e0172#353:    2.4776265 lnc_q8_c#225
                                 - 2.4776265 lnc_q8h2_c#226
                                 + 2.4776265 lnc_succ_c#234
                                 - 2.4776265 lnc_fum_c#258 + dG_err_q8_c#297
                                 - dG_err_q8h2_c#298 + dG_err_succ_c#306
                                 - dG_err_fum_c#330 + dG_FRD7_reverse_e0172#519
                                  = -73.6344007799448
 directionality_FRUpts2#354:     FRUpts2#136 - 1000 indicator_FRUpts2#524 <= 0
 directionality_FRUpts2_reverse_58d3a#355: FRUpts2_reverse_58d3a#137
                                 - 1000 indicator_FRUpts2_reverse_58d3a#525
                                 <= 0
 ind_FRUpts2#356:                dG_FRUpts2#522
                                 + 100000000 indicator_FRUpts2#524 <= 100000000
 ind_FRUpts2_reverse_58d3a#357:  dG_FRUpts2_reverse_58d3a#523
                                 + 100000000 indicator_FRUpts2_reverse_58d3a#525
                                 <= 100000000
 delG_FRUpts2#358:               2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_f6p_c#253
                                 + 2.4776265 lnc_fru_e#257 + dG_err_pep_c#289
                                 - dG_err_pyr_c#295 - dG_err_f6p_c#325
                                 + dG_err_fru_e#329 + dG_FRUpts2#522
                                  = -12.511678579821
 delG_FRUpts2_reverse_58d3a#359: - 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_f6p_c#253
                                 - 2.4776265 lnc_fru_e#257 - dG_err_pep_c#289
                                 + dG_err_pyr_c#295 + dG_err_f6p_c#325
                                 - dG_err_fru_e#329
                                 + dG_FRUpts2_reverse_58d3a#523
                                  = 12.511678579821
 directionality_FUM#360:         FUM#138 - 1000 indicator_FUM#528 <= 0
 directionality_FUM_reverse_d3642#361: FUM_reverse_d3642#139
                                 - 1000 indicator_FUM_reverse_d3642#529 <= 0
 ind_FUM#362:                    dG_FUM#526 + 100000000 indicator_FUM#528
                                 <= 100000000
 ind_FUM_reverse_d3642#363:      dG_FUM_reverse_d3642#527
                                 + 100000000 indicator_FUM_reverse_d3642#529
                                 <= 100000000
 delG_FUM#364:                   2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_fum_c#258 + dG_err_h2o_c#268
                                 - dG_err_mal__L_c#275 + dG_err_fum_c#330
                                 + dG_FUM#526  = -3.42313126752583
 delG_FUM_reverse_d3642#365:     - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_fum_c#258 - dG_err_h2o_c#268
                                 + dG_err_mal__L_c#275 - dG_err_fum_c#330
                                 + dG_FUM_reverse_d3642#527  = 3.42313126752583
 directionality_FUMt2_2#366:     FUMt2_2#140 - 1000 indicator_FUMt2_2#532 <= 0
 directionality_FUMt2_2_reverse_7b4b6#367: FUMt2_2_reverse_7b4b6#141
                                 - 1000 indicator_FUMt2_2_reverse_7b4b6#533
                                 <= 0
 ind_FUMt2_2#368:                dG_FUMt2_2#530
                                 + 100000000 indicator_FUMt2_2#532 <= 100000000
 ind_FUMt2_2_reverse_7b4b6#369:  dG_FUMt2_2_reverse_7b4b6#531
                                 + 100000000 indicator_FUMt2_2_reverse_7b4b6#533
                                 <= 100000000
 delG_FUMt2_2#370:               - 2.4776265 lnc_fum_c#258
                                 + 2.4776265 lnc_fum_e#259 - dG_err_fum_c#330
                                 + dG_err_fum_e#331 + dG_FUMt2_2#530
                                  = -0.515849591595289
 delG_FUMt2_2_reverse_7b4b6#371: 2.4776265 lnc_fum_c#258
                                 - 2.4776265 lnc_fum_e#259 + dG_err_fum_c#330
                                 - dG_err_fum_e#331
                                 + dG_FUMt2_2_reverse_7b4b6#531
                                  = 0.515849591595289
 directionality_G6PDH2r#372:     G6PDH2r#142 - 1000 indicator_G6PDH2r#536 <= 0
 directionality_G6PDH2r_reverse_19ddf#373: G6PDH2r_reverse_19ddf#143
                                 - 1000 indicator_G6PDH2r_reverse_19ddf#537
                                 <= 0
 ind_G6PDH2r#374:                dG_G6PDH2r#534
                                 + 100000000 indicator_G6PDH2r#536 <= 100000000
 ind_G6PDH2r_reverse_19ddf#375:  dG_G6PDH2r_reverse_19ddf#535
                                 + 100000000 indicator_G6PDH2r_reverse_19ddf#537
                                 <= 100000000
 delG_G6PDH2r#376:               2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_6pgl_c#220
                                 + 2.4776265 lnc_g6p_c#261 + dG_err_nadp_c#279
                                 - dG_err_nadph_c#280 - dG_err_6pgl_c#292
                                 + dG_err_g6p_c#333 + dG_G6PDH2r#534
                                  = -6.68660901301746
 delG_G6PDH2r_reverse_19ddf#377: - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_6pgl_c#220
                                 - 2.4776265 lnc_g6p_c#261 - dG_err_nadp_c#279
                                 + dG_err_nadph_c#280 + dG_err_6pgl_c#292
                                 - dG_err_g6p_c#333
                                 + dG_G6PDH2r_reverse_19ddf#535
                                  = 6.68660901301746
 directionality_GAPD#378:        GAPD#144 - 1000 indicator_GAPD#540 <= 0
 directionality_GAPD_reverse_459c1#379: GAPD_reverse_459c1#145
                                 - 1000 indicator_GAPD_reverse_459c1#541 <= 0
 ind_GAPD#380:                   dG_GAPD#538 + 100000000 indicator_GAPD#540
                                 <= 100000000
 ind_GAPD_reverse_459c1#381:     dG_GAPD_reverse_459c1#539
                                 + 100000000 indicator_GAPD_reverse_459c1#541
                                 <= 100000000
 delG_GAPD#382:                  2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_13dpg_c#210
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_g3p_c#260 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_13dpg_c#282
                                 + dG_err_pi_c#291 + dG_err_g3p_c#332
                                 + dG_GAPD#538  = 2.34265549966563
 delG_GAPD_reverse_459c1#383:    - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_13dpg_c#210
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_g3p_c#260 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_13dpg_c#282
                                 - dG_err_pi_c#291 - dG_err_g3p_c#332
                                 + dG_GAPD_reverse_459c1#539
                                  = -2.34265549966563
 directionality_GLCpts#384:      GLCpts#146 - 1000 indicator_GLCpts#544 <= 0
 directionality_GLCpts_reverse_a52ae#385: GLCpts_reverse_a52ae#147
                                 - 1000 indicator_GLCpts_reverse_a52ae#545 <= 0
 ind_GLCpts#386:                 dG_GLCpts#542 + 100000000 indicator_GLCpts#544
                                 <= 100000000
 ind_GLCpts_reverse_a52ae#387:   dG_GLCpts_reverse_a52ae#543
                                 + 100000000 indicator_GLCpts_reverse_a52ae#545
                                 <= 100000000
 delG_GLCpts#388:                2.4776265 lnc_glc__D_e#190
                                 + 2.4776265 lnc_pep_c#217
                                 - 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_g6p_c#261
                                 + dG_err_glc__D_e#262 + dG_err_pep_c#289
                                 - dG_err_pyr_c#295 - dG_err_g6p_c#333
                                 + dG_GLCpts#542  = -15.5528259976698
 delG_GLCpts_reverse_a52ae#389:  - 2.4776265 lnc_glc__D_e#190
                                 - 2.4776265 lnc_pep_c#217
                                 + 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_g6p_c#261
                                 - dG_err_glc__D_e#262 - dG_err_pep_c#289
                                 + dG_err_pyr_c#295 + dG_err_g6p_c#333
                                 + dG_GLCpts_reverse_a52ae#543
                                  = 15.5528259976698
 directionality_GLNS#390:        GLNS#148 - 1000 indicator_GLNS#548 <= 0
 directionality_GLNS_reverse_59581#391: GLNS_reverse_59581#149
                                 - 1000 indicator_GLNS_reverse_59581#549 <= 0
 ind_GLNS#392:                   dG_GLNS#546 + 100000000 indicator_GLNS#548
                                 <= 100000000
 ind_GLNS_reverse_59581#393:     dG_GLNS_reverse_59581#547
                                 + 100000000 indicator_GLNS_reverse_59581#549
                                 <= 100000000
 delG_GLNS#394:                  - 2.4776265 lnc_gln__L_c#191
                                 + 2.4776265 lnc_glu__L_c#193
                                 + 2.4776265 lnc_nh4_c#209
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244
                                 - dG_err_gln__L_c#263 + dG_err_glu__L_c#265
                                 + dG_err_nh4_c#281 - dG_err_pi_c#291
                                 - dG_err_adp_c#312 + dG_err_atp_c#316
                                 + dG_GLNS#546  = -18.0548609709597
 delG_GLNS_reverse_59581#395:    2.4776265 lnc_gln__L_c#191
                                 - 2.4776265 lnc_glu__L_c#193
                                 - 2.4776265 lnc_nh4_c#209
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244
                                 + dG_err_gln__L_c#263 - dG_err_glu__L_c#265
                                 - dG_err_nh4_c#281 + dG_err_pi_c#291
                                 + dG_err_adp_c#312 - dG_err_atp_c#316
                                 + dG_GLNS_reverse_59581#547
                                  = 18.0548609709597
 directionality_GLNabc#396:      GLNabc#150 - 1000 indicator_GLNabc#552 <= 0
 directionality_GLNabc_reverse_1d82a#397: GLNabc_reverse_1d82a#151
                                 - 1000 indicator_GLNabc_reverse_1d82a#553 <= 0
 ind_GLNabc#398:                 dG_GLNabc#550 + 100000000 indicator_GLNabc#552
                                 <= 100000000
 ind_GLNabc_reverse_1d82a#399:   dG_GLNabc_reverse_1d82a#551
                                 + 100000000 indicator_GLNabc_reverse_1d82a#553
                                 <= 100000000
 delG_GLNabc#400:                - 2.4776265 lnc_gln__L_c#191
                                 + 2.4776265 lnc_gln__L_e#192
                                 + 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_pi_c#219
                                 - 2.4776265 lnc_adp_c#240
                                 + 2.4776265 lnc_atp_c#244
                                 - dG_err_gln__L_c#263 + dG_err_gln__L_e#264
                                 + dG_err_h2o_c#268 - dG_err_pi_c#291
                                 - dG_err_adp_c#312 + dG_err_atp_c#316
                                 + dG_GLNabc#550  = -7.32402354125946
 delG_GLNabc_reverse_1d82a#401:  2.4776265 lnc_gln__L_c#191
                                 - 2.4776265 lnc_gln__L_e#192
                                 - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_pi_c#219
                                 + 2.4776265 lnc_adp_c#240
                                 - 2.4776265 lnc_atp_c#244
                                 + dG_err_gln__L_c#263 - dG_err_gln__L_e#264
                                 - dG_err_h2o_c#268 + dG_err_pi_c#291
                                 + dG_err_adp_c#312 - dG_err_atp_c#316
                                 + dG_GLNabc_reverse_1d82a#551
                                  = 7.32402354125946
 directionality_GLUDy#402:       GLUDy#152 - 1000 indicator_GLUDy#556 <= 0
 directionality_GLUDy_reverse_fa4e7#403: GLUDy_reverse_fa4e7#153
                                 - 1000 indicator_GLUDy_reverse_fa4e7#557 <= 0
 ind_GLUDy#404:                  dG_GLUDy#554 + 100000000 indicator_GLUDy#556
                                 <= 100000000
 ind_GLUDy_reverse_fa4e7#405:    dG_GLUDy_reverse_fa4e7#555
                                 + 100000000 indicator_GLUDy_reverse_fa4e7#557
                                 <= 100000000
 delG_GLUDy#406:                 2.4776265 lnc_glu__L_c#193
                                 + 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_nh4_c#209
                                 - 2.4776265 lnc_akg_c#241
                                 + dG_err_glu__L_c#265 + dG_err_h2o_c#268
                                 + dG_err_nadp_c#279 - dG_err_nadph_c#280
                                 - dG_err_nh4_c#281 - dG_err_akg_c#313
                                 + dG_GLUDy#554  = 33.6313854677301
 delG_GLUDy_reverse_fa4e7#407:   - 2.4776265 lnc_glu__L_c#193
                                 - 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_nh4_c#209
                                 + 2.4776265 lnc_akg_c#241
                                 - dG_err_glu__L_c#265 - dG_err_h2o_c#268
                                 - dG_err_nadp_c#279 + dG_err_nadph_c#280
                                 + dG_err_nh4_c#281 + dG_err_akg_c#313
                                 + dG_GLUDy_reverse_fa4e7#555
                                  = -33.6313854677301
 directionality_GLUN#408:        GLUN#154 - 1000 indicator_GLUN#560 <= 0
 directionality_GLUN_reverse_4ccdb#409: GLUN_reverse_4ccdb#155
                                 - 1000 indicator_GLUN_reverse_4ccdb#561 <= 0
 ind_GLUN#410:                   dG_GLUN#558 + 100000000 indicator_GLUN#560
                                 <= 100000000
 ind_GLUN_reverse_4ccdb#411:     dG_GLUN_reverse_4ccdb#559
                                 + 100000000 indicator_GLUN_reverse_4ccdb#561
                                 <= 100000000
 delG_GLUN#412:                  2.4776265 lnc_gln__L_c#191
                                 - 2.4776265 lnc_glu__L_c#193
                                 + 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_nh4_c#209
                                 + dG_err_gln__L_c#263 - dG_err_glu__L_c#265
                                 + dG_err_h2o_c#268 - dG_err_nh4_c#281
                                 + dG_GLUN#558  = -14.3411431230676
 delG_GLUN_reverse_4ccdb#413:    - 2.4776265 lnc_gln__L_c#191
                                 + 2.4776265 lnc_glu__L_c#193
                                 - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_nh4_c#209
                                 - dG_err_gln__L_c#263 + dG_err_glu__L_c#265
                                 - dG_err_h2o_c#268 + dG_err_nh4_c#281
                                 + dG_GLUN_reverse_4ccdb#559
                                  = 14.3411431230676
 directionality_GLUSy#414:       GLUSy#156 - 1000 indicator_GLUSy#564 <= 0
 directionality_GLUSy_reverse_6a00f#415: GLUSy_reverse_6a00f#157
                                 - 1000 indicator_GLUSy_reverse_6a00f#565 <= 0
 ind_GLUSy#416:                  dG_GLUSy#562 + 100000000 indicator_GLUSy#564
                                 <= 100000000
 ind_GLUSy_reverse_6a00f#417:    dG_GLUSy_reverse_6a00f#563
                                 + 100000000 indicator_GLUSy_reverse_6a00f#565
                                 <= 100000000
 delG_GLUSy#418:                 2.4776265 lnc_gln__L_c#191
                                 - 4.955253 lnc_glu__L_c#193
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_akg_c#241
                                 + dG_err_gln__L_c#263 - 2 dG_err_glu__L_c#265
                                 - dG_err_nadp_c#279 + dG_err_nadph_c#280
                                 + dG_err_akg_c#313 + dG_GLUSy#562
                                  = -47.9725285907978
 delG_GLUSy_reverse_6a00f#419:   - 2.4776265 lnc_gln__L_c#191
                                 + 4.955253 lnc_glu__L_c#193
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_akg_c#241
                                 - dG_err_gln__L_c#263 + 2 dG_err_glu__L_c#265
                                 + dG_err_nadp_c#279 - dG_err_nadph_c#280
                                 - dG_err_akg_c#313
                                 + dG_GLUSy_reverse_6a00f#563
                                  = 47.9725285907978
 directionality_GLUt2r#420:      GLUt2r#158 - 1000 indicator_GLUt2r#568 <= 0
 directionality_GLUt2r_reverse_3e88a#421: GLUt2r_reverse_3e88a#159
                                 - 1000 indicator_GLUt2r_reverse_3e88a#569 <= 0
 ind_GLUt2r#422:                 dG_GLUt2r#566 + 100000000 indicator_GLUt2r#568
                                 <= 100000000
 ind_GLUt2r_reverse_3e88a#423:   dG_GLUt2r_reverse_3e88a#567
                                 + 100000000 indicator_GLUt2r_reverse_3e88a#569
                                 <= 100000000
 delG_GLUt2r#424:                - 2.4776265 lnc_glu__L_c#193
                                 + 2.4776265 lnc_glu__L_e#194
                                 - dG_err_glu__L_c#265 + dG_err_glu__L_e#266
                                 + dG_GLUt2r#566  = 18.0548523661902
 delG_GLUt2r_reverse_3e88a#425:  2.4776265 lnc_glu__L_c#193
                                 - 2.4776265 lnc_glu__L_e#194
                                 + dG_err_glu__L_c#265 - dG_err_glu__L_e#266
                                 + dG_GLUt2r_reverse_3e88a#567
                                  = -18.0548523661902
 directionality_GND#426:         GND#160 - 1000 indicator_GND#572 <= 0
 directionality_GND_reverse_eec5c#427: GND_reverse_eec5c#161
                                 - 1000 indicator_GND_reverse_eec5c#573 <= 0
 ind_GND#428:                    dG_GND#570 + 100000000 indicator_GND#572
                                 <= 100000000
 ind_GND_reverse_eec5c#429:      dG_GND_reverse_eec5c#571
                                 + 100000000 indicator_GND_reverse_eec5c#573
                                 <= 100000000
 delG_GND#430:                   2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_6pgc_c#218
                                 - 2.4776265 lnc_ru5p__D_c#228
                                 - 2.4776265 lnc_co2_c#246 + dG_err_nadp_c#279
                                 - dG_err_nadph_c#280 + dG_err_6pgc_c#290
                                 - dG_err_ru5p__D_c#300 - dG_err_co2_c#318
                                 + dG_GND#570  = 10.4135511702652
 delG_GND_reverse_eec5c#431:     - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_6pgc_c#218
                                 + 2.4776265 lnc_ru5p__D_c#228
                                 + 2.4776265 lnc_co2_c#246 - dG_err_nadp_c#279
                                 + dG_err_nadph_c#280 - dG_err_6pgc_c#290
                                 + dG_err_ru5p__D_c#300 + dG_err_co2_c#318
                                 + dG_GND_reverse_eec5c#571
                                  = -10.4135511702652
 directionality_ICDHyr#432:      ICDHyr#164 - 1000 indicator_ICDHyr#576 <= 0
 directionality_ICDHyr_reverse_7f84b#433: ICDHyr_reverse_7f84b#165
                                 - 1000 indicator_ICDHyr_reverse_7f84b#577 <= 0
 ind_ICDHyr#434:                 dG_ICDHyr#574 + 100000000 indicator_ICDHyr#576
                                 <= 100000000
 ind_ICDHyr_reverse_7f84b#435:   dG_ICDHyr_reverse_7f84b#575
                                 + 100000000 indicator_ICDHyr_reverse_7f84b#577
                                 <= 100000000
 delG_ICDHyr#436:                2.4776265 lnc_icit_c#200
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_akg_c#241
                                 - 2.4776265 lnc_co2_c#246 + dG_err_icit_c#272
                                 + dG_err_nadp_c#279 - dG_err_nadph_c#280
                                 - dG_err_akg_c#313 - dG_err_co2_c#318
                                 + dG_ICDHyr#574  = 5.59903791423471
 delG_ICDHyr_reverse_7f84b#437:  - 2.4776265 lnc_icit_c#200
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_akg_c#241
                                 + 2.4776265 lnc_co2_c#246 - dG_err_icit_c#272
                                 - dG_err_nadp_c#279 + dG_err_nadph_c#280
                                 + dG_err_akg_c#313 + dG_err_co2_c#318
                                 + dG_ICDHyr_reverse_7f84b#575
                                  = -5.59903791423471
 directionality_ICL#438:         ICL#166 - 1000 indicator_ICL#580 <= 0
 directionality_ICL_reverse_2f27e#439: ICL_reverse_2f27e#167
                                 - 1000 indicator_ICL_reverse_2f27e#581 <= 0
 ind_ICL#440:                    dG_ICL#578 + 100000000 indicator_ICL#580
                                 <= 100000000
 ind_ICL_reverse_2f27e#441:      dG_ICL_reverse_2f27e#579
                                 + 100000000 indicator_ICL_reverse_2f27e#581
                                 <= 100000000
 delG_ICL#442:                   - 2.4776265 lnc_glx_c#195
                                 + 2.4776265 lnc_icit_c#200
                                 - 2.4776265 lnc_succ_c#234 - dG_err_glx_c#267
                                 + dG_err_icit_c#272 - dG_err_succ_c#306
                                 + dG_ICL#578  = 9.54566792621745
 delG_ICL_reverse_2f27e#443:     2.4776265 lnc_glx_c#195
                                 - 2.4776265 lnc_icit_c#200
                                 + 2.4776265 lnc_succ_c#234 + dG_err_glx_c#267
                                 - dG_err_icit_c#272 + dG_err_succ_c#306
                                 + dG_ICL_reverse_2f27e#579
                                  = -9.54566792621745
 directionality_LDH_D#444:       LDH_D#168 - 1000 indicator_LDH_D#584 <= 0
 directionality_LDH_D_reverse_f8507#445: LDH_D_reverse_f8507#169
                                 - 1000 indicator_LDH_D_reverse_f8507#585 <= 0
 ind_LDH_D#446:                  dG_LDH_D#582 + 100000000 indicator_LDH_D#584
                                 <= 100000000
 ind_LDH_D_reverse_f8507#447:    dG_LDH_D_reverse_f8507#583
                                 + 100000000 indicator_LDH_D_reverse_f8507#585
                                 <= 100000000
 delG_LDH_D#448:                 2.4776265 lnc_lac__D_c#201
                                 + 2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_pyr_c#223
                                 + dG_err_lac__D_c#273 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_pyr_c#295
                                 + dG_LDH_D#582  = 20.7139323940812
 delG_LDH_D_reverse_f8507#449:   - 2.4776265 lnc_lac__D_c#201
                                 - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_pyr_c#223
                                 - dG_err_lac__D_c#273 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_pyr_c#295
                                 + dG_LDH_D_reverse_f8507#583
                                  = -20.7139323940812
 directionality_MALS#450:        MALS#170 - 1000 indicator_MALS#588 <= 0
 directionality_MALS_reverse_d7382#451: MALS_reverse_d7382#171
                                 - 1000 indicator_MALS_reverse_d7382#589 <= 0
 ind_MALS#452:                   dG_MALS#586 + 100000000 indicator_MALS#588
                                 <= 100000000
 ind_MALS_reverse_d7382#453:     dG_MALS_reverse_d7382#587
                                 + 100000000 indicator_MALS_reverse_d7382#589
                                 <= 100000000
 delG_MALS#454:                  2.4776265 lnc_glx_c#195
                                 + 2.4776265 lnc_h2o_c#196
                                 - 2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_coa_c#248 + dG_err_glx_c#267
                                 + dG_err_h2o_c#268 - dG_err_mal__L_c#275
                                 + dG_err_accoa_c#305 - dG_err_coa_c#320
                                 + dG_MALS#586  = -39.0717649751016
 delG_MALS_reverse_d7382#455:    - 2.4776265 lnc_glx_c#195
                                 - 2.4776265 lnc_h2o_c#196
                                 + 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_coa_c#248 - dG_err_glx_c#267
                                 - dG_err_h2o_c#268 + dG_err_mal__L_c#275
                                 - dG_err_accoa_c#305 + dG_err_coa_c#320
                                 + dG_MALS_reverse_d7382#587
                                  = 39.0717649751016
 directionality_MALt2_2#456:     MALt2_2#172 - 1000 indicator_MALt2_2#592 <= 0
 directionality_MALt2_2_reverse_a635f#457: MALt2_2_reverse_a635f#173
                                 - 1000 indicator_MALt2_2_reverse_a635f#593
                                 <= 0
 ind_MALt2_2#458:                dG_MALt2_2#590
                                 + 100000000 indicator_MALt2_2#592 <= 100000000
 ind_MALt2_2_reverse_a635f#459:  dG_MALt2_2_reverse_a635f#591
                                 + 100000000 indicator_MALt2_2_reverse_a635f#593
                                 <= 100000000
 delG_MALt2_2#460:               - 2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_mal__L_e#204
                                 - dG_err_mal__L_c#275 + dG_err_mal__L_e#276
                                 + dG_MALt2_2#590  = 4.53248234665645
 delG_MALt2_2_reverse_a635f#461: 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_mal__L_e#204
                                 + dG_err_mal__L_c#275 - dG_err_mal__L_e#276
                                 + dG_MALt2_2_reverse_a635f#591
                                  = -4.53248234665645
 directionality_MDH#462:         MDH#174 - 1000 indicator_MDH#596 <= 0
 directionality_MDH_reverse_ee52c#463: MDH_reverse_ee52c#175
                                 - 1000 indicator_MDH_reverse_ee52c#597 <= 0
 ind_MDH#464:                    dG_MDH#594 + 100000000 indicator_MDH#596
                                 <= 100000000
 ind_MDH_reverse_ee52c#465:      dG_MDH_reverse_ee52c#595
                                 + 100000000 indicator_MDH_reverse_ee52c#597
                                 <= 100000000
 delG_MDH#466:                   2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_oaa_c#216
                                 + dG_err_mal__L_c#275 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_oaa_c#288
                                 + dG_MDH#594  = 26.4730436873124
 delG_MDH_reverse_ee52c#467:     - 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_oaa_c#216
                                 - dG_err_mal__L_c#275 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_oaa_c#288
                                 + dG_MDH_reverse_ee52c#595
                                  = -26.4730436873124
 directionality_ME1#468:         ME1#176 - 1000 indicator_ME1#600 <= 0
 directionality_ME1_reverse_9736c#469: ME1_reverse_9736c#177
                                 - 1000 indicator_ME1_reverse_9736c#601 <= 0
 ind_ME1#470:                    dG_ME1#598 + 100000000 indicator_ME1#600
                                 <= 100000000
 ind_ME1_reverse_9736c#471:      dG_ME1_reverse_9736c#599
                                 + 100000000 indicator_ME1_reverse_9736c#601
                                 <= 100000000
 delG_ME1#472:                   2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_co2_c#246
                                 + dG_err_mal__L_c#275 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_pyr_c#295
                                 - dG_err_co2_c#318 + dG_ME1#598
                                  = 12.1942141184761
 delG_ME1_reverse_9736c#473:     - 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_co2_c#246
                                 - dG_err_mal__L_c#275 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_pyr_c#295
                                 + dG_err_co2_c#318 + dG_ME1_reverse_9736c#599
                                  = -12.1942141184761
 directionality_ME2#474:         ME2#178 - 1000 indicator_ME2#604 <= 0
 directionality_ME2_reverse_2b0a2#475: ME2_reverse_2b0a2#179
                                 - 1000 indicator_ME2_reverse_2b0a2#605 <= 0
 ind_ME2#476:                    dG_ME2#602 + 100000000 indicator_ME2#604
                                 <= 100000000
 ind_ME2_reverse_2b0a2#477:      dG_ME2_reverse_2b0a2#603
                                 + 100000000 indicator_ME2_reverse_2b0a2#605
                                 <= 100000000
 delG_ME2#478:                   2.4776265 lnc_mal__L_c#203
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208
                                 - 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_co2_c#246
                                 + dG_err_mal__L_c#275 + dG_err_nadp_c#279
                                 - dG_err_nadph_c#280 - dG_err_pyr_c#295
                                 - dG_err_co2_c#318 + dG_ME2#602
                                  = 12.4611552531544
 delG_ME2_reverse_2b0a2#479:     - 2.4776265 lnc_mal__L_c#203
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208
                                 + 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_co2_c#246
                                 - dG_err_mal__L_c#275 - dG_err_nadp_c#279
                                 + dG_err_nadph_c#280 + dG_err_pyr_c#295
                                 + dG_err_co2_c#318 + dG_ME2_reverse_2b0a2#603
                                  = -12.4611552531544
 directionality_NADH16#480:      NADH16#180 - 1000 indicator_NADH16#608 <= 0
 directionality_NADH16_reverse_330a2#481: NADH16_reverse_330a2#181
                                 - 1000 indicator_NADH16_reverse_330a2#609 <= 0
 ind_NADH16#482:                 dG_NADH16#606 + 100000000 indicator_NADH16#608
                                 <= 100000000
 ind_NADH16_reverse_330a2#483:   dG_NADH16_reverse_330a2#607
                                 + 100000000 indicator_NADH16_reverse_330a2#609
                                 <= 100000000
 delG_NADH16#484:                - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_q8_c#225
                                 - 2.4776265 lnc_q8h2_c#226 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_q8_c#297
                                 - dG_err_q8h2_c#298 + dG_NADH16#606
                                  = -128.220601732739
 delG_NADH16_reverse_330a2#485:  2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_q8_c#225
                                 + 2.4776265 lnc_q8h2_c#226 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_q8_c#297
                                 + dG_err_q8h2_c#298
                                 + dG_NADH16_reverse_330a2#607
                                  = 128.220601732739
 directionality_NADTRHD#486:     NADTRHD#182 - 1000 indicator_NADTRHD#612 <= 0
 directionality_NADTRHD_reverse_49725#487: NADTRHD_reverse_49725#183
                                 - 1000 indicator_NADTRHD_reverse_49725#613
                                 <= 0
 ind_NADTRHD#488:                dG_NADTRHD#610
                                 + 100000000 indicator_NADTRHD#612 <= 100000000
 ind_NADTRHD_reverse_49725#489:  dG_NADTRHD_reverse_49725#611
                                 + 100000000 indicator_NADTRHD_reverse_49725#613
                                 <= 100000000
 delG_NADTRHD#490:               2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_nadp_c#207
                                 + 2.4776265 lnc_nadph_c#208 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 - dG_err_nadp_c#279
                                 + dG_err_nadph_c#280 + dG_NADTRHD#610
                                  = -0.266941134678291
 delG_NADTRHD_reverse_49725#491: - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_nadp_c#207
                                 - 2.4776265 lnc_nadph_c#208 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 + dG_err_nadp_c#279
                                 - dG_err_nadph_c#280
                                 + dG_NADTRHD_reverse_49725#611
                                  = 0.266941134678291
 directionality_NH4t#492:        NH4t#184 - 1000 indicator_NH4t#616 <= 0
 directionality_NH4t_reverse_551ee#493: NH4t_reverse_551ee#185
                                 - 1000 indicator_NH4t_reverse_551ee#617 <= 0
 ind_NH4t#494:                   dG_NH4t#614 + 100000000 indicator_NH4t#616
                                 <= 100000000
 ind_NH4t_reverse_551ee#495:     dG_NH4t_reverse_551ee#615
                                 + 100000000 indicator_NH4t_reverse_551ee#617
                                 <= 100000000
 delG_NH4t#496:                  - 2.4776265 lnc_nh4_c#209
                                 + 2.4776265 lnc_nh4_e#211 - dG_err_nh4_c#281
                                 + dG_err_nh4_e#283 + dG_NH4t#614
                                  = 9.16190424097423
 delG_NH4t_reverse_551ee#497:    2.4776265 lnc_nh4_c#209
                                 - 2.4776265 lnc_nh4_e#211 + dG_err_nh4_c#281
                                 - dG_err_nh4_e#283 + dG_NH4t_reverse_551ee#615
                                  = -9.16190424097423
 directionality_PDH#498:         PDH#188 - 1000 indicator_PDH#620 <= 0
 directionality_PDH_reverse_ca160#499: PDH_reverse_ca160#189
                                 - 1000 indicator_PDH_reverse_ca160#621 <= 0
 ind_PDH#500:                    dG_PDH#618 + 100000000 indicator_PDH#620
                                 <= 100000000
 ind_PDH_reverse_ca160#501:      dG_PDH_reverse_ca160#619
                                 + 100000000 indicator_PDH_reverse_ca160#621
                                 <= 100000000
 delG_PDH#502:                   2.4776265 lnc_nad_c#205
                                 - 2.4776265 lnc_nadh_c#206
                                 + 2.4776265 lnc_pyr_c#223
                                 - 2.4776265 lnc_accoa_c#233
                                 - 2.4776265 lnc_co2_c#246
                                 + 2.4776265 lnc_coa_c#248 + dG_err_nad_c#277
                                 - dG_err_nadh_c#278 + dG_err_pyr_c#295
                                 - dG_err_accoa_c#305 - dG_err_co2_c#318
                                 + dG_err_coa_c#320 + dG_PDH#618
                                  = -34.0653161132366
 delG_PDH_reverse_ca160#503:     - 2.4776265 lnc_nad_c#205
                                 + 2.4776265 lnc_nadh_c#206
                                 - 2.4776265 lnc_pyr_c#223
                                 + 2.4776265 lnc_accoa_c#233
                                 + 2.4776265 lnc_co2_c#246
                                 - 2.4776265 lnc_coa_c#248 - dG_err_nad_c#277
                                 + dG_err_nadh_c#278 - dG_err_pyr_c#295
                                 + dG_err_accoa_c#305 + dG_err_co2_c#318
                                 - dG_err_coa_c#320 + dG_PDH_reverse_ca160#619
                                  = 34.0653161132366
Bounds
 0 <= PFK#0 <= 1000
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
-11.5129254649702 <= lnc_glc__D_e#190 <= -3.91202300542815
-11.5129254649702 <= lnc_gln__L_c#191 <= -3.91202300542815
-11.5129254649702 <= lnc_gln__L_e#192 <= -3.91202300542815
-11.5129254649702 <= lnc_glu__L_c#193 <= -3.91202300542815
-11.5129254649702 <= lnc_glu__L_e#194 <= -3.91202300542815
-11.5129254649702 <= lnc_glx_c#195 <= -3.91202300542815
-11.5129254649702 <= lnc_h2o_c#196 <= -3.91202300542815
-11.5129254649702 <= lnc_h2o_e#197 <= -3.91202300542815
-11.5129254649702 <= lnc_h_c#198 <= -3.91202300542815
-11.5129254649702 <= lnc_h_e#199 <= -3.91202300542815
-11.5129254649702 <= lnc_icit_c#200 <= -3.91202300542815
-11.5129254649702 <= lnc_lac__D_c#201 <= -3.91202300542815
-11.5129254649702 <= lnc_lac__D_e#202 <= -3.91202300542815
-11.5129254649702 <= lnc_mal__L_c#203 <= -3.91202300542815
-11.5129254649702 <= lnc_mal__L_e#204 <= -3.91202300542815
-11.5129254649702 <= lnc_nad_c#205 <= -3.91202300542815
-11.5129254649702 <= lnc_nadh_c#206 <= -3.91202300542815
-11.5129254649702 <= lnc_nadp_c#207 <= -3.91202300542815
-11.5129254649702 <= lnc_nadph_c#208 <= -3.91202300542815
-11.5129254649702 <= lnc_nh4_c#209 <= -3.91202300542815
-11.5129254649702 <= lnc_13dpg_c#210 <= -3.91202300542815
-11.5129254649702 <= lnc_nh4_e#211 <= -3.91202300542815
-11.5129254649702 <= lnc_o2_c#212 <= -3.91202300542815
-11.5129254649702 <= lnc_2pg_c#213 <= -3.91202300542815
-11.5129254649702 <= lnc_o2_e#214 <= -3.91202300542815
-11.5129254649702 <= lnc_3pg_c#215 <= -3.91202300542815
-11.5129254649702 <= lnc_oaa_c#216 <= -3.91202300542815
-11.5129254649702 <= lnc_pep_c#217 <= -3.91202300542815
-11.5129254649702 <= lnc_6pgc_c#218 <= -3.91202300542815
-11.5129254649702 <= lnc_pi_c#219 <= -3.91202300542815
-11.5129254649702 <= lnc_6pgl_c#220 <= -3.91202300542815
-11.5129254649702 <= lnc_pi_e#221 <= -3.91202300542815
-11.5129254649702 <= lnc_ac_c#222 <= -3.91202300542815
-11.5129254649702 <= lnc_pyr_c#223 <= -3.91202300542815
-11.5129254649702 <= lnc_pyr_e#224 <= -3.91202300542815
-11.5129254649702 <= lnc_q8_c#225 <= -3.91202300542815
-11.5129254649702 <= lnc_q8h2_c#226 <= -3.91202300542815
-11.5129254649702 <= lnc_r5p_c#227 <= -3.91202300542815
-11.5129254649702 <= lnc_ru5p__D_c#228 <= -3.91202300542815
-11.5129254649702 <= lnc_ac_e#229 <= -3.91202300542815
-11.5129254649702 <= lnc_acald_c#230 <= -3.91202300542815
-11.5129254649702 <= lnc_s7p_c#231 <= -3.91202300542815
-11.5129254649702 <= lnc_acald_e#232 <= -3.91202300542815
-11.5129254649702 <= lnc_accoa_c#233 <= -3.91202300542815
-11.5129254649702 <= lnc_succ_c#234 <= -3.91202300542815
-11.5129254649702 <= lnc_succ_e#235 <= -3.91202300542815
-11.5129254649702 <= lnc_succoa_c#236 <= -3.91202300542815
-11.5129254649702 <= lnc_acon_C_c#237 <= -3.91202300542815
-11.5129254649702 <= lnc_xu5p__D_c#238 <= -3.91202300542815
-11.5129254649702 <= lnc_actp_c#239 <= -3.91202300542815
-11.5129254649702 <= lnc_adp_c#240 <= -3.91202300542815
-11.5129254649702 <= lnc_akg_c#241 <= -3.91202300542815
-11.5129254649702 <= lnc_akg_e#242 <= -3.91202300542815
-11.5129254649702 <= lnc_amp_c#243 <= -3.91202300542815
-11.5129254649702 <= lnc_atp_c#244 <= -3.91202300542815
-11.5129254649702 <= lnc_cit_c#245 <= -3.91202300542815
-11.5129254649702 <= lnc_co2_c#246 <= -3.91202300542815
-11.5129254649702 <= lnc_co2_e#247 <= -3.91202300542815
-11.5129254649702 <= lnc_coa_c#248 <= -3.91202300542815
-11.5129254649702 <= lnc_dhap_c#249 <= -3.91202300542815
-11.5129254649702 <= lnc_e4p_c#250 <= -3.91202300542815
-11.5129254649702 <= lnc_etoh_c#251 <= -3.91202300542815
-11.5129254649702 <= lnc_etoh_e#252 <= -3.91202300542815
-11.5129254649702 <= lnc_f6p_c#253 <= -3.91202300542815
-11.5129254649702 <= lnc_fdp_c#254 <= -3.91202300542815
-11.5129254649702 <= lnc_for_c#255 <= -3.91202300542815
-11.5129254649702 <= lnc_for_e#256 <= -3.91202300542815
-11.5129254649702 <= lnc_fru_e#257 <= -3.91202300542815
-11.5129254649702 <= lnc_fum_c#258 <= -3.91202300542815
-11.5129254649702 <= lnc_fum_e#259 <= -3.91202300542815
-11.5129254649702 <= lnc_g3p_c#260 <= -3.91202300542815
-11.5129254649702 <= lnc_g6p_c#261 <= -3.91202300542815
-1.56574136370723 <= dG_err_glc__D_e#262 <= 1.56574136370723
-2.66718638217477 <= dG_err_gln__L_c#263 <= 2.66718638217477
-2.66718638217477 <= dG_err_gln__L_e#264 <= 2.66718638217477
-2.18985333115078 <= dG_err_glu__L_c#265 <= 2.18985333115078
-2.18985333115078 <= dG_err_glu__L_e#266 <= 2.18985333115078
-2.33208841551842 <= dG_err_glx_c#267 <= 2.33208841551842
-1.70760216396246 <= dG_err_h2o_c#268 <= 1.70760216396246
-1.70760216396246 <= dG_err_h2o_e#269 <= 1.70760216396246
      dG_err_h_c#270 = 0
      dG_err_h_e#271 = 0
-2.09939473274859 <= dG_err_icit_c#272 <= 2.09939473274859
-2.93059906776013 <= dG_err_lac__D_c#273 <= 2.93059906776013
-2.93059906776013 <= dG_err_lac__D_e#274 <= 2.93059906776013
-1.94756943070797 <= dG_err_mal__L_c#275 <= 1.94756943070797
-1.94756943070797 <= dG_err_mal__L_e#276 <= 1.94756943070797
-5.07102950069713 <= dG_err_nad_c#277 <= 5.07102950069713
-5.07228918737574 <= dG_err_nadh_c#278 <= 5.07228918737574
-5.12302347253818 <= dG_err_nadp_c#279 <= 5.12302347253818
-5.12194663053225 <= dG_err_nadph_c#280 <= 5.12194663053225
-2.03595671876732 <= dG_err_nh4_c#281 <= 2.03595671876732
-2.02737482279183 <= dG_err_13dpg_c#282 <= 2.02737482279183
-2.03595671876732 <= dG_err_nh4_e#283 <= 2.03595671876732
-3.35605045954032 <= dG_err_o2_c#284 <= 3.35605045954032
-1.75744556214876 <= dG_err_2pg_c#285 <= 1.75744556214876
-3.35605045954032 <= dG_err_o2_e#286 <= 3.35605045954032
-1.7187756477977 <= dG_err_3pg_c#287 <= 1.7187756477977
-2.01194071866202 <= dG_err_oaa_c#288 <= 2.01194071866202
-1.7764908130702 <= dG_err_pep_c#289 <= 1.7764908130702
-2.15149902927523 <= dG_err_6pgc_c#290 <= 2.15149902927523
-1.69596377974741 <= dG_err_pi_c#291 <= 1.69596377974741
-2.3753570326077 <= dG_err_6pgl_c#292 <= 2.3753570326077
-1.69596377974741 <= dG_err_pi_e#293 <= 1.69596377974741
-1.94388952649076 <= dG_err_ac_c#294 <= 1.94388952649076
-1.7070517707298 <= dG_err_pyr_c#295 <= 1.7070517707298
-1.7070517707298 <= dG_err_pyr_e#296 <= 1.7070517707298
-9.35163520981352 <= dG_err_q8_c#297 <= 9.35163520981352
-9.31857408402917 <= dG_err_q8h2_c#298 <= 9.31857408402917
-2.28408808818179 <= dG_err_r5p_c#299 <= 2.28408808818179
-1.9269434183073 <= dG_err_ru5p__D_c#300 <= 1.9269434183073
-1.94388952649076 <= dG_err_ac_e#301 <= 1.94388952649076
-2.16958210539232 <= dG_err_acald_c#302 <= 2.16958210539232
-2.6881202169194 <= dG_err_s7p_c#303 <= 2.6881202169194
-2.16958210539232 <= dG_err_acald_e#304 <= 2.16958210539232
-5.48338162621174 <= dG_err_accoa_c#305 <= 5.48338162621174
-2.66644244297168 <= dG_err_succ_c#306 <= 2.66644244297168
-2.66644244297168 <= dG_err_succ_e#307 <= 2.66644244297168
-5.55560169406521 <= dG_err_succoa_c#308 <= 5.55560169406521
-2.43199413123243 <= dG_err_acon_C_c#309 <= 2.43199413123243
-2.160798181301 <= dG_err_xu5p__D_c#310 <= 2.160798181301
-2.05626203594538 <= dG_err_actp_c#311 <= 2.05626203594538
-2.17002266209277 <= dG_err_adp_c#312 <= 2.17002266209277
-2.23896462728279 <= dG_err_akg_c#313 <= 2.23896462728279
-2.23896462728279 <= dG_err_akg_e#314 <= 2.23896462728279
-2.1444649259964 <= dG_err_amp_c#315 <= 2.1444649259964
-2.39247329428105 <= dG_err_atp_c#316 <= 2.39247329428105
-2.09442517972286 <= dG_err_cit_c#317 <= 2.09442517972286
-3.35605045953849 <= dG_err_co2_c#318 <= 3.35605045953849
-3.35605045953849 <= dG_err_co2_e#319 <= 3.35605045953849
-5.4975078304955 <= dG_err_coa_c#320 <= 5.4975078304955
-1.46858899178543 <= dG_err_dhap_c#321 <= 1.46858899178543
-2.57354281877347 <= dG_err_e4p_c#322 <= 2.57354281877347
-2.32577961931657 <= dG_err_etoh_c#323 <= 2.32577961931657
-2.32577961931657 <= dG_err_etoh_e#324 <= 2.32577961931657
-1.60148353447451 <= dG_err_f6p_c#325 <= 1.60148353447451
-2.0382757719938 <= dG_err_fdp_c#326 <= 2.0382757719938
-2.28357709516953 <= dG_err_for_c#327 <= 2.28357709516953
-2.28357709516953 <= dG_err_for_e#328 <= 2.28357709516953
-1.601103762682 <= dG_err_fru_e#329 <= 1.601103762682
-2.13451774607038 <= dG_err_fum_c#330 <= 2.13451774607038
-2.13451774607038 <= dG_err_fum_e#331 <= 2.13451774607038
-1.60021355329842 <= dG_err_g3p_c#332 <= 1.60021355329842
-1.57477704800091 <= dG_err_g6p_c#333 <= 1.57477704800091
-1000000 <= dG_PFK#334 <= 100000
-1000000 <= dG_PFK_reverse_d24a6#335 <= 100000
 0 <= indicator_PFK#336 <= 1
 0 <= indicator_PFK_reverse_d24a6#337 <= 1
-1000000 <= dG_PFL#338 <= 100000
-1000000 <= dG_PFL_reverse_af9ec#339 <= 100000
 0 <= indicator_PFL#340 <= 1
 0 <= indicator_PFL_reverse_af9ec#341 <= 1
-1000000 <= dG_PGI#342 <= 100000
-1000000 <= dG_PGI_reverse_27efc#343 <= 100000
 0 <= indicator_PGI#344 <= 1
 0 <= indicator_PGI_reverse_27efc#345 <= 1
-1000000 <= dG_PGK#346 <= 100000
-1000000 <= dG_PGK_reverse_02696#347 <= 100000
 0 <= indicator_PGK#348 <= 1
 0 <= indicator_PGK_reverse_02696#349 <= 1
-1000000 <= dG_PGL#350 <= 100000
-1000000 <= dG_PGL_reverse_2bb6b#351 <= 100000
 0 <= indicator_PGL#352 <= 1
 0 <= indicator_PGL_reverse_2bb6b#353 <= 1
-1000000 <= dG_ACALD#354 <= 100000
-1000000 <= dG_ACALD_reverse_fda2b#355 <= 100000
 0 <= indicator_ACALD#356 <= 1
 0 <= indicator_ACALD_reverse_fda2b#357 <= 1
-1000000 <= dG_AKGt2r#358 <= 100000
-1000000 <= dG_AKGt2r_reverse_5d500#359 <= 100000
 0 <= indicator_AKGt2r#360 <= 1
 0 <= indicator_AKGt2r_reverse_5d500#361 <= 1
-1000000 <= dG_PGM#362 <= 100000
-1000000 <= dG_PGM_reverse_fc9af#363 <= 100000
 0 <= indicator_PGM#364 <= 1
 0 <= indicator_PGM_reverse_fc9af#365 <= 1
-1000000 <= dG_PIt2r#366 <= 100000
-1000000 <= dG_PIt2r_reverse_1cd61#367 <= 100000
 0 <= indicator_PIt2r#368 <= 1
 0 <= indicator_PIt2r_reverse_1cd61#369 <= 1
-1000000 <= dG_ALCD2x#370 <= 100000
-1000000 <= dG_ALCD2x_reverse_5d107#371 <= 100000
 0 <= indicator_ALCD2x#372 <= 1
 0 <= indicator_ALCD2x_reverse_5d107#373 <= 1
-1000000 <= dG_ACALDt#374 <= 100000
-1000000 <= dG_ACALDt_reverse_858fa#375 <= 100000
 0 <= indicator_ACALDt#376 <= 1
 0 <= indicator_ACALDt_reverse_858fa#377 <= 1
-1000000 <= dG_ACKr#378 <= 100000
-1000000 <= dG_ACKr_reverse_b49c0#379 <= 100000
 0 <= indicator_ACKr#380 <= 1
 0 <= indicator_ACKr_reverse_b49c0#381 <= 1
-1000000 <= dG_PPC#382 <= 100000
-1000000 <= dG_PPC_reverse_e854a#383 <= 100000
 0 <= indicator_PPC#384 <= 1
 0 <= indicator_PPC_reverse_e854a#385 <= 1
-1000000 <= dG_ACONTa#386 <= 100000
-1000000 <= dG_ACONTa_reverse_cad6d#387 <= 100000
 0 <= indicator_ACONTa#388 <= 1
 0 <= indicator_ACONTa_reverse_cad6d#389 <= 1
-1000000 <= dG_ACONTb#390 <= 100000
-1000000 <= dG_ACONTb_reverse_e198a#391 <= 100000
 0 <= indicator_ACONTb#392 <= 1
 0 <= indicator_ACONTb_reverse_e198a#393 <= 1
-1000000 <= dG_ATPM#394 <= 100000
-1000000 <= dG_ATPM_reverse_5b752#395 <= 100000
 0 <= indicator_ATPM#396 <= 1
 0 <= indicator_ATPM_reverse_5b752#397 <= 1
-1000000 <= dG_PPCK#398 <= 100000
-1000000 <= dG_PPCK_reverse_2557d#399 <= 100000
 0 <= indicator_PPCK#400 <= 1
 0 <= indicator_PPCK_reverse_2557d#401 <= 1
-1000000 <= dG_ACt2r#402 <= 100000
-1000000 <= dG_ACt2r_reverse_64e05#403 <= 100000
 0 <= indicator_ACt2r#404 <= 1
 0 <= indicator_ACt2r_reverse_64e05#405 <= 1
-1000000 <= dG_PPS#406 <= 100000
-1000000 <= dG_PPS_reverse_1c319#407 <= 100000
 0 <= indicator_PPS#408 <= 1
 0 <= indicator_PPS_reverse_1c319#409 <= 1
-1000000 <= dG_ADK1#410 <= 100000
-1000000 <= dG_ADK1_reverse_a6f90#411 <= 100000
 0 <= indicator_ADK1#412 <= 1
 0 <= indicator_ADK1_reverse_a6f90#413 <= 1
-1000000 <= dG_AKGDH#414 <= 100000
-1000000 <= dG_AKGDH_reverse_08bdc#415 <= 100000
 0 <= indicator_AKGDH#416 <= 1
 0 <= indicator_AKGDH_reverse_08bdc#417 <= 1
-1000000 <= dG_ATPS4r#418 <= 100000
-1000000 <= dG_ATPS4r_reverse_64306#419 <= 100000
 0 <= indicator_ATPS4r#420 <= 1
 0 <= indicator_ATPS4r_reverse_64306#421 <= 1
-1000000 <= dG_PTAr#422 <= 100000
-1000000 <= dG_PTAr_reverse_fce15#423 <= 100000
 0 <= indicator_PTAr#424 <= 1
 0 <= indicator_PTAr_reverse_fce15#425 <= 1
-1000000 <= dG_PYK#426 <= 100000
-1000000 <= dG_PYK_reverse_bc8ff#427 <= 100000
 0 <= indicator_PYK#428 <= 1
 0 <= indicator_PYK_reverse_bc8ff#429 <= 1
-1000000 <= dG_PYRt2#430 <= 100000
-1000000 <= dG_PYRt2_reverse_a6e24#431 <= 100000
 0 <= indicator_PYRt2#432 <= 1
 0 <= indicator_PYRt2_reverse_a6e24#433 <= 1
-1000000 <= dG_CO2t#434 <= 100000
-1000000 <= dG_CO2t_reverse_7c42f#435 <= 100000
 0 <= indicator_CO2t#436 <= 1
 0 <= indicator_CO2t_reverse_7c42f#437 <= 1
-1000000 <= dG_RPE#438 <= 100000
-1000000 <= dG_RPE_reverse_a1b04#439 <= 100000
 0 <= indicator_RPE#440 <= 1
 0 <= indicator_RPE_reverse_a1b04#441 <= 1
-1000000 <= dG_CS#442 <= 100000
-1000000 <= dG_CS_reverse_8d7e9#443 <= 100000
 0 <= indicator_CS#444 <= 1
 0 <= indicator_CS_reverse_8d7e9#445 <= 1
-1000000 <= dG_RPI#446 <= 100000
-1000000 <= dG_RPI_reverse_853a1#447 <= 100000
 0 <= indicator_RPI#448 <= 1
 0 <= indicator_RPI_reverse_853a1#449 <= 1
-1000000 <= dG_SUCCt2_2#450 <= 100000
-1000000 <= dG_SUCCt2_2_reverse_206ca#451 <= 100000
 0 <= indicator_SUCCt2_2#452 <= 1
 0 <= indicator_SUCCt2_2_reverse_206ca#453 <= 1
-1000000 <= dG_CYTBD#454 <= 100000
-1000000 <= dG_CYTBD_reverse_611ba#455 <= 100000
 0 <= indicator_CYTBD#456 <= 1
 0 <= indicator_CYTBD_reverse_611ba#457 <= 1
-1000000 <= dG_D_LACt2#458 <= 100000
-1000000 <= dG_D_LACt2_reverse_d0dc8#459 <= 100000
 0 <= indicator_D_LACt2#460 <= 1
 0 <= indicator_D_LACt2_reverse_d0dc8#461 <= 1
-1000000 <= dG_ENO#462 <= 100000
-1000000 <= dG_ENO_reverse_40eea#463 <= 100000
 0 <= indicator_ENO#464 <= 1
 0 <= indicator_ENO_reverse_40eea#465 <= 1
-1000000 <= dG_SUCCt3#466 <= 100000
-1000000 <= dG_SUCCt3_reverse_7d82a#467 <= 100000
 0 <= indicator_SUCCt3#468 <= 1
 0 <= indicator_SUCCt3_reverse_7d82a#469 <= 1
-1000000 <= dG_ETOHt2r#470 <= 100000
-1000000 <= dG_ETOHt2r_reverse_b3d52#471 <= 100000
 0 <= indicator_ETOHt2r#472 <= 1
 0 <= indicator_ETOHt2r_reverse_b3d52#473 <= 1
-1000000 <= dG_SUCDi#474 <= 100000
-1000000 <= dG_SUCDi_reverse_480f4#475 <= 100000
 0 <= indicator_SUCDi#476 <= 1
 0 <= indicator_SUCDi_reverse_480f4#477 <= 1
-1000000 <= dG_SUCOAS#478 <= 100000
-1000000 <= dG_SUCOAS_reverse_22958#479 <= 100000
 0 <= indicator_SUCOAS#480 <= 1
 0 <= indicator_SUCOAS_reverse_22958#481 <= 1
-1000000 <= dG_TALA#482 <= 100000
-1000000 <= dG_TALA_reverse_adfda#483 <= 100000
 0 <= indicator_TALA#484 <= 1
 0 <= indicator_TALA_reverse_adfda#485 <= 1
-1000000 <= dG_THD2#486 <= 100000
-1000000 <= dG_THD2_reverse_f65dd#487 <= 100000
 0 <= indicator_THD2#488 <= 1
 0 <= indicator_THD2_reverse_f65dd#489 <= 1
-1000000 <= dG_TKT1#490 <= 100000
-1000000 <= dG_TKT1_reverse_a1021#491 <= 100000
 0 <= indicator_TKT1#492 <= 1
 0 <= indicator_TKT1_reverse_a1021#493 <= 1
-1000000 <= dG_TKT2#494 <= 100000
-1000000 <= dG_TKT2_reverse_7ebc7#495 <= 100000
 0 <= indicator_TKT2#496 <= 1
 0 <= indicator_TKT2_reverse_7ebc7#497 <= 1
-1000000 <= dG_TPI#498 <= 100000
-1000000 <= dG_TPI_reverse_c2c3b#499 <= 100000
 0 <= indicator_TPI#500 <= 1
 0 <= indicator_TPI_reverse_c2c3b#501 <= 1
-1000000 <= dG_FBA#502 <= 100000
-1000000 <= dG_FBA_reverse_84806#503 <= 100000
 0 <= indicator_FBA#504 <= 1
 0 <= indicator_FBA_reverse_84806#505 <= 1
-1000000 <= dG_FBP#506 <= 100000
-1000000 <= dG_FBP_reverse_bf2c9#507 <= 100000
 0 <= indicator_FBP#508 <= 1
 0 <= indicator_FBP_reverse_bf2c9#509 <= 1
-1000000 <= dG_FORt2#510 <= 100000
-1000000 <= dG_FORt2_reverse_89839#511 <= 100000
 0 <= indicator_FORt2#512 <= 1
 0 <= indicator_FORt2_reverse_89839#513 <= 1
-1000000 <= dG_FORt#514 <= 100000
-1000000 <= dG_FORt_reverse_40f9f#515 <= 100000
 0 <= indicator_FORt#516 <= 1
 0 <= indicator_FORt_reverse_40f9f#517 <= 1
-1000000 <= dG_FRD7#518 <= 100000
-1000000 <= dG_FRD7_reverse_e0172#519 <= 100000
 0 <= indicator_FRD7#520 <= 1
 0 <= indicator_FRD7_reverse_e0172#521 <= 1
-1000000 <= dG_FRUpts2#522 <= 100000
-1000000 <= dG_FRUpts2_reverse_58d3a#523 <= 100000
 0 <= indicator_FRUpts2#524 <= 1
 0 <= indicator_FRUpts2_reverse_58d3a#525 <= 1
-1000000 <= dG_FUM#526 <= 100000
-1000000 <= dG_FUM_reverse_d3642#527 <= 100000
 0 <= indicator_FUM#528 <= 1
 0 <= indicator_FUM_reverse_d3642#529 <= 1
-1000000 <= dG_FUMt2_2#530 <= 100000
-1000000 <= dG_FUMt2_2_reverse_7b4b6#531 <= 100000
 0 <= indicator_FUMt2_2#532 <= 1
 0 <= indicator_FUMt2_2_reverse_7b4b6#533 <= 1
-1000000 <= dG_G6PDH2r#534 <= 100000
-1000000 <= dG_G6PDH2r_reverse_19ddf#535 <= 100000
 0 <= indicator_G6PDH2r#536 <= 1
 0 <= indicator_G6PDH2r_reverse_19ddf#537 <= 1
-1000000 <= dG_GAPD#538 <= 100000
-1000000 <= dG_GAPD_reverse_459c1#539 <= 100000
 0 <= indicator_GAPD#540 <= 1
 0 <= indicator_GAPD_reverse_459c1#541 <= 1
-1000000 <= dG_GLCpts#542 <= 100000
-1000000 <= dG_GLCpts_reverse_a52ae#543 <= 100000
 0 <= indicator_GLCpts#544 <= 1
 0 <= indicator_GLCpts_reverse_a52ae#545 <= 1
-1000000 <= dG_GLNS#546 <= 100000
-1000000 <= dG_GLNS_reverse_59581#547 <= 100000
 0 <= indicator_GLNS#548 <= 1
 0 <= indicator_GLNS_reverse_59581#549 <= 1
-1000000 <= dG_GLNabc#550 <= 100000
-1000000 <= dG_GLNabc_reverse_1d82a#551 <= 100000
 0 <= indicator_GLNabc#552 <= 1
 0 <= indicator_GLNabc_reverse_1d82a#553 <= 1
-1000000 <= dG_GLUDy#554 <= 100000
-1000000 <= dG_GLUDy_reverse_fa4e7#555 <= 100000
 0 <= indicator_GLUDy#556 <= 1
 0 <= indicator_GLUDy_reverse_fa4e7#557 <= 1
-1000000 <= dG_GLUN#558 <= 100000
-1000000 <= dG_GLUN_reverse_4ccdb#559 <= 100000
 0 <= indicator_GLUN#560 <= 1
 0 <= indicator_GLUN_reverse_4ccdb#561 <= 1
-1000000 <= dG_GLUSy#562 <= 100000
-1000000 <= dG_GLUSy_reverse_6a00f#563 <= 100000
 0 <= indicator_GLUSy#564 <= 1
 0 <= indicator_GLUSy_reverse_6a00f#565 <= 1
-1000000 <= dG_GLUt2r#566 <= 100000
-1000000 <= dG_GLUt2r_reverse_3e88a#567 <= 100000
 0 <= indicator_GLUt2r#568 <= 1
 0 <= indicator_GLUt2r_reverse_3e88a#569 <= 1
-1000000 <= dG_GND#570 <= 100000
-1000000 <= dG_GND_reverse_eec5c#571 <= 100000
 0 <= indicator_GND#572 <= 1
 0 <= indicator_GND_reverse_eec5c#573 <= 1
-1000000 <= dG_ICDHyr#574 <= 100000
-1000000 <= dG_ICDHyr_reverse_7f84b#575 <= 100000
 0 <= indicator_ICDHyr#576 <= 1
 0 <= indicator_ICDHyr_reverse_7f84b#577 <= 1
-1000000 <= dG_ICL#578 <= 100000
-1000000 <= dG_ICL_reverse_2f27e#579 <= 100000
 0 <= indicator_ICL#580 <= 1
 0 <= indicator_ICL_reverse_2f27e#581 <= 1
-1000000 <= dG_LDH_D#582 <= 100000
-1000000 <= dG_LDH_D_reverse_f8507#583 <= 100000
 0 <= indicator_LDH_D#584 <= 1
 0 <= indicator_LDH_D_reverse_f8507#585 <= 1
-1000000 <= dG_MALS#586 <= 100000
-1000000 <= dG_MALS_reverse_d7382#587 <= 100000
 0 <= indicator_MALS#588 <= 1
 0 <= indicator_MALS_reverse_d7382#589 <= 1
-1000000 <= dG_MALt2_2#590 <= 100000
-1000000 <= dG_MALt2_2_reverse_a635f#591 <= 100000
 0 <= indicator_MALt2_2#592 <= 1
 0 <= indicator_MALt2_2_reverse_a635f#593 <= 1
-1000000 <= dG_MDH#594 <= 100000
-1000000 <= dG_MDH_reverse_ee52c#595 <= 100000
 0 <= indicator_MDH#596 <= 1
 0 <= indicator_MDH_reverse_ee52c#597 <= 1
-1000000 <= dG_ME1#598 <= 100000
-1000000 <= dG_ME1_reverse_9736c#599 <= 100000
 0 <= indicator_ME1#600 <= 1
 0 <= indicator_ME1_reverse_9736c#601 <= 1
-1000000 <= dG_ME2#602 <= 100000
-1000000 <= dG_ME2_reverse_2b0a2#603 <= 100000
 0 <= indicator_ME2#604 <= 1
 0 <= indicator_ME2_reverse_2b0a2#605 <= 1
-1000000 <= dG_NADH16#606 <= 100000
-1000000 <= dG_NADH16_reverse_330a2#607 <= 100000
 0 <= indicator_NADH16#608 <= 1
 0 <= indicator_NADH16_reverse_330a2#609 <= 1
-1000000 <= dG_NADTRHD#610 <= 100000
-1000000 <= dG_NADTRHD_reverse_49725#611 <= 100000
 0 <= indicator_NADTRHD#612 <= 1
 0 <= indicator_NADTRHD_reverse_49725#613 <= 1
-1000000 <= dG_NH4t#614 <= 100000
-1000000 <= dG_NH4t_reverse_551ee#615 <= 100000
 0 <= indicator_NH4t#616 <= 1
 0 <= indicator_NH4t_reverse_551ee#617 <= 1
-1000000 <= dG_PDH#618 <= 100000
-1000000 <= dG_PDH_reverse_ca160#619 <= 100000
 0 <= indicator_PDH#620 <= 1
 0 <= indicator_PDH_reverse_ca160#621 <= 1
Binaries
 indicator_PFK#336  indicator_PFK_reverse_d24a6#337  indicator_PFL#340 
 indicator_PFL_reverse_af9ec#341  indicator_PGI#344 
 indicator_PGI_reverse_27efc#345  indicator_PGK#348 
 indicator_PGK_reverse_02696#349  indicator_PGL#352 
 indicator_PGL_reverse_2bb6b#353  indicator_ACALD#356 
 indicator_ACALD_reverse_fda2b#357  indicator_AKGt2r#360 
 indicator_AKGt2r_reverse_5d500#361  indicator_PGM#364 
 indicator_PGM_reverse_fc9af#365  indicator_PIt2r#368 
 indicator_PIt2r_reverse_1cd61#369  indicator_ALCD2x#372 
 indicator_ALCD2x_reverse_5d107#373  indicator_ACALDt#376 
 indicator_ACALDt_reverse_858fa#377  indicator_ACKr#380 
 indicator_ACKr_reverse_b49c0#381  indicator_PPC#384 
 indicator_PPC_reverse_e854a#385  indicator_ACONTa#388 
 indicator_ACONTa_reverse_cad6d#389  indicator_ACONTb#392 
 indicator_ACONTb_reverse_e198a#393  indicator_ATPM#396 
 indicator_ATPM_reverse_5b752#397  indicator_PPCK#400 
 indicator_PPCK_reverse_2557d#401  indicator_ACt2r#404 
 indicator_ACt2r_reverse_64e05#405  indicator_PPS#408 
 indicator_PPS_reverse_1c319#409  indicator_ADK1#412 
 indicator_ADK1_reverse_a6f90#413  indicator_AKGDH#416 
 indicator_AKGDH_reverse_08bdc#417  indicator_ATPS4r#420 
 indicator_ATPS4r_reverse_64306#421  indicator_PTAr#424 
 indicator_PTAr_reverse_fce15#425  indicator_PYK#428 
 indicator_PYK_reverse_bc8ff#429  indicator_PYRt2#432 
 indicator_PYRt2_reverse_a6e24#433  indicator_CO2t#436 
 indicator_CO2t_reverse_7c42f#437  indicator_RPE#440 
 indicator_RPE_reverse_a1b04#441  indicator_CS#444 
 indicator_CS_reverse_8d7e9#445  indicator_RPI#448 
 indicator_RPI_reverse_853a1#449  indicator_SUCCt2_2#452 
 indicator_SUCCt2_2_reverse_206ca#453  indicator_CYTBD#456 
 indicator_CYTBD_reverse_611ba#457  indicator_D_LACt2#460 
 indicator_D_LACt2_reverse_d0dc8#461  indicator_ENO#464 
 indicator_ENO_reverse_40eea#465  indicator_SUCCt3#468 
 indicator_SUCCt3_reverse_7d82a#469  indicator_ETOHt2r#472 
 indicator_ETOHt2r_reverse_b3d52#473  indicator_SUCDi#476 
 indicator_SUCDi_reverse_480f4#477  indicator_SUCOAS#480 
 indicator_SUCOAS_reverse_22958#481  indicator_TALA#484 
 indicator_TALA_reverse_adfda#485  indicator_THD2#488 
 indicator_THD2_reverse_f65dd#489  indicator_TKT1#492 
 indicator_TKT1_reverse_a1021#493  indicator_TKT2#496 
 indicator_TKT2_reverse_7ebc7#497  indicator_TPI#500 
 indicator_TPI_reverse_c2c3b#501  indicator_FBA#504 
 indicator_FBA_reverse_84806#505  indicator_FBP#508 
 indicator_FBP_reverse_bf2c9#509  indicator_FORt2#512 
 indicator_FORt2_reverse_89839#513  indicator_FORt#516 
 indicator_FORt_reverse_40f9f#517  indicator_FRD7#520 
 indicator_FRD7_reverse_e0172#521  indicator_FRUpts2#524 
 indicator_FRUpts2_reverse_58d3a#525  indicator_FUM#528 
 indicator_FUM_reverse_d3642#529  indicator_FUMt2_2#532 
 indicator_FUMt2_2_reverse_7b4b6#533  indicator_G6PDH2r#536 
 indicator_G6PDH2r_reverse_19ddf#537  indicator_GAPD#540 
 indicator_GAPD_reverse_459c1#541  indicator_GLCpts#544 
 indicator_GLCpts_reverse_a52ae#545  indicator_GLNS#548 
 indicator_GLNS_reverse_59581#549  indicator_GLNabc#552 
 indicator_GLNabc_reverse_1d82a#553  indicator_GLUDy#556 
 indicator_GLUDy_reverse_fa4e7#557  indicator_GLUN#560 
 indicator_GLUN_reverse_4ccdb#561  indicator_GLUSy#564 
 indicator_GLUSy_reverse_6a00f#565  indicator_GLUt2r#568 
 indicator_GLUt2r_reverse_3e88a#569  indicator_GND#572 
 indicator_GND_reverse_eec5c#573  indicator_ICDHyr#576 
 indicator_ICDHyr_reverse_7f84b#577  indicator_ICL#580 
 indicator_ICL_reverse_2f27e#581  indicator_LDH_D#584 
 indicator_LDH_D_reverse_f8507#585  indicator_MALS#588 
 indicator_MALS_reverse_d7382#589  indicator_MALt2_2#592 
 indicator_MALt2_2_reverse_a635f#593  indicator_MDH#596 
 indicator_MDH_reverse_ee52c#597  indicator_ME1#600 
 indicator_ME1_reverse_9736c#601  indicator_ME2#604 
 indicator_ME2_reverse_2b0a2#605  indicator_NADH16#608 
 indicator_NADH16_reverse_330a2#609  indicator_NADTRHD#612 
 indicator_NADTRHD_reverse_49725#613  indicator_NH4t#616 
 indicator_NH4t_reverse_551ee#617  indicator_PDH#620 
 indicator_PDH_reverse_ca160#621 
End
