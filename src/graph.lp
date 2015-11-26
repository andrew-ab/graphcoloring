\ENCODING=ISO-8859-1
\Problem name: Instance of partitioned graph coloring.

Minimize
 obj: w_1 + w_2 + w_3 + w_4 + w_5
Subject To
 Blue:         x_11 + x_21 <= 1
 Red:          x_12 + x_22 <= 1
 AliceBlue:    x_13 + x_23 <= 1
 AntiqueWhite: x_14 + x_24 <= 1
 Aqua:         x_15 + x_25 <= 1
 Blue:         x_11 + x_41 <= 1
 Red:          x_12 + x_42 <= 1
 AliceBlue:    x_13 + x_43 <= 1
 AntiqueWhite: x_14 + x_44 <= 1
 Aqua:         x_15 + x_45 <= 1
 Blue:         x_11 + x_71 <= 1
 Red:          x_12 + x_72 <= 1
 AliceBlue:    x_13 + x_73 <= 1
 AntiqueWhite: x_14 + x_74 <= 1
 Aqua:         x_15 + x_75 <= 1
 Blue:         x_11 + x_91 <= 1
 Red:          x_12 + x_92 <= 1
 AliceBlue:    x_13 + x_93 <= 1
 AntiqueWhite: x_14 + x_94 <= 1
 Aqua:         x_15 + x_95 <= 1
 Blue:         x_21 + x_31 <= 1
 Red:          x_22 + x_32 <= 1
 AliceBlue:    x_23 + x_33 <= 1
 AntiqueWhite: x_24 + x_34 <= 1
 Aqua:         x_25 + x_35 <= 1
 Blue:         x_21 + x_61 <= 1
 Red:          x_22 + x_62 <= 1
 AliceBlue:    x_23 + x_63 <= 1
 AntiqueWhite: x_24 + x_64 <= 1
 Aqua:         x_25 + x_65 <= 1
 Blue:         x_21 + x_81 <= 1
 Red:          x_22 + x_82 <= 1
 AliceBlue:    x_23 + x_83 <= 1
 AntiqueWhite: x_24 + x_84 <= 1
 Aqua:         x_25 + x_85 <= 1
 Blue:         x_31 + x_51 <= 1
 Red:          x_32 + x_52 <= 1
 AliceBlue:    x_33 + x_53 <= 1
 AntiqueWhite: x_34 + x_54 <= 1
 Aqua:         x_35 + x_55 <= 1
 Blue:         x_31 + x_71 <= 1
 Red:          x_32 + x_72 <= 1
 AliceBlue:    x_33 + x_73 <= 1
 AntiqueWhite: x_34 + x_74 <= 1
 Aqua:         x_35 + x_75 <= 1
 Blue:         x_31 + x_101 <= 1
 Red:          x_32 + x_102 <= 1
 AliceBlue:    x_33 + x_103 <= 1
 AntiqueWhite: x_34 + x_104 <= 1
 Aqua:         x_35 + x_105 <= 1
 Blue:         x_41 + x_51 <= 1
 Red:          x_42 + x_52 <= 1
 AliceBlue:    x_43 + x_53 <= 1
 AntiqueWhite: x_44 + x_54 <= 1
 Aqua:         x_45 + x_55 <= 1
 Blue:         x_41 + x_61 <= 1
 Red:          x_42 + x_62 <= 1
 AliceBlue:    x_43 + x_63 <= 1
 AntiqueWhite: x_44 + x_64 <= 1
 Aqua:         x_45 + x_65 <= 1
 Blue:         x_41 + x_101 <= 1
 Red:          x_42 + x_102 <= 1
 AliceBlue:    x_43 + x_103 <= 1
 AntiqueWhite: x_44 + x_104 <= 1
 Aqua:         x_45 + x_105 <= 1
 Blue:         x_51 + x_81 <= 1
 Red:          x_52 + x_82 <= 1
 AliceBlue:    x_53 + x_83 <= 1
 AntiqueWhite: x_54 + x_84 <= 1
 Aqua:         x_55 + x_85 <= 1
 Blue:         x_51 + x_91 <= 1
 Red:          x_52 + x_92 <= 1
 AliceBlue:    x_53 + x_93 <= 1
 AntiqueWhite: x_54 + x_94 <= 1
 Aqua:         x_55 + x_95 <= 1
 Blue:         x_61 + x_111 <= 1
 Red:          x_62 + x_112 <= 1
 AliceBlue:    x_63 + x_113 <= 1
 AntiqueWhite: x_64 + x_114 <= 1
 Aqua:         x_65 + x_115 <= 1
 Blue:         x_71 + x_111 <= 1
 Red:          x_72 + x_112 <= 1
 AliceBlue:    x_73 + x_113 <= 1
 AntiqueWhite: x_74 + x_114 <= 1
 Aqua:         x_75 + x_115 <= 1
 Blue:         x_81 + x_111 <= 1
 Red:          x_82 + x_112 <= 1
 AliceBlue:    x_83 + x_113 <= 1
 AntiqueWhite: x_84 + x_114 <= 1
 Aqua:         x_85 + x_115 <= 1
 Blue:         x_91 + x_111 <= 1
 Red:          x_92 + x_112 <= 1
 AliceBlue:    x_93 + x_113 <= 1
 AntiqueWhite: x_94 + x_114 <= 1
 Aqua:         x_95 + x_115 <= 1
 Blue:         x_101 + x_111 <= 1
 Red:          x_102 + x_112 <= 1
 AliceBlue:    x_103 + x_113 <= 1
 AntiqueWhite: x_104 + x_114 <= 1
 Aqua:         x_105 + x_115 <= 1
 partition_1:  x_61 + x_62 + x_63 + x_64 + x_65 <= 1
 partition_2:  x_31 + x_32 + x_33 + x_34 + x_35 + x_111 + x_112 + x_113 + x_114
               + x_115 <= 1
 partition_3:  x_81 + x_82 + x_83 + x_84 + x_85 + x_101 + x_102 + x_103 + x_104
               + x_105 <= 1
 partition_4:  x_41 + x_42 + x_43 + x_44 + x_45 + x_51 + x_52 + x_53 + x_54
               + x_55 <= 1
 partition_5:  x_11 + x_12 + x_13 + x_14 + x_15 + x_21 + x_22 + x_23 + x_24
               + x_25 + x_71 + x_72 + x_73 + x_74 + x_75 + x_91 + x_92 + x_93
               + x_94 + x_95 <= 1
 color_res:    - w_1 + x_11 <= 0
 color_res:    - w_2 + x_12 <= 0
 color_res:    - w_3 + x_13 <= 0
 color_res:    - w_4 + x_14 <= 0
 color_res:    - w_5 + x_15 <= 0
 color_res:    - w_1 + x_21 <= 0
 color_res:    - w_2 + x_22 <= 0
 color_res:    - w_3 + x_23 <= 0
 color_res:    - w_4 + x_24 <= 0
 color_res:    - w_5 + x_25 <= 0
 color_res:    - w_1 + x_31 <= 0
 color_res:    - w_2 + x_32 <= 0
 color_res:    - w_3 + x_33 <= 0
 color_res:    - w_4 + x_34 <= 0
 color_res:    - w_5 + x_35 <= 0
 color_res:    - w_1 + x_41 <= 0
 color_res:    - w_2 + x_42 <= 0
 color_res:    - w_3 + x_43 <= 0
 color_res:    - w_4 + x_44 <= 0
 color_res:    - w_5 + x_45 <= 0
 color_res:    - w_1 + x_51 <= 0
 color_res:    - w_2 + x_52 <= 0
 color_res:    - w_3 + x_53 <= 0
 color_res:    - w_4 + x_54 <= 0
 color_res:    - w_5 + x_55 <= 0
 color_res:    - w_1 + x_61 <= 0
 color_res:    - w_2 + x_62 <= 0
 color_res:    - w_3 + x_63 <= 0
 color_res:    - w_4 + x_64 <= 0
 color_res:    - w_5 + x_65 <= 0
 color_res:    - w_1 + x_71 <= 0
 color_res:    - w_2 + x_72 <= 0
 color_res:    - w_3 + x_73 <= 0
 color_res:    - w_4 + x_74 <= 0
 color_res:    - w_5 + x_75 <= 0
 color_res:    - w_1 + x_81 <= 0
 color_res:    - w_2 + x_82 <= 0
 color_res:    - w_3 + x_83 <= 0
 color_res:    - w_4 + x_84 <= 0
 color_res:    - w_5 + x_85 <= 0
 color_res:    - w_1 + x_91 <= 0
 color_res:    - w_2 + x_92 <= 0
 color_res:    - w_3 + x_93 <= 0
 color_res:    - w_4 + x_94 <= 0
 color_res:    - w_5 + x_95 <= 0
 color_res:    - w_1 + x_101 <= 0
 color_res:    - w_2 + x_102 <= 0
 color_res:    - w_3 + x_103 <= 0
 color_res:    - w_4 + x_104 <= 0
 color_res:    - w_5 + x_105 <= 0
 color_res:    - w_1 + x_111 <= 0
 color_res:    - w_2 + x_112 <= 0
 color_res:    - w_3 + x_113 <= 0
 color_res:    - w_4 + x_114 <= 0
 color_res:    - w_5 + x_115 <= 0
Bounds
 0 <= w_1 <= 1
 0 <= w_2 <= 1
 0 <= w_3 <= 1
 0 <= w_4 <= 1
 0 <= w_5 <= 1
 0 <= x_11 <= 1
 0 <= x_12 <= 1
 0 <= x_13 <= 1
 0 <= x_14 <= 1
 0 <= x_15 <= 1
 0 <= x_21 <= 1
 0 <= x_22 <= 1
 0 <= x_23 <= 1
 0 <= x_24 <= 1
 0 <= x_25 <= 1
 0 <= x_31 <= 1
 0 <= x_32 <= 1
 0 <= x_33 <= 1
 0 <= x_34 <= 1
 0 <= x_35 <= 1
 0 <= x_41 <= 1
 0 <= x_42 <= 1
 0 <= x_43 <= 1
 0 <= x_44 <= 1
 0 <= x_45 <= 1
 0 <= x_51 <= 1
 0 <= x_52 <= 1
 0 <= x_53 <= 1
 0 <= x_54 <= 1
 0 <= x_55 <= 1
 0 <= x_61 <= 1
 0 <= x_62 <= 1
 0 <= x_63 <= 1
 0 <= x_64 <= 1
 0 <= x_65 <= 1
 0 <= x_71 <= 1
 0 <= x_72 <= 1
 0 <= x_73 <= 1
 0 <= x_74 <= 1
 0 <= x_75 <= 1
 0 <= x_81 <= 1
 0 <= x_82 <= 1
 0 <= x_83 <= 1
 0 <= x_84 <= 1
 0 <= x_85 <= 1
 0 <= x_91 <= 1
 0 <= x_92 <= 1
 0 <= x_93 <= 1
 0 <= x_94 <= 1
 0 <= x_95 <= 1
 0 <= x_101 <= 1
 0 <= x_102 <= 1
 0 <= x_103 <= 1
 0 <= x_104 <= 1
 0 <= x_105 <= 1
 0 <= x_111 <= 1
 0 <= x_112 <= 1
 0 <= x_113 <= 1
 0 <= x_114 <= 1
 0 <= x_115 <= 1
Binaries
 w_1  w_2  w_3  w_4  w_5  x_11  x_12  x_13  x_14  x_15  x_21  x_22  x_23  x_24 
 x_25  x_31  x_32  x_33  x_34  x_35  x_41  x_42  x_43  x_44  x_45  x_51  x_52 
 x_53  x_54  x_55  x_61  x_62  x_63  x_64  x_65  x_71  x_72  x_73  x_74  x_75 
 x_81  x_82  x_83  x_84  x_85  x_91  x_92  x_93  x_94  x_95  x_101  x_102 
 x_103  x_104  x_105  x_111  x_112  x_113  x_114  x_115 
End
