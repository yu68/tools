
#grep -v "random" ~/MNase-seq/May_2013_MNase-seq_data/expression/result/pooled/Pn_2_GCCAAT_pairend_rmdup.bed | head -n +1000 > example/bed1.bed
#grep -v "random" ~/MNase-seq/May_2013_MNase-seq_data/expression/result/pooled/Pn_3_CTTGTA_pairend_rmdup.bed | head -n +1000 > example/bed2.bed
#grep -v "random" ~/MNase-seq/May_2013_MNase-seq_data/expression/result/pooled/PnLuci_ACAGTG_pairend_rmdup.bed | head -n +1000 > example/bed3.bed
#grep -v "random" ~/MNase-seq/May_2013_MNase-seq_data/expression/result/pooled/PnpSuper_GTGAAA_pairend_rmdup.bed | head -n +1000 > example/bed4.bed

Rscript Venn_diagram_BED.R example/bed1.bed example/bed2.bed example/bed4.bed -w -m 5
