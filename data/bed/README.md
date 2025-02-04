# Create BED files for coverage analysis outside of the VNTR region

# for the TWIST v2.0 exome
bedtools intersect -a data/bed/_muc1_regions/muc1_region_hg38.bed -b data/bed/_original_files/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed | bedtools subtract -a - -b data/bed/_muc1_regions/muc1_vntr_hg38.bed > data/bed/twist_v2_muc1_no_vntr_targets.bed

# for the IDT xGen Research Panel V1
bedtools intersect -a data/bed/_muc1_regions/muc1_region_hg38.bed -b data/bed/_original_files/IDT_xGen_Research_Targets_V1_hg38.bed | bedtools subtract -a - -b data/bed/muc1_vntr_hg38.bed > data/bed/_muc1_regions/idt_v1_muc1_no_vntr_targets.bed

# for the Agilent SureSelectXT Human All Exon V8
bedtools intersect -a data/bed/_muc1_regions/muc1_region_hg38.bed -b data/bed/_original_files/SureSelectXT_Human_All_Exon_V8_S33266436_Regions_hg38.bed | bedtools subtract -a - -b data/bed/_muc1_regions/muc1_vntr_hg38.bed > data/bed/agilent_v8_muc1_no_vntr_targets.bed

## TWIST target files
Downloaded from [https://www.twistbioscience.com/ on 2023-10-18](https://www.twistbioscience.com/resources/data-files/twist-exome-20-bed-files) on 2025-02-04

## checksums:
29b4f306472b88604973f0e3c798cfd1  Twist_Exome_v2.0.2_Targets_hg19.bed
c3a7cea67f992e0412db4b596730d276  Twist_Exome_v2.0.2_Targets_hg38.bed

## IDT target files

Files have been downloaded from UCSC Table browser (https://genome.ucsc.edu/cgi-bin/hgTables>) using Galaxy on 2023-10-18.
(Mammal - Human -Mapping and Sequencing - Exome Probesets - IDT xGen V1.0 (xGen_Research_Targets))

## checksums:
139b7d428fee47a014205dd32937a5b0  IDT_xGen_Research_Targets_V1_hg38.bed
f3b253a0cf1fb39e91ce25b9fe5bc0f2  IDT_xGen_Research_Targets_V1_hg19.bed

## Agilent target files
downloaded from https://earray.chem.agilent.com/suredesign/search/entity.htm on 2025-02-04

## checksums:
a8f5165475442ae6165a434f0c37ce27  SureSelectXT_Human_All_Exon_V8_S33266436_Regions_hg19.bed
c592685888f54bae1cfdfc808f8c5d58  SureSelectXT_Human_All_Exon_V8_S33266436_Regions_hg38.bed
