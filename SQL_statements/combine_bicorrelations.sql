/** Obtain Bicorrelations with Gene Symbol **/

select cor.probeid1, cor.probeid2, cor.value as bicor, p.value as pvalue 
	into [HMDP].[dbo].[bicor_all_reactive_liver_male] 
	from [HMDP].[dbo].[bicor_expression_differences_liver_male] as cor
  inner join
  [HMDP].[dbo].[bicor_pval_expression_differences_liver_male] as p
  on cor.probeid1=p.probeid1 and cor.probeid2 = p.probeid2
;


select cor.probeid1, affy.gene_symbol as gene_symbol1, affy.gene_title as gene_title1, affy.chr as chr1, affy.start_bp as start_bp1, affy.end_bp as end_bp1, affy.strand as strand1,
	cor.probeid2, affy2.gene_symbol as gene_symbol2, affy2.gene_title as gene_title2, affy2.chr as chr2, affy2.start_bp as start_bp2, affy2.end_bp as end_bp2, affy2.strand as strand2,
	cor.bicor, cor.pvalue
	into [HMDP].[dbo].[bicor_combined_reactive_liver_male]
  from bicor_all_reactive_liver_male as cor
  left join 
  (select probesetID, gene_symbol, gene_title, chr, start_bp, end_bp, strand from [HMDP].[dbo].[affy_microarray_annotation]) as affy
  on cor.probeid1 = affy.probesetID
  left join
  (select probesetID, gene_symbol, gene_title, chr, start_bp, end_bp, strand from [HMDP].[dbo].[affy_microarray_annotation]) as affy2
  on cor.probeid2 = affy2.probesetID
 ;


insert into [HMDP].[dbo].[bicorrelations_reactive_liver_male]
 select * from [HMDP].[dbo].[bicor_combined_reactive_liver_male]
;
