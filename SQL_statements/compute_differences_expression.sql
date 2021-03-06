/** Compute Differences between RMA Expression Across Highfat & Chow **/

create view expression_differences_liver_male as
select highfat.gene_symbol, highfat.chr, highfat.probe_id, highfat.strain, (highfat.mean_rma_expression - chow.mean_rma_expression) as expression_difference
	from dbo.avg_expression_by_strain_liver_highfat_male as highfat
	inner join 
	dbo.avg_expression_by_strain_liver_chow_male as chow
	on highfat.strain = chow.strain and highfat.probe_id = chow.probe_id
;
