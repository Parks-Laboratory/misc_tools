/** Aggregate RMA Expression by Strain and Add Gene Annotations **/

/*Note: to change diet or tissue do a replace all*/

create view avg_expression_by_strain_liver_highfat_female as
/*inner join aggregated expression to affy annotations*/
select gene_symbol, chr, agg_expression.probe_id, strain, mean_rma_expression 
	from (
		/*inner join microarray to annotation AND aggregate expression by probe and strain*/
		select annotate.strain, microarray.probe_id, AVG(microarray.rma_expression) as mean_rma_expression 
		from dbo.microarray_highfat_liver as microarray
		inner join
		(select * from dbo.annotation_highfat where sex = 'Female') as annotate
		on annotate.mouse_id = microarray.mouse_id
		group by strain, probe_id
	) as agg_expression
	left join
	(select probesetID, gene_symbol, chr from dbo.affy_microarray_annotation) as affy_annotate
	on agg_expression.probe_id = affy_annotate.probesetID
;
