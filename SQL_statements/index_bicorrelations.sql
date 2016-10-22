/** Index on BicorrelationsGene Symbol **/

USE [HMDP]
GO
CREATE NONCLUSTERED INDEX [index_on_gene_symbol]
ON [dbo].[bicorrelations_adipose_highfat_female] ([gene_symbol1])
INCLUDE ([probeid1],[gene_title1],[chr1],[start_bp1],[end_bp1],[strand1],[probeid2],[gene_symbol2],[gene_title2],[chr2],[start_bp2],[end_bp2],[strand2],[bicor],[pvalue])
GO

