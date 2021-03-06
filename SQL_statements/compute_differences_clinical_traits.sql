/****** Compute Differences Between HF and Chow Clinical Traits ******/

/* initialize local vars */
DECLARE @sql AS NVARCHAR(2000);
DECLARE @col AS NVARCHAR(2000);

/* using common columns, write part of sql statements regarding differences */
SELECT @col = ISNULL(@col + ', ', '') + test
FROM (SELECT concat('(hf.', highfat, ' - ', 'c.', chow, ') as ', common_name) as test FROM [HMDP].[dbo].[clinical_traits_in_common] where notes != 'ignore' and common_name != '') AS p;

/* concatenate the sql statement */
SET @sql =
N'create view clinical_trait_differences_male as
SELECT hf.strain, ' + @col + 
' 
from dbo.avg_clinical_traits_by_strain_chow_male as c
inner join 
dbo.avg_clinical_traits_by_strain_highfat_male as hf
on c.strain = hf.strain' 
;

/* execute sql statement */
EXEC sp_executesql @sql;

