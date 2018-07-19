# Getting started with Tachyon

## Table of contents
- [Workflow example: using the CLI](#workflow-example-using-the-cli)
    - [`import`: Importing `VCF`/`BCF`](#import-importing-vcfbcf)
    - [`view`: Viewing, converting, and slicing `YON` files](#view-viewing-converting-and-slicing-yon-files)
    - [Field-slicing](#field-slicing)
    - [Searching for genomic regions](#searching-for-genomic-regions)
    - [Annotating meta-data](#annotating-meta-data)
- [C++ API Examples](#c-api-examples)

---

## Workflow example: using the CLI
### `import`: Importing `VCF`/`BCF`
Import a `bcf` file to `yon` with a block-size of `-c` number of variants and/or `-C` number of base-pairs. If both `-c` and `-C` are set then the block breaks whenever either condition is satisfied. **Please note that importing VCF files are currently disabled**
```bash
tachyon import -i examples/example_dataset.bcf -o example_dataset.yon -c 2000
```

Tachyon can protect your sensitive identifying information with high-grade encryption. By default, each data field is encrypted separately in each block with different keys using AES-256. Simply pass the `-e` flag and the best practices will be used.
```bash
tachyon import -i examples/example_dataset.bcf -o example_dataset.yon -c 2000 -e
```
This will produce two output files:
* example_dataset.yon
* example_dataset.kyon

### `view`: Viewing, converting, and slicing `YON` files
Printing a `yon` file as a bit-exact copy of the input `VCF`
```bash
tachyon view -i example_dataset.yon -H
```
Output
```
Contig110_arrow	672	.	A	T	525.07	basic_filtering	AC=10;AF=0.217;AN=46;BaseQRankSum=0.967;DP=72;ExcessHet=0.8113;FS=54.73;InbreedingCoeff=-0.0525;MLEAC=11;MLEAF=0.239;MQ=31.05;MQRankSum=1.38;QD=18.11;ReadPosRankSum=-0.431;SOR=5.889	GT:AD:DP:GQ:PL	0/0:1,0:1:3:0,3,38	1/1:0,2:.:6:76,6,0	./.:0,0:0:.:0,0,0	0/1:2,2:.:43:58,0,43	0/0:1,0:1:3:0,3,24	0/0:4,0:4:12:0,12,141	./.:0,0:0:.:0,0,0	0/0:1,0:1:3:0,3,29	0/1:1,4:.:19:147,0,19	0/1:3,2:.:49:71,0,49	0/0:5,0:5:0:0,0,81	./.:1,0:1:.:0,0,0	./.:0,0:0:.:0,0,0	0/0:6,0:6:18:0,18,192	1/1:0,2:.:6:58,6,0	0/0:1,0:1:3:0,3,11	./.:1,0:1:.:0,0,0	./.:0,0:0:.:0,0,0	0/1:3,2:.:58:58,0,63	./.:0,0:0:.:0,0,0	0/0:4,0:4:12:0,12,134	0/0:3,0:3:0:0,0,44	0/0:3,0:3:9:0,9,90	./.:0,0:0:.:0,0,0	./.:0,0:0:.:0,0,0	./.:0,0:0:.:0,0,0	0/0:2,0:2:6:0,6,53	0/1:1,2:.:19:62,0,19	0/0:3,0:3:9:0,9,84	0/0:2,0:2:6:0,6,49	0/1:1,2:.:19:74,0,19	0/0:1,0:1:3:0,3,38	./.:2,0:2:.:0,0,0	./.:0,0:0:.:0,0,0	0/0:2,0:2:6:0,6,65	./.:0,0:0:.:0,0,0
```
We can check for bit-exact output from `tachyon` by comparing the output of the cryptographic hash function `SHA512` for `bcftools`. We drop the header `-H` as these two are different: both tools inject a timestamp and library versions each time a command is executed.
```bash
tachyon view -i example_dataset.yon -H | openssl dgst -sha512
```
```
4c94ee35fa3509935e5ea63f6da9b39dc94b1073b551c7d4d56bca7666a6872ad629b6f91f43a8dc45b306c0b0bbb2f414fb811ed45c7e6434c3570b2e448c68
```
```bash
bcftools view example_dataset.bcf -H | openssl dgst -sha512
```
```
4c94ee35fa3509935e5ea63f6da9b39dc94b1073b551c7d4d56bca7666a6872ad629b6f91f43a8dc45b306c0b0bbb2f414fb811ed45c7e6434c3570b2e448c68
```

Listing only site-specific information and `INFO` fields:
```bash
tachyon view -i example_dataset.yon -GH
```
Output
```
Contig110_arrow	672	.	A	T	525.07	basic_filtering	AC=10;AF=0.217;AN=46;BaseQRankSum=0.967;DP=72;ExcessHet=0.8113;FS=54.73;InbreedingCoeff=-0.0525;MLEAC=11;MLEAF=0.239;MQ=31.05;MQRankSum=1.38;QD=18.11;ReadPosRankSum=-0.431;SOR=5.889
```

### Field-slicing
Listing a specific `INFO` field with output data still adhering to the `VCF` specification:
```bash
tachyon view -i example_dataset.yon -GH -f "INFO=AC"
```
Output
```
Contig110_arrow	672	.	.	.	.	basic_filtering	AC=10
```

Add `REF` and `ALT` to the output
```bash
tachyon view -i example_dataset.yon -GH -f "INFO=AC;REF;ALT"
```
Output
```
Contig110_arrow	19575	.	A	C	.	basic_filtering	AC=2
```

Listing this output in a custom tab-delimited format
```bash
tachyon view -i example_dataset.yon -GH -f "INFO=AC;REF;ALT" -y -d'\t'
```
Output (first five lines)
```
A	T	10
A	T	36
TG	T	34
G	A	2
G	C	2
```

Listing `CHROM`, `POS`, and all `INFO` fields
```bash
tachyon view -i example_dataset.yon -GH -f "CHROM;POS;INFO" -y -d'|'
```
Output
```
Contig110_arrow|672|10|0.217|46|0.967|72|0.8113|54.73|-0.0525|11|0.239|31.05|1.38|18.11|-0.431|5.889
```

Listing all available `INFO` fields and the `FORMAT` fields `DP` and `PL` in VCF 
```bash
tachyon view -i example_dataset.yon -f "chrom;pos;ref;alt;info;format=dp,pl" -H -O vcf
```
Output
```
Contig110_arrow	672	.	A	T	.	basic_filtering	AC=10;AF=0.217;AN=46;BaseQRankSum=0.967;DP=72;ExcessHet=0.8113;FS=54.73;InbreedingCoeff=-0.0525;MLEAC=11;MLEAF=0.239;MQ=31.05;MQRankSum=1.38;QD=18.11;ReadPosRankSum=-0.431;SOR=5.889	DP:PL	1:0,3,38	.:76,6,0	0:0,0,0	.:58,0,43	1:0,3,24	4:0,12,141	0:0,0,0	1:0,3,29	.:147,0,19	.:71,0,49	5:0,0,81	1:0,0,0	0:0,0,0	6:0,18,192	.:58,6,0	1:0,3,11	1:0,0,0	0:0,0,0	.:58,0,63	0:0,0,0	4:0,12,134	3:0,0,44	3:0,9,90	0:0,0,0	0:0,0,0	0:0,0,0	2:0,6,53	.:62,0,19	3:0,9,84	2:0,6,49	.:74,0,19	1:0,3,38	2:0,0,0	0:0,0,0	2:0,6,65	0:0,0,0
```

Listing all available `INFO` fields and the `FORMAT` fields `DP` and `PL` in JSON
```bash
tachyon view -i example_dataset.yon -f "chrom;pos;ref;alt;info;format=dp,pl" -O json
```
Output
```json
"block": {
	"obj-0": {
		"contig": "Contig110_arrow",
		"position": 672,
		"ref": "A",
		"alt": "T",
		"INFO-AC": 10,
		"INFO-AF": 0.217,
		"INFO-AN": 46,
		"INFO-BaseQRankSum": 0.967,
		"INFO-DP": 72,
		"INFO-ExcessHet": 0.8113,
		"INFO-FS": 54.73,
		"INFO-InbreedingCoeff": -0.0525,
		"INFO-MLEAC": 11,
		"INFO-MLEAF": 0.239,
		"INFO-MQ": 31.05,
		"INFO-MQRankSum": 1.38,
		"INFO-QD": 18.11,
		"INFO-ReadPosRankSum": -0.431,
		"INFO-SOR": 5.889,
		"FORMAT-DP": [1, null, 0, null, 1, 4, 0, 1, null, null, 5, 1, 0, 6, null, 1, 1, 0, null, 0, 4, 3, 3, 0, 0, 0, 2, null, 3, 2, null, 1, 2, 0, 2, 0],
		"FORMAT-PL": [
			[0, 3, 38],
			[76, 6, 0],
			[0, 0, 0],
			[58, 0, 43],
			[0, 3, 24],
			[0, 12, 141],
			[0, 0, 0],
			[0, 3, 29],
			[147, 0, 19],
			[71, 0, 49],
			[0, 0, 81],
			[0, 0, 0],
			[0, 0, 0],
			[0, 18, 192],
			[58, 6, 0],
			[0, 3, 11],
			[0, 0, 0],
			[0, 0, 0],
			[58, 0, 63],
			[0, 0, 0],
			[0, 12, 134],
			[0, 0, 44],
			[0, 9, 90],
			[0, 0, 0],
			[0, 0, 0],
			[0, 0, 0],
			[0, 6, 53],
			[62, 0, 19],
			[0, 9, 84],
			[0, 6, 49],
			[74, 0, 19],
			[0, 3, 38],
			[0, 0, 0],
			[0, 0, 0],
			[0, 6, 65],
			[0, 0, 0]
		]
	},
```

Listing all available `INFO` fields and the `FORMAT` fields `DP` and `PL` in a custom output format with a tab-delimiter. No restrictions are placed on the output format
```bash
tachyon view -i example_dataset.yon -f "pos;info;format=dp,pl" -O custom -yd'\t'
```
Output
```
672	10	0.217	46	0.967	72	0.8113	54.73	-0.0525	11	0.239	31.05	1.38	18.11	-0.431	5.889	DP:PL	1:0,3,38	.:76,6,0	0:0,0,0	.:58,0,43	1:0,3,24	4:0,12,141	0:0,0,01:0,3,29	.:147,0,19	.:71,0,49	5:0,0,81	1:0,0,0	0:0,0,0	6:0,18,192	.:58,6,0	1:0,3,11	1:0,0,0	0:0,0,0	.:58,0,63	0:0,0,0	4:0,12,134	3:0,0,44	3:0,9,90	0:0,0,0	0:0,0,00:0,0,0	2:0,6,53	.:62,0,19	3:0,9,84	2:0,6,49	.:74,0,19	1:0,3,38	2:0,0,0	0:0,0,0	2:0,6,65	0:0,0,0
```

Listing all available `INFO` fields and the `FORMAT` fields `DP` and `PL` in a custom output format with a tab-delimiter and
all `FORMAT` data are printed as vectors of samples instead of per-sample
```bash
tachyon view -i example_dataset.yon -f "chrom;pos;ref;alt;info;format=dp,pl" -O custom -yd'\t' -V
```
Output
```
Contig110_arrow	672	A	T	10	0.217	46	0.967	72	0.8113	54.73	-0.0525	11	0.239	31.05	1.38	18.11	-0.431	5.889	DP:PL	1,.,0,.,1,4,0,1,.,.,5,1,0,6,.,1,1,0,.,0,4,3,3,0,0,0,2,.,3,2,.,1,2,0,2,0	0,3,38,76,6,0,0,0,0,58,0,43,0,3,24,0,12,141,0,0,0,0,3,29,147,0,19,71,0,49,0,0,81,0,0,0,0,0,0,0,18,192,58,6,0,0,3,11,0,0,0,0,0,0,58,0,63,0,0,0,0,12,134,0,0,44,0,9,90,0,0,0,0,0,0,0,0,0,0,6,53,62,0,19,0,9,84,0,6,49,74,0,19,0,3,38,0,0,0,0,0,0,0,6,65,0,0,0
```

Output all available `INFO` fields and the `FORMAT` field `DP` and `FILTERS`
```bash
tachyon view -i example_dataset.yon -f "chrom;pos;ref;alt;info;format=dp;filter" -O custom -yd';' -V
```

### Searching for genomic regions
Slicing intervals either as a contig, contig with a single position, or interval with a contig:
```bash
tachyon view -i example_dataset.yon -r "Contig110_arrow"
```
```bash
tachyon view -i example_dataset.yon -r "Contig110_arrow:672"
```
```bash
tachyon view -i example_dataset.yon -r "Contig110_arrow:672-1500"
```

### Annotating meta-data
It is possible to annotate data with a series of `INFO` fields computed directly from the genotypic vectors or from the reference/alternative allele data:

| Field           | Length | Type    | Description                                 |
|-----------------|--------|---------|---------------------------------------------|
| `FS_A`          | `A`      | `Float`   | PHRED-scaled Fisher's exact test P-value for allelic strand bias |
| `AN`            | `A`      | `Integer` | Total number of alleles in called genotypes |
| `NM`            | `A`      | `Integer` | Total number of missing alleles in called genotypes |
| `NPM`           | `A`      | `Integer` | Total number of samples with non-reference (compared to largest) ploidy |
| `AC`            | `A`      | `Integer` | Total number of alleles |
| `AC_FWD`        | `A`      | `Integer` | Total number of alleles on the *forward* strand |
| `AC_REV`        | `A`      | `Integer` | Total number of alleles on the *reverse* strand |
| `AF`            | `A`      | `Float`   | Allele frequency of allele |
| `HWE_P`         | `A`      | `Float`   | Hardy-Weinberg equilibrium P-value |
| `VT`            | `A`      | `String`  | Variant classification (SNP, MNP, INDEL, CLUMPED, SV, UNKNOWN) |
| `MULTI_ALLELIC` | 0        | `Flag`    | Indicates if a site is multi-allelic (number of alternative alleles > 1) |
| `F_PIC`         | `A`      | `Float`   | Population inbreeding coefficient (F-statistic) |

The contingency table, or matrix, for the Fisher's exact test (`FS_A` ) for strand bias looks like this:

|                | Target allele | *Not* target allele |
|----------------|---------------|-------------------|
| Forward strand | A             | B                 |
| Reverse strand | C             | D                 |

In the biallelic case, only one P-value is reported becuase of symmetry. If the site is not biallelic then each individual allele is computed separately. 

Using the example dataset, we can compute those fields that are not already available by passing the flag `-X`

```bash
tachyon view -i example_dataset.yon -GHX
```
Output
```
Contig110_arrow	672	.	A	T	525.07	basic_filtering	AC=10;AF=0.217;AN=46;BaseQRankSum=0.967;DP=72;ExcessHet=0.8113;FS=54.73;InbreedingCoeff=-0.0525;MLEAC=11;MLEAF=0.239;MQ=31.05;MQRankSum=1.38;QD=18.11;ReadPosRankSum=-0.431;SOR=5.889;FS_A=11.5091,0;NM=26;AC_FWD=21,2;AC_REV=15,8;HWE_P=0.25072;VT=SNP
```

## C++ API Examples
We provide several API examples in the `lib_example` directory. Get started by `make examples`: this requires you to have compiled the shared library first.
