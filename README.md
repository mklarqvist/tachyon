[![Release](https://img.shields.io/badge/Release-beta_0.1-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

# Tachyon
Tachyon is an efficient API to store, query, and handle big data resequencing data. Tachyon efficiently stores data tables by column and compress genotypic data by exploiting intrinsic genetic properties. We describe algorithms to directly query, manipulate, and explore this jointly compressed genotypic representation in-place.

Marcus D. R. Klarqvist (<mk819@cam.ac.uk>)

### Installation instructions
Building requires [zstd][zstd] and [openssl][openssl]
```bash
git clone --recursive https://github.com/mklarqvist/Tachyon
cd Tachyon/build
make
```
## <a name="notes"></a>Further Notes
### <a name="others"></a>Other formats

* Tachyon vs [PBWT][pbwt]. Compression of genotype data in Tachyon is inspired by the
  2-bin radix sort in PBWT. Tachyon, unlike PBWT, keeps all multi-allelic variants
  and all INFO/FORMAT/FILTER fields. PBWT supports advanced queries for 
  haplotype matching, phasing and imputation. Tachyon has no such functionality.

* Tachyon vs [BCF2][vcf]. Tachyon supports all functionality in BCF. BCF keeps genotype data internally in roughly
  the same space as raw VCF. This makes querying BCF extremely slow. Tachyon stores
  genotypic data as (dynamic) fixed-width run-length encoded objects that we can
  query directly. 

* Tachyon vs [GQT][gqt]. GQT is designed with speed in mind and completely sacrifice
  compressibility. Because of this, GQT archives are >20-fold larger than Tachyon and
  many hundred-fold larger than PBWT, BGF, and GTC. Tachyon is just as fast at processing
  genotype data while compressing much better.

* Tachyon vs [GTC][gtc]. GTC provides superior compressibility of genotypes but is considerably
  less expressible as it throws away all other fields. GTC will struggle to rapidly 
  query data because of the bit-vector representation it employs.
  
### <a name="comp"></a>Compression evaluation
The test is run on the first release of [Haplotype Reference Consortium][hrc]
(HRC) data. There are ~39 million phased SNPs in 32,488 samples. We have
generated the `YON` for chromosome 11. The following table shows the file sizes
across a number of tools.

| File-size         | Command line                                                          |
| ----------------: | :-------------------------------------------------------------------- |
| 3322M             | `bcftools view HRC_chr11.vcf -o hrc_chr11.bcf`                        |
| 323M + 965K       | `tachyon import -i HRC_chr11.bcf -o hrc_chr11`                       |
| 185M + 30M + 679K | `pbwt -readVcfGT HRC_chr11.bcf -checkpoint 10000 -writeAll hrc_chr11` |
| 188M + 12M + 62K  | `gtc compress -b HRC_chr11.bcf -o hrc_chr11`                          |
| 349M + 12M + 68K  | `bgt import hrc_chr11 HRC_chr11.bcf`                                  |
| 4110M + 58M + 7M  | `gqt convert bcf -i HRC_chr11.bcf`                                    |
| 4110M + 58M + 7M  | `seqarray convert bcf -i HRC_chr11.bcf`                                    |

Listed files represent the primary archive plus any auxiliary files (1M=1024\*1024 bytes).

[hrc]: http://www.haplotype-reference-consortium.org
[gqt]: https://github.com/ryanlayer/gqt
[pbwt]: https://github.com/richarddurbin/pbwt
[gtc]: https://github.com/refresh-bio/GTC
[vcf]: https://samtools.github.io/hts-specs/
[sa]: https://github.com/zhengxwen/SeqArray
[openssl]: https://www.openssl.org/
[zstd]: https://github.com/facebook/zstd

### License
[MIT](LICENSE)
