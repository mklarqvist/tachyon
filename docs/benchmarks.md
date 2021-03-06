## Benchmarks
### Simulations
We simulated haplotypes using [msprime][msprime] for a 50 megabase region for varying number of individuals (fixed parameters: recombination rate 2e-8, mutation rate 2e-8, effective population size 1e4) and concatened pairs of haplotypes together to form diploid genotypes. Decompression times were measured by `time zcat <file.bcf> > /dev/null` for `bcf` and using the API for Tachyon and timings for printing fixed-fields were measured as `time bcftools view <file.bcf> -GH > /dev/null` for `bcf` and `time tachyon view -GH > /dev/null` for Tachyon. File sizes are listed in gigabytes (1 GB = 10e9 b) and timings in seconds. All experiments were run on a single CPU.  

| Variants  | Samples  | Filesize (BCF) | Filesize (YON) | Decomp. (BCF) | Decomp. (YON) | Print sites (BCF) | Print sites (YON) |
|-----------|----------|----------------|----------------|---------------|---------------|-------------------|-------------------|
| 523,842   | 10,000   | 0.53655        | 0.16043        | 44.295        | 1.859         | 28.098            | 0.680             |
| 604,487   | 50,000   | 2.57511        | 0.59032        | 251.682       | 7.077         | 137.200           | 1.609             |
| 639,666   | 100,000  | 5.12134        | 0.99690        | 526.733       | 12.473        | 286.413           | 1.936             |
| 685,363   | 250,000  | 12.64719       | 1.90681        | 1927.011      | 24.671        | 1028.046          | 4.745             |
| 719,754   | 500,000  | 25.08209       | 3.04404        | 4139.424      | 45.247        | 2241.756          | 11.706            |  

![screenshot](sim_50gbp.jpeg)


### Real datasets
The following table shows data for the 1000 Genomes Project Phase 3 release (2,504 samples, in megabytes; 1 MB = 1E6 bytes). The uncompressed file size represents the amount of bytes needed to be parsed internally   

| Contig | BCF-compressed | BCF-uncompressed | YON-compressed | YON-uncompressed | YON-fold | BCF-fold | Uncompressed-fold |
|--------|----------------|------------------|----------------|------------------|----------|----------|-------------------|
| 1      | 1004.15        | 33140.9          | 258.842        | 1263.27          | 128.04   | 33.004   | 26.234            |
| 2      | 1084.63        | 36284.2          | 275.986        | 1347.19          | 131.47   | 33.453   | 26.933            |
| 3      | 914.54         | 29883.5          | 229.379        | 1132.73          | 130.28   | 32.676   | 26.382            |
| 4      | 922.15         | 29373.0          | 225.035        | 1146.94          | 130.53   | 31.853   | 25.610            |
| 5      | 817.85         | 26980.8          | 204.648        | 1014.88          | 131.84   | 32.990   | 26.585            |
| 6      | 826.64         | 25743.2          | 201.314        | 1016.07          | 127.88   | 31.142   | 25.336            |
| 7      | 750.82         | 24167.4          | 192.220        | 937.06           | 125.73   | 32.188   | 25.791            |
| 8      | 711.66         | 23554.5          | 180.243        | 884.54           | 130.68   | 33.098   | 26.629            |
| 9      | 556.59         | 18244.0          | 149.945        | 704.78           | 121.67   | 32.778   | 25.886            |
| 10     | 639.80         | 20455.4          | 164.439        | 799.29           | 124.40   | 31.971   | 25.592            |
| 11     | 633.70         | 20728.9          | 159.191        | 783.82           | 130.21   | 32.711   | 26.446            |
| 12     | 613.10         | 19821.2          | 157.542        | 766.12           | 125.82   | 32.330   | 25.872            |
| 13     | 460.82         | 14643.6          | 116.473        | 575.91           | 125.72   | 31.777   | 25.427            |
| 14     | 418.90         | 13604.1          | 108.813        | 525.03           | 125.02   | 32.476   | 25.911            |
| 15     | 378.46         | 12423.5          | 104.407        | 479.29           | 118.99   | 32.826   | 25.921            |
| 16     | 408.44         | 13823.2          | 116.201        | 525.16           | 118.96   | 33.843   | 26.322            |
| 17     | 358.58         | 11934.5          | 100.589        | 459.74           | 118.65   | 33.283   | 25.959            |
| 18     | 360.95         | 11616.7          | 96.252         | 462.71           | 120.69   | 32.184   | 25.106            |
| 19     | 295.92         | 9389.2           | 83.711         | 387.25           | 112.16   | 31.728   | 24.246            |
| 20     | 282.46         | 9288.6           | 77.330         | 364.36           | 120.12   | 32.885   | 25.493            |
| 21     | 180.55         | 5664.7           | 50.263         | 236.30           | 112.70   | 31.374   | 23.972            |
| 22     | 177.18         | 5654.2           | 51.468         | 229.95           | 109.86   | 31.912   | 24.589            |

![screenshot](1kgp3_yon_bcf.jpeg)  
The following table shows data for the Haplotype Reference Consortium (32,488 whole-genome sequenced samples)  

| Contig | BCF-compressed | BCF-uncompressed | YON-compressed | YON-uncompressed | YON_fold | BCF_fold | Uncompressed_fold |
|--------|----------------|------------------|----------------|------------------|----------|----------|-------------------|
| 1      | 5359.20        | 199628           | 895.46         | 5865.5           | 222.93   | 37.250   | 34.034            |
| 2      | 5855.74        | 220586           | 942.84         | 6242.9           | 233.96   | 37.670   | 35.334            |
| 3      | 4991.04        | 183499           | 796.51         | 5355.0           | 230.38   | 36.766   | 34.266            |
| 4      | 5078.87        | 181268           | 763.12         | 5488.2           | 237.53   | 35.691   | 33.029            |
| 5      | 4533.42        | 168300           | 714.14         | 4861.6           | 235.67   | 37.124   | 34.619            |
| 6      | 4588.08        | 159973           | 692.21         | 4807.6           | 231.10   | 34.867   | 33.275            |
| 7      | 4101.22        | 148866           | 676.84         | 4536.9           | 219.94   | 36.298   | 32.812            |
| 8      | 3902.68        | 145836           | 635.00         | 4224.8           | 229.66   | 37.368   | 34.519            |
| 9      | 3000.07        | 109666           | 555.14         | 3456.6           | 197.55   | 36.554   | 31.727            |
| 10     | 3536.07        | 125972           | 593.09         | 3935.3           | 212.40   | 35.625   | 32.010            |
| 11     | 3476.07        | 125957           | 561.41         | 3689.2           | 224.36   | 36.235   | 34.142            |
| 12     | 3330.20        | 120177           | 556.28         | 3661.9           | 216.04   | 36.087   | 32.819            |
| 13     | 2571.01        | 90091            | 421.10         | 2859.1           | 213.94   | 35.041   | 31.510            |
| 14     | 2270.78        | 82613            | 386.20         | 2496.6           | 213.91   | 36.381   | 33.090            |
| 15     | 1996.39        | 74080            | 379.79         | 2278.8           | 195.05   | 37.107   | 32.508            |
| 16     | 2181.47        | 83319            | 451.09         | 2612.2           | 184.70   | 38.194   | 31.896            |
| 17     | 1874.08        | 70884            | 365.79         | 2153.1           | 193.78   | 37.823   | 32.921            |
| 18     | 2001.76        | 71839            | 364.10         | 2321.1           | 197.31   | 35.888   | 30.950            |
| 19     | 1559.36        | 56480            | 312.13         | 1887.2           | 180.95   | 36.220   | 29.928            |
| 20     | 1573.38        | 57548            | 303.62         | 1878.9           | 189.54   | 36.576   | 30.628            |
| 21     | 986.49         | 34548            | 186.18         | 1221.8           | 185.56   | 35.021   | 28.275            |
| 22     | 953.16         | 34110            | 196.84         | 1156.9           | 173.29   | 35.786   | 29.484            |

![screenshot](examples/hrc_bcf_yon.jpeg)  


[msprime]:  https://github.com/jeromekelleher/msprime