# Python code related to networks from Daniel Himmelstein's PhD in the Sergio Baranzini Lab

This repository contains Python code, written by Daniel Himmelstein from 2012-2015 during his PhD in the [Baranzini Lab](https://baranzinilab.ucsf.edu/) at UCSF.
This includes code that contributed to the following studies:

<!-- 
manubot cite --format=markdown \
  doi:10.1371/journal.pcbi.1004259 \
  doi:10.12688/f1000research.6836.2 \
  doi:10.6084/m9.figshare.4724797
-->

1. **Heterogeneous Network Edge Prediction: A Data Integration Approach to Prioritize Disease-Associated Genes**   
Daniel S. Himmelstein, Sergio E. Baranzini  
*PLOS Computational Biology* (2015-07-09) <https://doi.org/98q>   
DOI: [10.1371/journal.pcbi.1004259](https://doi.org/10.1371/journal.pcbi.1004259) · PMID: [26158728](https://www.ncbi.nlm.nih.gov/pubmed/26158728) · PMCID: [PMC4497619](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4497619)

2. **iCTNet2: integrating heterogeneous biological interactions to understand complex traits**   
Lili Wang, Daniel S. Himmelstein, Adam Santaniello, Mousavi Parvin, Sergio E. Baranzini  
*F1000Research* (2015-09-28) <https://doi.org/ghfv3v>   
DOI: [10.12688/f1000research.6836.2](https://doi.org/10.12688/f1000research.6836.2) · PMID: [26834985](https://www.ncbi.nlm.nih.gov/pubmed/26834985) · PMCID: [PMC4706053](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4706053)

3. **The hetnet awakens: understanding complex diseases through data integration and open science** [Thesis]   
Daniel Himmelstein  
*figshare* (2017) <https://doi.org/b2nz>   
DOI: [10.6084/m9.figshare.4724797](https://doi.org/10.6084/m9.figshare.4724797) · ISBN: [9781339919881](https://worldcat.org/isbn/9781339919881)

Code related to [Project Rephetio](https://doi.org/10.7554/eLife.26726 "Systematic integration of biomedical knowledge prioritizes drugs for repurposing. eLife. 2017") is released separately, and largely post-dates this repository.

The code to create the PLOS Comp Bio (2015) network, i.e. `het.io-dag`, is in [`projects/gene_disease_hetnet/createnet.py`](projects/gene_disease_hetnet/createnet.py) ([access data here](https://github.com/dhimmel/het.io-dag-data/tree/master/downloads#network)).

## History

This repository originally used mercurial and was hosted on BitBucket at `https://bitbucket.org/dhimmel/serg`.
Unfortunately, BitBucket [deleted](https://community.atlassian.com/t5/Bitbucket-articles/What-to-do-with-your-Mercurial-repos-when-Bitbucket-sunsets/ba-p/1155380) mercurial repos in 2020 without an automated migration path.
On 2020-10-15, dhimmel realized this source code was no longer online.

The repo could not be retrieved from BitBucket, but dhimmel had a local copy in `~/Documents/serg/pycode`.
The local copy matched the version archived by [Software Heritage](https://archive.softwareheritage.org/browse/directory/2361a3289af5749822df1c2b5addb5c73f36cab1/?branch=default&origin_url=https://bitbucket.org/dhimmel/serg&snapshot=179679a99b43c5549e9f89bafc6551200d3bb712).
Using the [hg-fast-export](https://github.com/frej/fast-export/tree/ead75895b058d16ffc7330dab78054c94a189377) utility, dhimmel exported the repository to git and uploaded it to <https://github.com/dhimmel/serg-pycode>.

This code is uploaded to GitHub for archiving and historical purposes.
The extent this python code was used to perform certain studies and analyses requires a bit of investigation.
The repository is not actively maintained nor accepting feature requests.
However, questions about existing code or analyses are welcome via GitHub Issues.
Furthermore, any documentation contributions are appreciated.
