# genomicduplicationr

### Short description
Cluster single gene duplications into genomic duplications for input rooted gene trees and a species tree.

### Quick install
Create virtual environment `python -m venv gd`

Activate the environment `source gd/bin/activate`

Install package `pip install git+https://github.com/j-paszek/genomicduplicationr.git`

Use program (-h for help) `gdscore -h`

### Sample usage

Download sample input file `https://raw.githubusercontent.com/j-paszek/genomicduplicationr/main/tests/data/RME/in.txt`

Above file consists of a set of rooted binary gene trees and a rooted binary species tree.

To compute genomic duplication score simply run `gdscore in.txt`

Default genomic duplications are clustered by minimum episodes (ME) clustering. 
To use episode clustering (EC) method run `gdscore -e in.txt` 

Default model of duplication mappings (intervals) is PG model from: \
_Paszek, J., Gorecki, P., "Efficient algorithms for genomic duplication models"
IEEE/ACM Trans. Comput. Biol. Bioinform. 15 (5), 1515-1524, 2018._\
To specify model choose (-L, -G, -P, -F) for LCA, GMS, PG, and FHS (unrestricted) model.\
In example: To compute EC score for FHS model use `gdscore -e -F in.txt`

### license

[![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].



[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

<!---
Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
-->
