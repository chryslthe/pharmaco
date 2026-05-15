# HG00367: comparison with PharmCAT

Сравнение выполнено на одном и том же VCF sample `HG00367`.

Главный вывод: для простых SNP/star-allele случаев результаты хорошо совпадают. Полное совпадение наблюдается для `CYP2C19`, `CYP3A5`, `CYP2B6`, `CYP2C9`, `CYP3A4`, `ABCG2`, `TPMT`. Направление совпадает для `VKORC1` и `NAT2`. Расхождения ожидаемы для сложных случаев: `CFTR`, `CYP2D6`, `HLA-B`, `G6PD`, потому что они требуют более полного panel/caller/CNV/HLA typing.

Файлы:

- `outputs/pharmcat_comparison_HG00367.csv`
- `outputs/pharmcat_comparison_HG00367.html`
- `figures/pharmcat_comparison_HG00367.svg`
