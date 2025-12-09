# How to run the NanoCliD pipeline

1. suffer
2. start interactive job
3. load miniconda 24.9.2?
4. activate nanoclid environment
5. activate venv 
    1. cd to `/hpc/local/Rocky8/pmc_kuiper/software/NanoCliD/NanoCliD`
    2. `source venv/bin/activate`
6. set the correct tmp dir
    1. `export TMPDIR=/tmp`



    you need a hg19.fa, the index, a dictionary (made with picardtools), and a genome file. maybe even more