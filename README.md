# Assignment 2

## Requirements

Python
```
pip install httpx tqdm
```

System (linux)
```
 $ apt install ncbi-blast+

```

## How to use

Create a folder to put all the files

```
mkdir assignment2
```

Go to the created folder

```
cd assignment2
```

Clone the repo

```
git clone git@github.com:atiliopereira/psi-blast-2024.git
```


## Results

### Task 1.2

[`2229317dabd1c78fc85da327dc3c46255402bf0e`](https://github.com/atiliopereira/psi-blast-2024/commit/2229317dabd1c78fc85da327dc3c46255402bf0e) contains the result of running:

```
python3 fetch_fasta_skeleton.py -i data/SCOP_selections.txt -db data/myDatabase.fasta
```

[`27e2625a0d18a58b7ec3131e5c9ada26c17aa6f4`](https://github.com/atiliopereira/psi-blast-2024/commit/27e2625a0d18a58b7ec3131e5c9ada26c17aa6f4) contains the result of running:

```
makeblastdb -in ./data/myDatabase.fasta -parse_seqids -dbtype prot
```

### Task 2.2

[`553c197202901719e074fcf4071bcc1b42bfb7be`](https://github.com/atiliopereira/psi-blast-2024/commit/553c197202901719e074fcf4071bcc1b42bfb7be) contains the result of running:

```
python3 run_local_skeleton.py -i data/SCOP_selections.txt -db data/myDatabase.fasta -o results/blastEvalues.tsv
```

[`b3c62068745fb38d8e347383e9cc099a6e7ec26b`](https://github.com/atiliopereira/psi-blast-2024/commit/b3c62068745fb38d8e347383e9cc099a6e7ec26b) contains the result of running:

```
python3 run_local_skeleton.py -i data/SCOP_selections.txt -db data/myDatabase.fasta -o results/psi-blastEvalues.tsv -psi 
```

### Question 3.1
[`b0c69ca3aeea8b8424d100a8593e4ed70c581ca2`](https://github.com/atiliopereira/psi-blast-2024/commit/b0c69ca3aeea8b8424d100a8593e4ed70c581ca2) contains the result of running:

```
python3 run_local_skeleton.py -r ./results/blastEvalues.tsv -o ./results/BLAST_histogram.png
```


[`dc6c702dbcbbcc558b06d539197f5ca52371ef67`](https://github.com/atiliopereira/psi-blast-2024/commit/dc6c702dbcbbcc558b06d539197f5ca52371ef67) contains the result of running:

```
python3 run_local_skeleton.py -r ./results/blastEvalues.tsv -o ./results/PSIBLAST_histogram.png -psi
```