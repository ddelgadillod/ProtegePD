# PROTÉGÉ PD

**PROTE**in coding **GE**ne for phylogenetic tag and identification - **P**rimer **D**esign tool

Based on [Phylotag approach (Caro-Quintero, 2015)](https://academic.oup.com/gbe/article/7/12/3416/2467318).

This Protégé version runs in a Docker container, including Python 3.10, with a minimal interactive graphic user interface deployed with Dash.

## Requisites
### Docker

Please follow Docker's official installation instructions according to your Operating System:

- [Linux](https://docs.docker.com/desktop/install/linux-install/)
- [Windows](https://docs.docker.com/desktop/install/windows-install/)
- [Mac](https://docs.docker.com/desktop/install/mac-install/)
 
Check the installation process by following the official [Get started tutorial](https://docs.docker.com/desktop/get-started/).

## Downloading Image

For this Protégé version, access to your Operating System Command Line Interface (CLI) is mandatory.

Push the image from Docker Hub cloud by running:

```
docker pull ddelgadillo/protege_base:v1.0.1
```

## Run Container

In the CLI, run the following command:

```
docker run --rm --mount type=bind,source=/your/files/path/,target=/root/. --name protege -p 127.0.0.1:8050:8050 --cpus 4 protege_base:latest protege-pd -s myseqs.fna
```

Edit the previous command according to your files:

1. Replace the `/your/files/path/` field with the full path where your fasta files are located.
2. Replace the `myseqs.fna` field according to the name of your fasta file.

### Running parameters
It is possible to check additional running parameters by displaying the help:

```
docker run --rm --mount type=bind,source=.,target=/root/. --name protege -p 127.0.0.1:8050:8050 --cpus 4 protege_base:v1.0.1 protege-pd --help

usage: protege-pd [-h] -s SEQ [-c CONS] [-g] [-d COD] [-v]

'PROTÉGÉ' PROTEin coding GEne for phylogenetic tag and identification. Useful design and visualization of primers.

options:
  -h, --help            show this help message and exit
  -s SEQ, --seq SEQ     Fasta file with genes sequences (default: None)
  -c CONS, --consensus CONS
                        Consensus percentage (default: 90)
  -g, --nogapconsensus  Do not consider consensus with gaps (default: True)
  -d COD, --codon COD   Codon primer length (default: 7)
  -v, --verbose         increase output verbosity (default: False)
```

## Contact
Suggestions, questions, comments, contributions, and issue reports are welcome. Please contact me via [email](mailto:diegodelgadilloduran@gmail.com) or feel free to open an issue.

