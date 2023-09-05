# MutAlign
Creating HTML visualizations of [VESPA](https://github.com/Rostlab/VESPA) mutation effect predictions based on sequence alignments

![Visualization of aligned mutation effect predictions](https://rostlab.org/~conpred/MutAlign.png)

This tool is useful, if you want to visually compare [VESPA](https://github.com/Rostlab/VESPA) mutation effect predictions for highly similar sequences (e.g. the same protein in different organisms).

# Requirements

MutAlign itself requires BioPython:
```
pip install biopython
```
Refer to the [VESPA](https://github.com/Rostlab/VESPA) documentation for installation instructions.
You will also need [MMSeqs2](https://github.com/soedinglab/MMseqs2) for sequence alignments.

# Usage example

Let's say you are interested in a comparison of mutation effects in the POLK gene (DNA polymerase) in human and mice. Obtain the two amino-acid sequences and store them in `POLK_HUMAN.fasta` and `POLK_MOUSE.fasta`:

```
>POLK_HUMAN
MDSTKEKCDSYKDDLLLRMGLNDNKAGMEGLDKEKINKIIMEATKGSRFYGNELKKEKQV
NQRIENMMQQKAQITSQQLRKAQLQVDRFAMELEQSRNLSNTIVHIDMDAFYAAVEMRDN
PELKDKPIAVGSMSMLSTSNYHARRFGVRAAMPGFIAKRLCPQLIIVPPNFDKYRAVSKE
VKEILADYDPNFMAMSLDEAYLNITKHLEERQNWPEDKRRYFIKMGSSVENDNPGKEVNK
LSEHERSISPLLFEESPSDVQPPGDPFQVNFEEQNNPQILQNSVVFGTSAQEVVKEIRFR
IEQKTTLTASAGIAPNTMLAKVCSDKNKPNGQYQILPNRQAVMDFIKDLPIRKVSGIGKV
TEKMLKALGIITCTELYQQRALLSLLFSETSWHYFLHISLGLGSTHLTRDGERKSMSVER
TFSEINKAEEQYSLCQELCSELAQDLQKERLKGRTVTIKLKNVNFEVKTRASTVSSVVST
AEEIFAIAKELLKTEIDADFPHPLRLRLMGVRISSFPNEEDRKHQQRSIIGFLQAGNQAL
SATECTLEKTDKDKFVKPLEMSHKKSFFDKKRSERKWSHQDTFKCEAVNKQSFQTSQPFQ
VLKKKMNENLEISENSDDCQILTCPVCFRAQGCISLEALNKHVDECLDGPSISENFKMFS
CSHVSATKVNKKENVPASSLCEKQDYEAHPKIKEISSVDCIALVDTIDNSSKAESIDALS
NKHSKEECSSLPSKSFNIEHCHQNSSSTVSLENEDVGSFRQEYRQPYLCEVKTGQALVCP
VCNVEQKTSDLTLFNVHVDVCLNKSFIQELRKDKFNPVNQPKESSRSTGSSSGVQKAVTR
TKRPGLMTKYSTSKKIKPNNPKHTLDIFFK
```

```
>POLK_MOUSE
MDNTKEKDNFKDDLLLRMGLNDNKAGMEGLDKEKINKIIMEATKGSRFYGNELKKEKQVN
QRIENMMQQKAQITSQQLRKAQLQVDKFAMELERNRNLNNTIVHVDMDAFYAAVEMRDNP
ELKDKPIAVGSMSMLATSNYHARRFGVRAAMPGFIAKRLCPQLIIVPPNFDKYRAVSKEV
KEILAEYDPNFMAMSLDEAYLNITQHLQERQDWPEDKRRYFIKMGNYLKIDTPRQEANEL
TEYERSISPLLFEDSPPDLQPQGSPFQLNSEEQNNPQIAQNSVVFGTSAEEVVKEIRFRI
EQKTTLTASAGIAPNTMLAKVCSDKNKPNGQYQILPSRSAVMDFIKDLPIRKVSGIGKVT
EKMLMALGIVTCTELYQQRALLSLLFSETSWHYFLHIALGLGSTDLARDGERKSMSVERT
FSEISKTEEQYSLCQELCAELAHDLQKEGLKGRTVTIKLKNVNFEVKTRASTVPAAISTA
EEIFAIAKELLRTEVNVGSPHPLRLRLMGVRMSTFSSEDDRKHQQRSIIGFLQAGNQALS
STGDSLDKTDKTELAKPLEMSHKKSFFDKKRSERISNCQDTSRCKTAGQQALQILEPSQA
LKKLSESFETSENSNDCQTFICPVCFREQEGVSLEAFNEHVDECLDGPSTSENSKISCYS
HASSADIGQKEDVHPSIPLCEKRGHENGEITLVDGVDLTGTEDRSLKAARMDTLENNRSK
EECPDIPDKSCPISLENETISTLSRQDSVQPCTDEVVTGRALVCPVCNLEQETSDLTLFN
IHVDICLNKGIIQELRNSEGNSVKQPKESSRSTDRLQKASGRTKRPGTKTKSSTLKKTKP
RDPRHTLDGFFK
```

## VESPA predictions
Run [VESPA](https://github.com/Rostlab/VESPA) on the individual sequences and store the results in the directory `vespa_predictions` as the files `POLK_HUMAN.csv` and `POLK_MOUSE.csv`.

## Sequence alignment
Use MMSeqs2 to align the sequences to each other:
```
mmseqs createdb POLK_HUMAN.fasta human
mmseqs createdb POLK_MOUSE.fasta mouse
mmseqs search human mouse result tmp
mmseqs result2msa human mouse result msa --msa-format-mode 6
head -n -1 msa > alignment_files/POLK.a3m
```

## Tell MutAlign the alignments it should process
Put the name of the identifiers you want to process (here: POLK) into the file `identifiers.txt`:
```
cat POLK > identifiers.txt
```
The identifiers should correspond to alignment file names. Respectively, the sequence identifiers inside of the alignment files should correspond to the file names of the VESPA outputs from above.

## Run MutAlign
```
python create_html.py
```
This will output HTML files in the directory `./html/`. If you want a different output directory, you can specify it as a parameter to the MutAlign script, e.g.:
```
python create_html.py /home/johndoe/public_html/
```
