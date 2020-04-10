# Agnostic Search

Summary

Tandem mass spectrometry enables high throughput peptide identification in complex biological specimens. In conventional method, peptide identification relies on database search which is limited by identifying only small number of post translational modifications (less than 3) and such method is unable to identify undefined PTMs. Here we developed SAMPEI, a searching method leveraging high quality query spectra within the same or different dataset to assign target spectra with peptide sequence and undefined modification (mass shift). Prior to SAMPEI, we utilized database search (X!tandem) to assign spectra with peptide sequences in each sample. Only spectra with unique modification setting and with the highest peptide identification confidence (lowest e value) were selected as queries. SAMPEI would then perform a series of orthogonal measures to evaluate the similarity between all unassigned spectra and query spectra within predefined mass difference window (default= +/- 200 Dalton). It firstly aligns discrete m/z ranges within unassigned spectra to the query spectra. The proportion of matched MS2 ion intensity from query over the total MS2 intensity defined as matched query intensity is used to pre-select candidate spectra. Then, two additional measures assess the quality of the assignment against query peptide sequence were determined to evaluate the goodness of the matching. Specifically, the proportion of MS2 intensity of target spectrum matched to the theoretical m/z of the query peptide sequence over the total MS2 intensity in the target scan is one of the measures. Finally, the proportion of largest consecutive b/y ion missing over the length of the peptide sequence defined as largest gap percentage is the last measure. We experimentally determined the optimal cutoff of these measures with two alkylation treated human cell lines. The output highlights all target scans neglected by the conventional search and assigned mass shift to the position with highest confidence. The cutoff of the measures can be adjusted and optimized by user’s own specimen. 


Clone the repository
```shell
git clone https://github.com/FenyoLab/AgnosticSearch.git
```

Navigate to the cloned repository
```shell
cd AgnosticSearch
```

If you have not installed Command Line Tools run this command
```shell
xcode-select --install
```

Create the virtual environment (requires make which comes with Command Line Tools)
```shell
make venv
```

Activate the virtual environment (Do this every time you open a new terminal)
```shell
source venv/bin/activate
```

Run agnostic search
```shell
python -m src.agnostic_search.cli <query_mgf_path> <target_mgf_path> <id_path>
```

Run this for help on the command line arguments for setting the output dir, errors etc
```shell
python -m src.agnostic_search.cli --help
```
