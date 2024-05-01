## Usage
### 1. Initialize the 2D track class:

```
track_2d = Characteristic2dWithLimit(
    res_fold="hic_result",
    loader_type="hard",
    hic_path=file_path_mcool)
```
### 2. Initialize the DNA sequence class:

```
dna_seq = DNASequenceBase(
    fai_path="helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai",
    fa_path="helper/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
```
### 3. Initialize the 1D track class:

```
tracks_1d = CharacteristicBigWigCSR(
    folder_res="/home/ojpochemy/dnaloader/resfold_hard",
    size_of_bw=24,
    type_of_loader="hard",
    path="filepaths.txt")
```
### 4. Create the DNA dataset:

```
dna_dataset = _DnaDataset(
    window_size=1_000,
    dna_seq=dna_seq,
    char_list=[tracks_1d, track_2d],
    seed=59)
```
### 5. Load the data using DataLoader:

```
datal = DataLoader(
    dna_dataset,
    num_workers=3,
    batch_size=10)
```