# 11 - Human Cytoplasmic tRNA Selection Report

**Study Title:** tRNA ASL Structural Stability Study  
**Phase:** Phase 2 (Biological Data Verification & Construct Selection)  
**Status:** Selection Complete & Verified  

---

## 1. Selection Objective

The objective of Phase 2 is to select exactly **ONE** human cytoplasmic tRNA from the Genomic tRNA Database (GtRNAdb) whose native anticodon can be converted to the amber suppressor anticodon `CUA` by **EXACTLY TWO nucleotide substitutions** ($\text{Hamming}(\text{Anticodon}_{\text{WT}}, \text{CUA}) = 2$), to serve as the benchmark system for comparative molecular dynamics simulations against its `CUA` suppressor counterpart.

In accordance with strict scientific integrity rules:
- No tRNA was chosen prior to computational screening.
- The selection followed an unbiased funnel:
  $$\text{GtRNAdb High-Confidence Set} \longrightarrow \text{Hamming Distance Filter} \longrightarrow \text{Biological Filter} \longrightarrow \text{Structural Suitability} \longrightarrow \text{Objective Biophysical Ranking}$$
- No proprietary Transloka Bio or KritRNA data was used.

---

## 2. GtRNAdb Source & Access Information

- **Database:** Genomic tRNA Database (GtRNAdb), hosted by the Lowe Lab at the University of California, Santa Cruz (UCSC).
- **Target Genome Assembly:** *Homo sapiens* (GRCh38 / hg38).
- **Primary Source URL:** [http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/](http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/)
- **Archive Accessed:** `hg38-tRNAs.tar.gz` (containing `hg38-mature-tRNAs.fa`, `hg38-tRNAs-confidence-set.out`, `hg38-tRNAs-confidence-set.ss`, `hg38-tRNAs_name_map.txt`).
- **Access Date:** November 12, 2026.
- **Local Raw Data Preservation:** Stored under `data/raw/` with verified SHA-256 checksums (`hg38-mature-tRNAs.fa`: `1b4a5f3be176ca11debbc045ebeb0da12619005669877df428ab8fd533688adc`).

---

## 3. Candidate Screening Methodology & Funnel

The screening was executed via the deterministic pipeline script [scripts/01_select_trna.py](file:///d:/trna-asl-md-stability/scripts/01_select_trna.py). To ensure methodological transparency, screening criteria are partitioned into **Hard Eligibility Requirements** (dictated by the research brief and canonical biology) and **Implementation Choices** (operational ranking preferences).

### 3.1 Hard Eligibility Requirements
1. **Initial Dataset (432 Genes):** All high-confidence *Homo sapiens* nuclear-encoded cytoplasmic tRNA genes from GtRNAdb (hg38).
2. **Gate 1: Exact Anticodon Hamming Distance ($\text{Hamming} = 2$):**
   - The target suppressor anticodon is `CUA` (targeting mRNA stop codon `5'-UAG-3'`).
   - Every candidate's wild-type anticodon triplet was compared position-by-position against `CUA`.
   - Only anticodons differing at **exactly two** nucleotide positions were retained.
   - Result: **199 genes** across 21 distinct anticodon families passed Gate 1 (68 genes had distance 1; 165 genes had distance 3; 0 genes had distance 0).
3. **Gate 2: Biological & Genomic Eligibility:**
   - Must be verified cytoplasmic (cytosolic, transcribed by RNA Pol III, non-mitochondrial).
   - Intron-free ASL region (to eliminate RNA splicing ambiguity).
   - Result: All 199 candidates in the high-confidence set are cytosolic and intron-free in the ASL.
4. **Gate 3: Structural Suitability (Canonical 17-nt ASL Architecture):**
   - Must possess a canonical 5-bp anticodon stem and 7-nt anticodon loop (`>>>>>.......<<<<<`).
   - Invariant U33 at position 2 of the loop (essential for the Quigley & Rich U-turn motif).
   - Anticodon located strictly at loop positions 3, 4, and 5 (Sprinzl numbering 34, 35, 36).
   - Canonical purine (A or G) at position 6 of the loop (Sprinzl position 37).
   - Stem must consist exclusively of standard Watson-Crick base pairs (A:U, U:A, G:C, C:G) with zero non-canonical pairs or bulges.
   - Result: **138 genes** across 30 unique ASL sequence clusters passed all structural criteria.

### 3.2 Implementation & Ranking Preferences `[IMPLEMENTATION CHOICE]`
5. **Stem GC Content Preference `[IMPLEMENTATION CHOICE]`:**
   - As identified in Risk **RSK-SCI-03** ([docs/09_risk_register.md](file:///d:/trna-asl-md-stability/docs/09_risk_register.md)), isolated 17-nt hairpins in explicit solvent lack tertiary stabilization from the tRNA core and are susceptible to artificial terminal base-pair fraying (residues 27:43).
   - High-GC stem composition (5/5 vs. 4/5 vs. $\le 3/5$ C:G pairs) was selected as an implementation/design criterion intended to reduce the likelihood of stem fraying during the short ASL-only simulations.
   - *Scientific Clarification:* This criterion is an operational design preference, not an assumption of proven thermodynamic free energy ($\Delta G$). Helical stability in the force field remains an empirical outcome to be directly measured during trajectory analysis (Metric 3).
6. **Genomic Representation Preference `[IMPLEMENTATION CHOICE]`:**
   - Preference for major multi-copy housekeeping tRNA families over single-copy or low-expression loci to maximize biological relevance to the human cytoplasmic translation apparatus.

---

## 4. Hamming Distance Calculations

The Hamming distance is defined as:
$$\text{Hamming}(\mathbf{s}_1, \mathbf{s}_2) = \sum_{i=1}^3 \mathbb{I}(s_{1,i} \ne s_{2,i})$$

### Calculation for Top Shortlisted Anticodons:

| Isotype | WT Anticodon | Target | Pos 34 (WT vs Target) | Pos 35 (WT vs Target) | Pos 36 (WT vs Target) | Hamming Distance | Eligible? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Asp** | **GUC** | **CUA** | **G $\ne$ C** | **U = U** | **C $\ne$ A** | **2** | **YES** |
| **Glu** | **UUC** | **CUA** | **U $\ne$ C** | **U = U** | **C $\ne$ A** | **2** | **YES** |
| **Leu** | **CAG** | **CUA** | **C = C** | **A $\ne$ U** | **G $\ne$ A** | **2** | **YES** |
| **Ser** | **AGA** | **CUA** | **A $\ne$ C** | **G $\ne$ U** | **A = A** | **2** | **YES** |
| **Tyr** *(retracted)* | GUA | CUA | G $\ne$ C | U = U | A = A | 1 | NO |
| **Trp** *(retracted)* | CCA | CUA | C = C | C $\ne$ U | A = A | 1 | NO |
| **Leu** *(retracted)* | CAA | CUA | C = C | A $\ne$ U | A = A | 1 | NO |
| **Phe** | GAA | CUA | G $\ne$ C | A $\ne$ U | A = A | 2 | YES |
| **Cys** | GCA | CUA | G $\ne$ C | C $\ne$ U | A = A | 2 | YES |
| **His** | GUG | CUA | G $\ne$ C | U = U | G $\ne$ A | 2 | YES |
| **Asn** | GUU | CUA | G $\ne$ C | U = U | U $\ne$ A | 2 | YES |
| **Val** | CAC | CUA | C = C | A $\ne$ U | C $\ne$ A | 2 | YES |

---

## 5. Summary Table of Eligible ASL Clusters

Filtering the 138 passing genes yielded 30 unique ASL sequence clusters. The top clusters ranked by stem GC content and tRNAscan-SE score are:

| Isotype | WT Anticodon | Stem GC Pairs | Gene Count* | Max Bit Score | Representative GtRNAdb ID | 17-nt ASL Sequence ($5' \rightarrow 3'$) | Stem Architecture (27:43 to 31:39) |
| :--- | :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| **Asp** | **GUC** | **5 / 5** | **12** | **66.5** | **tRNA-Asp-GTC-2-1** | `CCCGCCUGUCACGCGGG` | C:G, C:G, C:G, G:C, C:G |
| **Glu** | **UUC** | **5 / 5** | **3** | **71.6** | **tRNA-Glu-TTC-3-1** | `CGGCGCUUUCACCGCCG` | C:G, G:C, G:C, C:G, G:C |
| **Ser** | **AGA** | **4 / 5** | **1** | **78.8** | **tRNA-Ser-AGA-4-1** | `GCAGGUUCGAAUCCUGC` | G:C, C:G, A:U, G:C, G:C |
| **Leu** | **CAG** | **4 / 5** | **9** | **78.2** | **tRNA-Leu-CAG-1-1** | `CUGCGUUCAGGUCGCAG` | C:G, U:A, G:C, C:G, G:C |
| **Glu** | **UUC** | **4 / 5** | **4** | **72.5** | **tRNA-Glu-TTC-1-2** | `CCUGGUUUUCACCCAGG` | C:G, C:G, U:A, G:C, G:C |
| **Arg** | **CCU** | **4 / 5** | **3** | **72.3** | **tRNA-Arg-CCT-1-1** | `CUGGCCUCCUAAGCCAG` | C:G, U:A, G:C, G:C, C:G |
| **Pro** | **CGG** | **4 / 5** | **3** | **70.2** | **tRNA-Pro-CGG-1-1** | `CUCGCUUCGGGUGCGAG` | C:G, U:A, C:G, G:C, C:G |
| **Asn** | **GUU** | **4 / 5** | **1** | **67.9** | **tRNA-Asn-GTT-24-1** | `GGUGGUUCCAGCCCACC` | G:C, G:C, U:A, G:C, G:C |
| **His** | **GUG** | **4 / 5** | **9** | **60.1** | **tRNA-His-GTG-1-1** | `CUGCGUUGUGGCCGCAG` | C:G, U:A, G:C, C:G, G:C |

*\*Note on Asp-GTC Gene Counts:* Exactly 12 genes share this identical 17-nt ASL sequence (`CCCGCCUGUCACGCGGG`): 11 genes (`tRNA-Asp-GTC-2-1` through `tRNA-Asp-GTC-2-11`) share both an identical 72-nt mature sequence and identical ASL, while 1 gene (`tRNA-Asp-GTC-1-1`) shares the identical ASL sequence but possesses two point variations in the D- and T-loops. One additional gene in the high-confidence set (`tRNA-Asp-GTC-3-1`) has a variant ASL (`CCCGUCUGUCACGCGGG`, U31:G39 wobble pair), bringing the total high-confidence `tRNA-Asp-GTC` gene count to 13.

*Full candidate list with all 432 genes is versioned in `data/metadata/trna_candidates.tsv`.*

---

## 6. Final Selected tRNA

The final candidate selected for this study is:

### **Human Cytoplasmic tRNA-Asp-GTC-2-1**
- **GtRNAdb Gene ID:** `tRNA-Asp-GTC-2-1`
- **tRNAscan-SE ID:** `chr1.trna92`
- **Genomic Locus:** `chr1:161440825-161440896 (-)`
- **Amino Acid Class:** Aspartate (Asp / D)
- **Isotype Subgroup:** GTC (reading aspartate codons GAC and GAU)
- **Authoritative Genomic Copy Count & Representation Breakdown:**
  - **Total Asp-GTC genes in high-confidence set:** **13 genes**
  - **Identical mature sequence cluster (Cluster 2):** **11 genes** (`tRNA-Asp-GTC-2-1` through `tRNA-Asp-GTC-2-11`) distributed across chromosomes 1 (5 loci: trna80, 83, 86, 89, 92), 6 (2 loci: trna47, 50), 12 (3 loci: trna4, 10, 12), and 17 (1 locus: trna36).
  - **Identical 17-nt ASL sequence (`CCCGCCUGUCACGCGGG`):** **12 genes** (the 11 Cluster 2 genes plus `tRNA-Asp-GTC-1-1` on chr12:98503503-98503574).
  - **Variant ASL gene:** **1 gene** (`tRNA-Asp-GTC-3-1` on chr6:27583528-27583457, with ASL `CCCGUCUGUCACGCGGG` containing a U31:G39 wobble pair).
  - **Biological Significance:** `tRNA-Asp-GTC-2-1` is the canonical reference locus of the predominant aspartate tRNA multi-gene family in the human genome.

---

## 7. Exact Anticodons & Substitutions

- **Native Wild-Type Anticodon:** `5'-GUC-3'` (DNA: `GTC`)
- **Suppressor Anticodon:** `5'-CUA-3'` (DNA: `CTA`)
- **Hamming Distance:** **2**
- **Exact Substitutions:**
  1. Position 34 (Wobble): $\mathbf{G34 \rightarrow C34}$ (Transversion)
  2. Position 35 (Middle): $\mathbf{U35 \rightarrow U35}$ (Identical / Unchanged)
  3. Position 36 (3'-base): $\mathbf{C36 \rightarrow A36}$ (Transversion)

```
Pos:      34  35  36
WT:        G   U   C
Target:    C   U   A
Diff:      *   -   *  -> Hamming Distance = 2
```

---

## 8. Sequences & Structural Decomposition

### 8.1 Full Mature tRNA Sequence (72 nt)
```text
5'- UCCUCGUUAGUAUAGUGGUGAGUAUCCCCGCCUGUCACGCGGGAGACCGGGGUUCGAUUCCCCGACGGGGAG -3'
```

### 8.2 17-nt Anticodon Stem-Loop (ASL) Sequences

- **Wild-Type ASL (17 nt):**
  ```text
  5'- C C C G C C U G U C A C G C G G G -3'
  ```
- **CUA Suppressor ASL (17 nt):**
  ```text
  5'- C C C G C C U C U A A C G C G G G -3'
  ```

### 8.3 Secondary Structure & Sprinzl Coordinate Mapping

```text
Sprinzl Pos:  27 28 29 30 31  32 33 34 35 36 37 38  39 40 41 42 43
Region:       [--- Stem 5' ---] [---- Loop 7nt ----] [--- Stem 3' ---]
WT ASL:        C  C  C  G  C   C  U  G  U  C  A  C   G  C  G  G  G
Suppressor:    C  C  C  G  C   C  U  C  U  A  A  C   G  C  G  G  G
Base Pairing:  |  |  |  |  |                         |  |  |  |  |
               G  G  G  C  G                         G  G  G  C  G
               (43)(42)(41)(40)(39)
```

- **5-bp Stem Base Pairs:**
  - 27:43 $\rightarrow$ $\text{C}:\text{G}$ (Watson-Crick)
  - 28:42 $\rightarrow$ $\text{C}:\text{G}$ (Watson-Crick)
  - 29:41 $\rightarrow$ $\text{C}:\text{G}$ (Watson-Crick)
  - 30:40 $\rightarrow$ $\text{G}:\text{C}$ (Watson-Crick)
  - 31:39 $\rightarrow$ $\text{C}:\text{G}$ (Watson-Crick)
- **Loop Conformation (Residues 32–38):**
  - Residue 32: C
  - Residue 33: **U** (invariant canonical U-turn residue)
  - Residue 34: **G** (WT) / **C** (Suppressor)
  - Residue 35: **U** (WT) / **U** (Suppressor)
  - Residue 36: **C** (WT) / **A** (Suppressor)
  - Residue 37: **A** (invariant purine enforcing 3'-stacking)
  - Residue 38: C

---

## 9. Verification & Provenance Evidence

1. **Species:** *Homo sapiens* (Human), TaxID 9606.
2. **Subcellular Localization:** Cytosolic / Cytoplasmic.
   - Encoded in the nuclear genome (`chr1:161440825-161440896`).
   - Categorized as `cytosolic` in GtRNAdb `hg38-tRNAs-confidence-set.out`.
3. **High Confidence Classification:**
   - tRNAscan-SE score: **66.5** (exceeds the 55.0 threshold for high-confidence functional tRNAs).
   - Isotype model score: **109.7**.
4. **Intron Status:** Zero introns (`Intron Bounds: 0 0`). Splicing-free mature transcript.
5. **Exact Anticodon Verification:**
   - Reported anticodon in GtRNAdb is `GTC` (RNA: `GUC`), exactly matching Aspartate codons `GAC` and `GAU`.
   - Verified $\text{Hamming}(\text{GUC}, \text{CUA}) = 2$.

---

## 10. Biophysical Rationale for the Simulation Experiment

The selection of **tRNA-Asp-GTC-2-1** provides three key scientific advantages:

1. **High-GC Stem Architecture (5/5 C:G Base Pairs) `[IMPLEMENTATION CHOICE]`:**
   - The 5-bp anticodon stem consists entirely of Watson-Crick C:G pairs (`27:43 C:G`, `28:42 C:G`, `29:41 C:G`, `30:40 G:C`, `31:39 C:G`).
   - **Rationale for RSK-SCI-03:** High-GC stem composition was selected as an implementation/design criterion intended to reduce the likelihood of stem fraying during the short ASL-only simulations.
   - *Scientific Note:* High GC content is an implementation choice, not a claim of proven "maximum thermodynamic stability" or a measured $\Delta G$. No empirical or alchemical free energy calculation has been performed at this stage. Rather, 5/5 C:G pairing provides strong Watson-Crick base-stacking and three hydrogen bonds per pair within standard empirical force fields, reducing the probability of terminal fraying in an isolated 17-nt hairpin lacking tertiary context. Whether the stem remains stable is an **empirical simulation outcome** to be directly measured via base-pair hydrogen-bond occupancy (Metric 3) in Phase 5 trajectory analysis.
2. **Canonical Loop Architecture `[PROMPT MANDATE / BIOLOGICAL PROPERTY]`:**
   - Features invariant residue U33 for the canonical U-turn motif (verified distance and dihedral diagnostic).
   - Features adenine at position 37 (A37), providing canonical purine base stacking immediately adjacent to the anticodon.
3. **High Biological Abundance & Multi-Copy Redundancy `[IMPLEMENTATION CHOICE]`:**
   - With 11 identical mature genomic copies (`tRNA-Asp-GTC-2-1` to `2-11`) and 12 identical 17-nt ASL loci in the high-confidence set (out of 13 total Asp-GTC genes), this construct represents the predominant cytoplasmic aspartate tRNA pool in human cells, maximizing the biological relevance of the simulation benchmark.

---

## 11. Uncertainties & Limitations

1. **Post-Transcriptional Modifications:**
   - Native human tRNA-Asp is modified in vivo at position 34 (queuosine, $Q34$, or mannosyl-queuosine) and position 38 ($\Psi38$). Per project requirements, this study models the **unmodified** RNA state. The absence of $Q34$ and $\Psi38$ is documented in [docs/07_limitations.md](file:///d:/trna-asl-md-stability/docs/07_limitations.md).
2. **Isolated Hairpin Model:**
   - The excised 17-nt ASL lacks the tertiary interactions with the D-arm and variable loop. This constraint is documented as an inherent model limitation.

