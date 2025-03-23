"""
**num2hap**

`num2hap` is a Python package for converting numeric genotype data into hapmap diploid format. It ensures meaningful haplotype representation based on allele composition.

Function: process_numeric_to_hapmap
---------------------------------
This function converts numeric genotype data back to **HapMap diploid format** using SNP reference (`REF`) 
and alternate (`ALT`) alleles. The conversion is parallelized for efficiency, and the user can 
choose between **different numeric encoding formats** (012 or -101).

**Novelty of this function:**
Unlike a simple numeric-to-text conversion, this function ensures that the **reconstructed HapMap 
genotypes retain biological meaning** based on each marker's alleles (`REF` and `ALT`). This ensures 
that the original allele composition is preserved, making the data compatible with genetic analyses.

Input:
------
- The input file must be in **CSV format** and must contain the following five **mandatory columns** 
  **in this exact order**:
    1. `SNP`  (SNP marker name)
    2. `CHR`  (Chromosome number)
    3. `POS`  (Genomic position)
    4. `REF`  (Reference allele)
    5. `ALT`  (Alternate allele)
- **Genotype data** starts from **column 6 onward** and must be numeric (`0`, `1`, `2`, `-1`, or `-101`).
- The function assumes that missing values are represented as `-9`.

Parameters:
-----------
- input_file (str): Path to the input **numeric genotype** file in CSV format.
- output_file (str): Path to save the converted **HapMap genotype** file.
- format_type (str, optional): The numeric format used in the input file. Options:
    * `"012"`:
        - Homozygous reference (0) → `AA`
        - Heterozygous (1) → `AT`
        - Homozygous alternate (2) → `TT`
    * `"-101"`:
        - Homozygous reference (-1) → `AA`
        - Heterozygous (0) → `AT`
        - Homozygous alternate (1) → `TT`
    * `"-9"`: Missing values
- chunk_size (int, optional): Number of SNPs processed per batch. Default is **100**.
- num_processes (int, optional): Number of parallel CPU cores to use. Default is **10**.
Raises:
-------
- **ValueError**: If the input file does not contain the required columns in the correct order.
- **ValueError**: If `format_type` is not `"012"` or `"-101"`.

Value:
------
The function generates a new **HapMap CSV file (`output_file`)** containing:
- The same **SNP metadata** as the input (`SNP`, `CHR`, `POS`, `REF`, `ALT`).
- **Converted genotypic data**, where:
    - **Numeric values** are converted back to **biologically meaningful alleles** (`AA`, `AT`, `TT`).
    - **Missing values (`-9`)** are replaced with `"NA"`.
- The final dataset maintains the **original marker order**.

This function is essential for **GWAS, genomic selection, and population structure analyses**, where 
data must be in HapMap format for compatibility with genetic tools.

"""
## Usage
process_numeric_to_hapmap(
        input_file="/path_to_/Input_file.csv",
        output_file="/path_to_/Output_file.csv",
        format_type="012",  # Change to "-101" if needed
        chunk_size=100, # Change if needed
        num_processes=10 # Change if needed
    )
## Installation

To install from GitHub:
```bash
pip install git+https://github.com/FAkohoue/num2hap.git
