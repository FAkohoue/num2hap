# **num2hap**

`num2hap` is a Python package for converting numeric genotype data into **HapMap diploid format**.  
It ensures **biologically meaningful haplotype representation** based on allele composition.

## **Function: `process_numeric_to_hapmap`**
This function converts **numeric genotype data** back to **HapMap diploid format** using SNP reference (`REF`) 
and alternate (`ALT`) alleles. The conversion is **parallelized for efficiency**, and users can 
choose between **different numeric encoding formats** (`012` or `-101`).

### **🔹 Novelty of This Function**
Unlike a simple numeric-to-text conversion, this function ensures that the **reconstructed HapMap 
genotypes retain biological meaning** based on each marker's alleles (`REF` and `ALT`).  
This ensures that the original allele composition is preserved, making the data compatible with genetic analyses.

## **📥 Input Requirements**
- The input file must be in **CSV format** and must contain the following five **mandatory columns** 
  **in this exact order**:
    1. `SNP`  (SNP marker name)
    2. `CHR`  (Chromosome number)
    3. `POS`  (Genomic position)
    4. `REF`  (Reference allele)
    5. `ALT`  (Alternate allele)
- **Genotype data** starts from **column 6 onward** and must be **numeric** (`0`, `1`, `2`, `-1`, or `-101`).
- The function assumes that **missing values** are represented as `-9`.

## **⚙ Parameters**
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

## **🚀 Output**
The function generates a **HapMap CSV file (`output_file`)**, which includes:
- **SNP metadata** (`SNP`, `CHR`, `POS`, `REF`, `ALT`) identical to the input.
- **Converted genotypic data**, where:
    - **Numeric values** are converted back to **biologically meaningful alleles** (`AA`, `AT`, `TT`).
    - **Missing values (`-9`)** are replaced with `"NA"`.
- The dataset maintains the **original SNP order**, ensuring compatibility with genetic tools.

### **🔬 Why Use This?**
This function is essential for:
✅ **GWAS (Genome-Wide Association Studies)**  
✅ **Genomic Selection & Breeding**  
✅ **Population Structure Analysis**  

## 🛠 Requirements**
- Python 3.7+
- pandas
- numpy

## **🔧 Installation**
### **Install from GitHub**
```bash
pip install git+https://github.com/FAkohoue/num2hap.git

```
## ** Usage **
```python
from num2hap import process_numeric_to_hapmap

process_numeric_to_hapmap(
    input_file="Input_numeric_genotype.csv",
    output_file="Output_hapmap_genotype.csv",
    format_type="012",  # Use "-101" if needed
    chunk_size=100,  # Adjust based on dataset size
    num_processes=10  # Optimize based on system resources
)
