"""
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
# Remove warning
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

# Load required libraries
import pandas as pd
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial

def process_row(row, format_type):
    """
    Convert numeric genotypes in a single row (SNP) to HapMap format.
    
    Args:
        row (pd.Series): A row from the DataFrame representing a SNP.
        format_type (str): The numeric encoding format ('012' or '-101').
    
    Returns:
        pd.Series: The row with converted genotype values.
    """
    ref_allele = row['REF']
    alt_allele = row['ALT']
    genotype_columns = row.index[5:]  # Columns after the first five are genotypes
    
    # Determine the conversion mapping based on the format type
    if format_type == "012":
        conversion_map = {
            0: ref_allele + ref_allele,
            1: ref_allele + alt_allele,
            2: alt_allele + alt_allele,
            -9: 'NA'
        }
    elif format_type == "-101":
        conversion_map = {
            -1: ref_allele + ref_allele,
            0: ref_allele + alt_allele,
            1: alt_allele + alt_allele,
            -9: 'NA'
        }
    else:
        raise ValueError("Invalid format_type. Choose either '012' or '-101'.")
    
    # Convert the genotype values using the mapping
    converted = row[genotype_columns].replace(conversion_map)
    # Ensure any values not in the conversion map are set to 'NA'
    converted = converted.apply(lambda x: 'NA' if x not in conversion_map.values() and x != 'NA' else x)
    
    # Update the row with converted genotypes
    row[genotype_columns] = converted
    return row

def process_chunk(chunk, format_type):
    """
    Process a chunk of SNPs by converting each row's genotype data.
    
    Args:
        chunk (pd.DataFrame): A subset of the DataFrame (rows = SNPs).
        format_type (str): The numeric encoding format.
    
    Returns:
        pd.DataFrame: Processed chunk with genotypes converted.
    """
    return chunk.apply(process_row, axis=1, format_type=format_type)

def process_numeric_to_hapmap(input_file, output_file, format_type="012", chunk_size=100, num_processes=10):
    """
    Convert numeric genotype data to HapMap format with parallel processing.
    
    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to save the output HapMap CSV file.
        format_type (str): Numeric encoding format ('012' or '-101').
        chunk_size (int): Number of SNPs processed per batch.
        num_processes (int): Number of parallel processes.
    """
    # Load the input data
    print(f"Loading numeric genotype file: {input_file}")
    genotype_data = pd.read_csv(input_file, sep=",")
    
    # Validate the required columns
    required_columns = ["SNP", "CHR", "POS", "REF", "ALT"]
    if genotype_data.columns.tolist()[:5] != required_columns:
        raise ValueError(f"Input file must have the first five columns as: {required_columns}")
        
    # Validate format_type and check for invalid genotype values
    allowed_values = {
        "012": {0, 1, 2, -9},
        "-101": {-1, 0, 1, -9}
    }
    if format_type not in allowed_values:
        raise ValueError(f"Invalid format_type: {format_type}. Choose '012' or '-101'.")
    
    genotype_columns = genotype_data.columns[5:]
    invalid_snps = []
    
    for idx, row in genotype_data.iterrows():
        snp_name = row['SNP']
        genotypes = row[genotype_columns]
        
        # Check for invalid values
        unique_vals = set(pd.to_numeric(genotypes, errors='coerce').unique())
        invalid_vals = {v for v in unique_vals if v not in allowed_values[format_type] and not np.isnan(v)}
        if invalid_vals:
            invalid_snps.append((snp_name, invalid_vals))
    
    if invalid_snps:
        error_msg = (
            f"🚨 Format mismatch! Input data contains values inconsistent with `format_type='{format_type}'.\n"
            "Affected SNPs and invalid values:\n" +
            "\n".join([f"- {snp}: {', '.join(map(str, vals))}" for snp, vals in invalid_snps]) +
            "\n\n🛠️ Fix: Ensure your data matches the specified format or adjust `format_type`."
        )
        raise ValueError(error_msg)
        
    # Split the data into chunks of rows for parallel processing
    chunks = [genotype_data.iloc[i:i + chunk_size] for i in range(0, len(genotype_data), chunk_size)]
    
    # Process each chunk in parallel
    print("Starting genotype conversion...")
    with Pool(num_processes) as pool:
        process_chunk_partial = partial(process_chunk, format_type=format_type)
        processed_chunks = list(tqdm(pool.imap(process_chunk_partial, chunks), total=len(chunks), desc="Processing chunks"))
    # Combine all processed chunks
    converted_data = pd.concat(processed_chunks)
    
    # Save the converted data to the output file
    print(f"Saving HapMap-formatted data to {output_file}")
    converted_data.to_csv(output_file, index=False)
    
    print("Conversion completed successfully!")