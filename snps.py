#!/usr/bin/env python3

import subprocess
import os
import csv

'''
    This will process a .cram file from a genomic database and search for a specific variant (see line 57 for details)
    The output will print to stdout and then to a report csv file.

    IMPORTANT: The .cram file needs to be renamed to seq.cram and in the same working directory

    Usage: python snps.py seq.cram
'''

def run_command(command):
    """Run a shell command and return the output."""
    process = subprocess.run(command, shell=True, check=True, text=True, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return process.stdout.strip()

def main():
    # 1. Index the CRAM file
    print("Indexing CRAM file...")
    subprocess.run("samtools index seq.cram", 
                  shell=True, check=True)
    
    # 2. Extract the CYP1A2 region (region contains CYP1A2 gene on chromosome 15)
    print("Extracting CYP1A2 region...")
    subprocess.run("samtools view seq.cram chr15:74745879-74756607 -b > CYP1A2.bam", 
                  shell=True, check=True)

    # 3. Index the BAM file
    print("Indexing BAM file...")
    subprocess.run("samtools index CYP1A2.bam", shell=True, check=True)

    # 4. Use bcftools to call variants
    print("Calling variants...")
    cmd = ("bcftools mpileup -f ../../../GRCh38.d1.vd1.fa CYP1A2.bam | "
           "bcftools call -mv -Ov -o CYP1A2_variants.vcf.gz")
    subprocess.run(cmd, shell=True, check=True)

    # 5. Annotate the variants using SnpEff
    print("Annotating variants...")
    cmd = ("java -jar /Users/matthewkinder/snpEff/snpEff.jar -v GRCh38.99 "
           "CYP1A2_variants.vcf.gz > CYP1A2_annotated.vcf")
    subprocess.run(cmd, shell=True, check=True)

    # 6. Compress and index the annotated VCF file for bcftools compatibility
    print("Compressing and indexing annotated VCF file...")
    subprocess.run("bgzip -c CYP1A2_annotated.vcf > CYP1A2_annotated.vcf.gz", 
                  shell=True, check=True)
    subprocess.run("tabix -p vcf CYP1A2_annotated.vcf.gz", shell=True, check=True)

    # 7. Extract the specific SNPs of interest
    print("Extracting specific SNPs of interest...")
    
    # Known SNPs to analyze in json format
    # id: name of the varant
    # name: alternate name of the variant
    # position: location of the variant
    snps = [
        {"id": "rs762551", "name": "CYP1A2*1F", "position": "chr15:74749576"},
        {"id": "rs2069514", "name": "CYP1A2*1C", "position": "chr15:74754823"},
        {"id": "rs2472300", "name": "", "position": "chr15:74749234"},
        {"id": "rs2472304", "name": "", "position": "chr15:74751897"},
        {"id": "rs2470893", "name": "", "position": "chr15:74745879"}
    ]

    # Create CSV report file
    with open("snp_report.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["SNP", "Position", "Reference", "Alternative", "Genotype", "Annotation"])

        for snp in snps:
            snp_id = snp["id"]
            position = snp["position"]
            snp_name = snp["name"]
            
            description = f"{snp_id}"
            if snp_name:
                description += f" ({snp_name})"
                
            print(f"Checking {description} at position {position}...")
            
            # Run bcftools to extract the SNP information
            temp_file = f"temp_{snp_id}.txt"
            # Execute commands separately to handle the case where no variants are found
            bcftools_cmd = f"bcftools view CYP1A2_annotated.vcf.gz {position} > temp_bcftools_output.txt"
            subprocess.run(bcftools_cmd, shell=True, check=False)
            
            # Now grep from the temporary file, but don't check for errors
            # since grep will return exit code 1 if no lines match
            grep_cmd = f"grep -v '^#' temp_bcftools_output.txt > {temp_file}"
            subprocess.run(grep_cmd, shell=True, check=False)
            
            # Clean up the intermediate file
            if os.path.exists("temp_bcftools_output.txt"):
                os.remove("temp_bcftools_output.txt")
            
            # Check if any variants were found
            if os.path.getsize(temp_file) > 0:
                with open(temp_file, "r") as snp_file:
                    line = snp_file.readline().strip()
                    if line:
                        parts = line.split("\t")
                        if len(parts) >= 10:
                            genotype = parts[9].split(":")[0]
                            ref = parts[3]
                            alt = parts[4]
                            info = parts[7]
                            writer.writerow([snp_id, position, ref, alt, genotype, info])
                            print(f"  Found {snp_id} variant")
                        else:
                            writer.writerow([snp_id, position, "not found", "not found", "not found", "not found"])
                            print(f"  No variant found at {snp_id} position")
                    else:
                        writer.writerow([snp_id, position, "not found", "not found", "not found", "not found"])
                        print(f"  No variant found at {snp_id} position")
            else:
                writer.writerow([snp_id, position, "not found", "not found", "not found", "not found"])
                print(f"  No variant found at {snp_id} position")

    # 7. Clean up temporary files
    for snp in snps:
        temp_file = f"temp_{snp['id']}.txt"
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    # Also remove any other temporary files that might be left
    if os.path.exists("temp_bcftools_output.txt"):
        os.remove("temp_bcftools_output.txt")

    # 8. Convert genotype codes to actual alleles
    print("Processing report with actual alleles...")
    
    # Read the current CSV
    rows = []
    with open("snp_report.csv", "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) >= 5:  # Ensure there's a genotype column
                if row[4] == "0/0":
                    row[4] = "Homozygous Reference"
                elif row[4] == "0/1":
                    row[4] = "Heterozygous"
                elif row[4] == "1/1":
                    row[4] = "Homozygous Alternate"
            rows.append(row)
    
    # Write the updated CSV
    with open("snp_report.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    print("=== SNP Analysis Complete ===")
    print("Report saved to: snp_report.csv")
    print("")

if __name__ == "__main__":
    main()