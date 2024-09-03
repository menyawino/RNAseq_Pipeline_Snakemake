# RNA-seq Pipeline

**Folder Structure**

- {**date**}{**project_name**}
    - **samples (needs to be added manually)**
        
        
    - **analysis (created and managed by the workflow)**
        
        
    - **results (created and managed by the workflow)**
        
        
    - **workflow (needs to be added manually)**
        
        
    
    📝 Metadata
    
    📝 Snakefile
    
    📝 README.md


## Note for pseudoalign multiqc
MultiQC parses the standard out from Kallisto, not any of its output files (abundance.h5, abundance.tsv, and run_info.json). As such, you must capture the Kallisto stdout to a file when running to use the MultiQC module.