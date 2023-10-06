1. Rename perform_edta_annotation to FASTA_PERFORM_EDTA
2. Extract subworkflows
3. STAR ignores softmasking and, thus, should be fed the unmasked genome so that masking and mapping can run in parallel.