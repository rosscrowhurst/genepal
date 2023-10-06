1. Rename perform_edta_annotation to FASTA_PERFORM_EDTA
2. Fix braker 3 version reporting problem
3. Extract subworkflows
4. STAR ignores softmasking and, thus, should be fed the unmasked genome so that masking and mapping can run in parallel.