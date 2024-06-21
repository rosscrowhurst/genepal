def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}

def validateFastqMetadata(metas, fqs, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check if multiple runs of a sample have the same target assemblies
    if ( metas.collect { meta -> meta.target_assemblies }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must target same assemblies"
    }

    // Check if multiple runs of a sample have the same endedness
    if ( metas.collect { meta -> meta.single_end }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must have same endedness"
    }

    [ metas.first(), fqs ]
}


def validateBamMetadata(metas, bams, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check that when the first file is bam then the second file is absent
    bams.findAll { files ->
        files.first().extension == 'bam' && files.size() != 1
    }
    .each { error "Sample ${metas.first().id} contains both bam and fastq pairs. When a bam file is provided as file_1, a fastq for file_2 is not permitted" }

    // Check that a bam file only targets a single assembly
    bams.eachWithIndex { files, index ->
        if ( files.first().extension == 'bam' && metas[index].target_assemblies.size() > 1 ) {
            error "BAM file for sample ${metas.first().id} can only target one assembly: ${metas[index].target_assemblies}"
        }
    }

    metas.every { it.target_assemblies == metas.first().target_assemblies }
    ? [ [ metas.first(), bams.flatten() ] ]
    : metas.withIndex().collect { meta, index -> [ meta, bams[index].flatten() ] }
}
