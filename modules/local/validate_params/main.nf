nextflow.enable.dsl=2

// https://github.com/Plant-Food-Research-Open/assembly_qc
// GPL-3.0: https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE
def validateParams(params) {
    validateFastaTags(params)
}

def validateFastaTags(params) {
    def listOfFastaTuples   = params["target_assemblies"]

    if (isNotListOfLists(listOfFastaTuples, 2)) {
        error 'Error: target_assemblies must be a list of sublists, with each sublist containing 2 elements'
    }

    def fastaTags = listOfFastaTuples.collect { it[0] }

    fastaTags.each {
        if (!(it =~ /^\w+$/)) {
            error "Error: $it is not a valid tag in target_assemblies"
        }
    }

    if (fastaTags.size() != (fastaTags as Set).size()) {
        error "All the tags in target_assemblies should be unique"
    }
}

def isNotListOfLists(thisOne, subListSize) {
    return (!(thisOne instanceof List) || thisOne.isEmpty() || thisOne.any { !(it instanceof List) || it.size() != subListSize })
}