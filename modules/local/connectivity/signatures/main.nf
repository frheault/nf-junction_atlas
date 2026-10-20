process GENERATE_JUNCTION_SIGNATURES {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'scil.usherbrooke.ca/containers/scilus_2.0.2.sif':
        'scilus/scilus:2.0.2' }"

    input:
    tuple val(meta), path(trk), path(labels), path(wm), path(nufo), path(signatures)

    output:
    tuple val(meta), path("*__junction_labels.nii.gz"), emit: junction_labels, optional: true
    tuple val(meta), path("*.txt"), emit: junction_singatures, optional: true
    path ("split/")               , emit: split

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    decompose_simple.py ${trk} ${labels} split/

    generate_junctions.py ${signatures} split/ \
        ${wm} ${prefix}__junction_labels.nii.gz
    """
}
