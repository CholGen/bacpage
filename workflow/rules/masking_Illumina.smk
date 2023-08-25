rule mask_consensus:
    message: "Mask wholly recombinant regions of the consensus sequence of {wildcards.sample} using the mask provided in {input.recombinant_mask}"
    input:
        consensus=rules.call_consensus.output.consensus_sequence,
        recombinant_mask=config["recombinant_mask"]
    output:
        first_tmp_consensus=temp( "intermediates/illumina/consensus/{sample}.tmp1.fasta" ),
        second_tmp_consensus=temp( "intermediates/illumina/consensus/{sample}.tmp2.fasta" ),
        masked_consensus="results/consensus_sequences/{sample}.masked.fasta"
    shell:
        """
        HEADER=$(cut -f1 {input.recombinant_mask} | uniq | head -n1) && \
        sed "1s/.*/>${{HEADER}}/" {input.consensus} > {output.first_tmp_consensus} && \
        bedtools maskfasta \
            -fi {output.first_tmp_consensus} \
            -bed {input.recombinant_mask} \
            -fo {output.second_tmp_consensus} && \
        sed "1s/.*/>{wildcards.sample}/" {output.second_tmp_consensus} > {output.masked_consensus}
    """
