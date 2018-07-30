samtools quickcheck -v *.bam > bad_bams   && echo 'All ok' || echo 'Some files failed check, see bad_bams'
