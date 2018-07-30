panels   = Channel.fromPath("${params.panel_dir}/*.bed")
foci     = Channel.fromPath("${params.foci_dir}/P*.vcf").filter( ~/^.*P\d+-T.*/ ).map{file -> [((file =~ /(P\d+)/)[0][1]): file] }

process score{
    tag { focus[focus.keySet()[0]].baseName + " => " + panel.baseName.replaceAll(/_locations_hg19/,'') }
    echo true
    executor 'local'
    afterScript "sort -k3,3 -k1,1nr ${params.home_dir}/somatic_panel.tsv > ${params.home_dir}/somatic_panel_sorted.tsv"
    
    input:
    file panel from panels
    each focus from foci

    shell: 
    '''
    # Define some variables
    DIR="!{focus[focus.keySet()[0]].getParent()}"
    PATIENT="!{focus.keySet()[0]}"
    FOCUS="!{focus[focus.keySet()[0]].baseName}"
    PANEL="!{panel.baseName}"

    # Remove germline mutations from focus
    bedtools subtract -header -a $DIR/$FOCUS.vcf -b $DIR/$PATIENT-N*.vcf > somatic.vcf

    #Subset the panel
    head -n!{params.panel_size} !{params.panel_dir}/$PANEL.bed | sort-bed - > panel.bed

    # Intersect somatic mutations with the panel
    bedtools intersect -a somatic.vcf -b panel.bed > somatic_panel.vcf

    # Append output to the final results file
    echo -e `grep -v "#" somatic_panel.vcf | wc -l`'\t'`grep -v "#" somatic.vcf | wc -l`'\t'$FOCUS'\t'$PANEL >> !{params.home_dir}/somatic_panel.tsv
    '''
}
