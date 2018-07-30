panels   = Channel.fromPath("${params.panel_dir}/*.bed")
foci     = Channel.fromPath("${params.foci_dir}/P*.vcf").filter( ~/^.*P\d+-T.*/ ).map{file -> [((file =~ /(P\d+)/)[0][1]): file] }

process score{
    tag { focus[focus.keySet()[0]].baseName.split('_')[0] + " => " + panel.baseName.replaceAll(/_locations_hg19/,'') }
    executor 'local'

    input:
    file panel from panels
    each focus from foci

    output:
    stdout into result

    shell: 
    '''
    # Define some variables
    DIR="!{focus[focus.keySet()[0]].getParent()}"
    PATIENT="!{focus.keySet()[0].split('_')[0]}"
    FOCUS="!{focus[focus.keySet()[0]].baseName}"
    PANEL="!{panel.baseName}"

    # Somatic
    bedtools subtract -header -a $DIR/$FOCUS.vcf -b $DIR/$PATIENT-N*.vcf > somatic.vcf
    
    # cfDNA
    !{params.home_dir}/../remove_chr.awk !{params.cfDNA_dir}/$PATIENT-called.vcf > cfDNA.vcf

    # Panel
    head -n!{params.panel_size} !{params.panel_dir}/$PANEL.bed | sort-bed - > panel.bed

    # cfDNA^Focus
    bedtools intersect -header -a cfDNA.vcf -b $DIR/$FOCUS.vcf > cfDNA_focus.vcf

    # cfDNA^Somatic
    bedtools intersect -header -a cfDNA.vcf -b somatic.vcf > cfDNA_somatic.vcf

    # Panel^cfDNA^Focus
    bedtools intersect -a panel.bed -b cfDNA_focus.vcf > panel_cfDNA_focus.vcf

    # Panel^cfDNA^Somatic
    bedtools intersect -a panel.bed -b cfDNA_somatic.vcf > panel_cfDNA_somatic.vcf

    # Append output to the final results file
    echo "Focus Name\tPanel Name\tPanel Size\tPanel^cfDNA^Somatic\tPanel^cfDNA^Focus\tcfDNA^Somatic\tcfDNA^Focus\tcfDNA\tSomatic\tFocus"
    echo -e "!{focus[focus.keySet()[0]].baseName.split('_')[0]}"'\t'\
            $PANEL'\t'\
            `grep -v "#" panel.bed               | wc -l`'\t'\
            `grep -v "#" panel_cfDNA_somatic.vcf | wc -l`'\t'\
            `grep -v "#" panel_cfDNA_focus.vcf   | wc -l`'\t'\
            `grep -v "#" cfDNA_somatic.vcf       | wc -l`'\t'\
            `grep -v "#" cfDNA_focus.vcf         | wc -l`'\t'\
            `grep -v "#" cfDNA.vcf               | wc -l`'\t'\
            `grep -v "#" somatic.vcf             | wc -l`'\t'\
            `grep -v "#" $DIR/$FOCUS.vcf         | wc -l`'\t'\

    '''
}

result.collectFile(name: "${params.home_dir}/panel_performances.tsv", keepHeader: true, newLine: false, sort:true)