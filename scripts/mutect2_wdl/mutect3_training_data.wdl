version 1.0

import "https://raw.githubusercontent.com/gatk-workflows/gatk4-somatic-snvs-indels/2.6.0/mutect2.wdl" as m2

workflow Mutect3TrainingData {
    input {
        File? intervals
        File? masks
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        File tumor_bam
        File tumor_bai
        File? normal_bam
        File? normal_bai
        File? pon
        File? gnomad
        File? variants_for_contamination
        String ref_downsample
        Boolean? run_orientation_bias_mixture_model_filter
        File? realignment_index_bundle
        String? realignment_extra_args
        String? m2_extra_args
        String? m2_extra_filtering_args
        File truth_vcf
        File truth_vcf_idx
        Boolean? make_bamout

        # runtime
        String gatk_docker
        File? gatk_override
        Int? preemptible
        Int? max_retries
    }

    call m2.Mutect2 {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            tumor_reads = tumor_bam,
            tumor_reads_index = tumor_bai,
            normal_reads = normal_bam,
            normal_reads_index = normal_bai,
            intervals = intervals,
            pon = pon,
            gnomad = gnomad,
            variants_for_contamination = variants_for_contamination,
            run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
            realignment_index_bundle = realignment_index_bundle,
            realignment_extra_args = realignment_extra_args,
            preemptible = preemptible,
            max_retries = max_retries,
            m2_extra_args = select_first([m2_extra_args, ""]) + " --training-data-mode --training-data-mode-ref-downsample " + ref_downsample,
            m2_extra_filtering_args = m2_extra_filtering_args,
            make_bamout = make_bamout,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call Concordance  {
        input:
            intervals = intervals,
            masks = masks,
            truth_vcf = truth_vcf,
            truth_vcf_idx = truth_vcf_idx,
            eval_vcf = Mutect2.filtered_vcf,
            eval_vcf_idx = Mutect2.filtered_vcf_idx,
            preemptible = preemptible,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call MakeTable {
        input:
            tpfp = Concordance.tpfp,
            tpfp_idx = Concordance.tpfp_idx,
            ftnfn = Concordance.ftnfn,
            ftnfn_idx = Concordance.ftnfn_idx,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            preemptible = preemptible
    }


    output {
        File table = MakeTable.table
    }
}

task Concordance {
    input {
        File? intervals
        File? masks
        File truth_vcf
        File truth_vcf_idx
        File eval_vcf
        File eval_vcf_idx

        File? gatk_override

        # runtime
        String gatk_docker
        Int? preemptible
    }

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx2g" Concordance \
        ~{"-L " + intervals} \
        ~{"-XL " + masks} \
        -truth ~{truth_vcf} -eval ~{eval_vcf} \
        -tpfp "tpfp.vcf" \
        -ftnfn "ftnfn.vcf" \
        -summary "summary.txt"
    }

    runtime {
        memory: "5 GB"
        bootDiskSizeGb: 12
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible, 2])
    }

    output {
        File tpfp = "tpfp.vcf"
        File tpfp_idx = "tpfp.vcf.idx"
        File ftnfn = "ftnfn.vcf"
        File ftnfn_idx = "ftnfn.vcf.idx"
        File summary = "summary.txt"
    }
}

task MakeTable {
    input {
        File tpfp
        File tpfp_idx
        File ftnfn
        File ftnfn_idx

        File? gatk_override
        String gatk_docker
        Int? preemptible
    }

    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        touch output.table
        for file in ~{tpfp} ~{ftnfn}; do
        gatk --java-options "-Xmx2g" SelectVariants -V $file --restrict-alleles-to BIALLELIC -O biallelic.vcf
        gatk --java-options "-Xmx2g" VariantsToTable -V biallelic.vcf \
        -F CHROM -F POS -F REF -F ALT -F POPAF -F TLOD -F STATUS -F REF_BASES -GF FRS \
        --show-filtered \
        -O tmp.table

        # if it's the first table, copy it to the output; otherwise copy all but the header line
        if [ -f output.table ]; then
        mv tmp.table output.table
        else
        tail -n +2 tmp.table >> output.table
        fi
        done
    }

    runtime {
        memory: "5 GB"
        bootDiskSizeGb: 12
        docker: "${gatk_docker}"
        disks: "local-disk " + 100 + " HDD"
        preemptible: select_first([preemptible, 2])
    }

    output {
        File table = "output.table"
    }
}