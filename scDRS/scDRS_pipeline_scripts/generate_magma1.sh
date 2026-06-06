        #!/bin/bash
        # generate_magma1.sh: Creates a script (magma1.sh) containing MAGMA commands for gene-level association analysis.
        #################################
        output_file="magma1.sh"
        > "$output_file" # Clear or create the output script file

        # Read summary_Nsample.tsv to create a mapping from trait prefix to sample size (N)
        declare -A prefix_to_N
        while IFS=$'\t' read -r file N; do
            if [[ "$file" == *.tsv && "$file" != "file" ]]; then # Skip header
                prefix=${file%.tsv}  # Remove .tsv suffix
                prefix_to_N[$prefix]=$N
            fi
        done < summary_Nsample.tsv

        # Iterate through all .snpp files and generate MAGMA commands
        for snpp_file in *.snpp; do
            if [[ -f "$snpp_file" ]]; then
                prefix=${snpp_file%.snpp}  # Remove .snpp suffix
                if [[ -n "${prefix_to_N[$prefix]}" ]]; then
                    N=${prefix_to_N[$prefix]}
                    # Construct the MAGMA command
                    cmd="magma --bfile g1000_eur_qc --pval \"$snpp_file\" N=$N --gene-annot g1000_eur_qc.genes.annot --out \"$prefix\""
                    echo "$cmd" >> "$output_file"
                else
                    echo "Warning: No N value found for prefix $prefix in summary_Nsample.tsv" >&2
                fi
            fi
        done
        echo "MAGMA commands written to $output_file. Run 'bash $output_file' or 'parallel -j <num_jobs> < $output_file' to execute."

