version 1.0

task convert_allele_table {
	input {
		File allele_table_tsv
		String? output_tsv_path = "allele_data_filtered.tsv"
		String docker_image = "us-central1-docker.pkg.dev/operations-portal-427515/mad4hatter-drugres/python-scripts:latest"
	}

	command <<<
		python3 /opt/python-scripts/convert_allele_table.py \
			--input ~{allele_table_tsv} \
			--output ~{output_tsv_path}
	>>>

	output {
		File converted_allele_table = "~{output_tsv_path}"
	}

	runtime {
		docker: docker_image
	}
}