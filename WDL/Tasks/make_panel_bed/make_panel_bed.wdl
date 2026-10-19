version 1.0

task make_panel_bed {
	input {
		File amplicon_info_tsv
		File reference_fasta
		String? bed_output_path = "Mad4hatter_insert_locs_of_panel.bed"
		String docker_image = "us-central1-docker.pkg.dev/operations-portal-427515/mad4hatter-drugres/python-scripts:latest"
	}

	command <<<
		python3 /opt/python-scripts/make_panel_bed.py \
			--amplicon-info ~{amplicon_info_tsv} \
			--reference ~{reference_fasta} \
			--output ~{bed_output_path}
	>>>

	output {
		File panel_bed = "~{bed_output_path}"
	}

	runtime {
		docker: docker_image
	}
}